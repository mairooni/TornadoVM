#!/usr/bin/env python3
"""
Aggregate the probe results into the tables that go in the report.

Reports medians, not means: cold numbers have a long right tail (page cache, JIT, GPU wake-up) and
a single slow run would drag a mean somewhere no individual run ever was.

The two A/B pairs are printed separately and labelled, because they answer different questions and
only one of them is a controlled experiment:

  MECHANISM (JDK 21)  upstream/develop -> removal branch is a clean ancestor relationship, so the
                      delta is attributable to the JVMCI removal.
  PRODUCT   (JDK 25)  upstream/jdk25 has diverged (730 behind develop, 437 of its own commits), so
                      the delta is "what a JDK 25 user sees", not the cost of the removal.

Usage:  ./03-aggregate.py [results-dir]
"""
import json
import os
import re
import statistics
import sys
from collections import defaultdict

RESULTS = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
    os.environ.get("PERF_WORK", "~/perf-jvmci-work") + "/results")

PHASES = ["cold", "warm-new", "warm-same"]


def load():
    """samples[(label, jdk)][phase] -> [ms, ...]"""
    samples = defaultdict(lambda: defaultdict(list))
    for name in sorted(os.listdir(RESULTS)):
        m = re.match(r"(.+)__on-jdk(\d+)\.tsv$", name)
        if not m:
            continue
        label, jdk = m.group(1), int(m.group(2))
        with open(os.path.join(RESULTS, name)) as fh:
            for line in fh:
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                phase, _idx, nanos = parts[0], parts[1], parts[2]
                try:
                    samples[(label, jdk)][phase].append(int(nanos) / 1e6)
                except ValueError:
                    continue
    return samples


def med(xs):
    return statistics.median(xs) if xs else float("nan")


def fmt(x):
    return "—" if x != x else (f"{x:.2f}" if x < 10 else f"{x:.1f}")


def table(samples):
    print("## Wall clock, profiler OFF (median ms)\n")
    print("| SDK | run JDK | cold | warm-new | warm-same | n(cold) |")
    print("|---|---|---|---|---|---|")
    for (label, jdk) in sorted(samples, key=lambda k: (k[0], k[1])):
        s = samples[(label, jdk)]
        print(f"| {label} | {jdk} | {fmt(med(s['cold']))} | {fmt(med(s['warm-new']))} "
              f"| {fmt(med(s['warm-same']))} | {len(s['cold'])} |")
    print()


def ab(samples, base, new, jdk, kind, caveat):
    b, n = samples.get((base, jdk)), samples.get((new, jdk))
    if not b or not n:
        print(f"_({kind} pair on JDK {jdk} incomplete — missing "
              f"{base if not b else new}; skipped.)_\n")
        return
    print(f"## {kind} A/B — JDK {jdk}: `{base}` (JVMCI) vs `{new}` (reflection)\n")
    print(f"{caveat}\n")
    print("| phase | JVMCI | reflection | delta | delta % |")
    print("|---|---|---|---|---|")
    for ph in PHASES:
        mb, mn = med(b[ph]), med(n[ph])
        if mb != mb or mn != mn:
            continue
        d = mn - mb
        pct = (d / mb * 100) if mb else float("nan")
        print(f"| {ph} | {fmt(mb)} | {fmt(mn)} | {d:+.2f} | {pct:+.1f}% |")
    print()


def profiler(results):
    """Phase attribution. TOTAL_DRIVER_COMPILE_TIME is the control: the GPU driver cannot know how
    the source was produced, so if it moves, the measurement is wrong."""
    keys = ["TOTAL_GRAAL_COMPILE_TIME", "TOTAL_CODE_GENERATION_TIME",
            "TOTAL_DRIVER_COMPILE_TIME", "TOTAL_BYTE_CODE_GENERATION", "TOTAL_TASK_GRAPH_TIME"]
    rows = {}
    for name in sorted(os.listdir(results)):
        m = re.match(r"(.+)__on-jdk(\d+)\.profiler\.json$", name)
        if not m:
            continue
        totals = defaultdict(float)
        try:
            with open(os.path.join(results, name)) as fh:
                text = fh.read()
            # The dump is one JSON object per execution, appended. Sum each timer over all of them
            # rather than trusting a single record to be representative.
            for obj in re.finditer(r'"(' + "|".join(keys) + r')"\s*:\s*"?(\d+)"?', text):
                totals[obj.group(1)] += int(obj.group(2))
        except (OSError, ValueError):
            continue
        if totals:
            rows[(m.group(1), int(m.group(2)))] = totals
    if not rows:
        return
    print("## Phase attribution, profiler ON (summed ns — attribution only, not wall clock)\n")
    print("> A profiler-on run is a different workload; use this for *where*, never for *how much*.")
    print("> `TOTAL_DRIVER_COMPILE_TIME` is the control — it must not move between JVMCI and reflection.\n")
    print("| SDK | run JDK | " + " | ".join(k.replace("TOTAL_", "").replace("_TIME", "") for k in keys) + " |")
    print("|---" * (len(keys) + 2) + "|")
    for (label, jdk), t in sorted(rows.items()):
        print(f"| {label} | {jdk} | " + " | ".join(f"{t.get(k, 0)/1e6:.1f}" for k in keys) + " |")
    print()


def main():
    if not os.path.isdir(RESULTS):
        sys.exit(f"no results directory: {RESULTS}")
    samples = load()
    if not samples:
        sys.exit(f"no .tsv results in {RESULTS}")
    table(samples)
    ab(samples, "jvmci-jdk21", "reflect-jdk21", 21, "MECHANISM",
       "`upstream/develop` is a clean ancestor of the removal branch, so this delta **is** the cost "
       "of the JVMCI removal. This is the number to defend.")
    ab(samples, "jvmci-jdk25", "reflect-jdk22plus", 25, "PRODUCT",
       "`upstream/jdk25` has diverged from the removal branch by ~730 commits plus 437 of its own, "
       "so this delta includes unrelated change. Read it as *what a JDK 25 user sees*, **not** as "
       "the cost of the removal.")
    print("## Cross-JDK variation of the reflection path (no baseline needed)\n")
    print("| run JDK | cold | warm-new | warm-same |")
    print("|---|---|---|---|")
    for (label, jdk) in sorted(samples):
        if label != "reflect-jdk22plus":
            continue
        s = samples[(label, jdk)]
        print(f"| {jdk} | {fmt(med(s['cold']))} | {fmt(med(s['warm-new']))} | {fmt(med(s['warm-same']))} |")
    print()
    profiler(RESULTS)


if __name__ == "__main__":
    main()

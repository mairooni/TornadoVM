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

The report is written to a FILE as well as the console, because a number nobody can find again is
not evidence. Default: $PERF_WORK/report-<timestamp>.md, plus a report-latest.md symlink.

Usage:  ./03-aggregate.py [results-dir] [-o report.md]
"""
import json
import os
import re
import statistics
import sys
from collections import defaultdict

import datetime

_args = [a for a in sys.argv[1:]]
_out = None
if "-o" in _args:
    i = _args.index("-o")
    _out = _args[i + 1]
    del _args[i:i + 2]

PERF_WORK = os.path.expanduser(os.environ.get("PERF_WORK", "~/perf-jvmci-work"))
RESULTS = _args[0] if _args else os.path.join(PERF_WORK, "results")
REPORT = _out or os.path.join(
    PERF_WORK, "report-" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + ".md")

PHASES = ["cold", "warm-new", "warm-same"]
TIERS = ["S", "M", "L", "XL"]
TIER_DESC = {"S": "elementwise", "M": "math-heavy", "L": "DFT", "XL": "n-body"}


def load():
    """samples[(label, jdk)][phase] -> [ms]; tiered[(label, jdk)][(phase, tier)] -> [ms]"""
    samples = defaultdict(lambda: defaultdict(list))
    tiered = defaultdict(lambda: defaultdict(list))
    for name in sorted(os.listdir(RESULTS)):
        m = re.match(r"(.+)__on-jdk(\d+)\.tsv$", name)
        if not m:
            continue
        label, jdk = m.group(1), int(m.group(2))
        with open(os.path.join(RESULTS, name)) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                phase, tier, _idx, nanos = parts[0], parts[1], parts[2], parts[3]
                try:
                    ms = int(nanos) / 1e6
                except ValueError:
                    continue
                samples[(label, jdk)][phase].append(ms)
                tiered[(label, jdk)][(phase, tier)].append(ms)
    return samples, tiered


_buf = []


def print(*a, **kw):  # noqa: A001 - deliberately shadowed so every table lands in the report too
    import builtins
    line = " ".join(str(x) for x in a)
    _buf.append(line)
    builtins.print(*a, **kw)


def environment():
    """Reproduce the run's context at the top of the report: which machine, GPU, driver and, above
    all, which commits each SDK was built from."""
    env = os.path.join(RESULTS, "environment.txt")
    if not os.path.isfile(env):
        return
    print("## Environment\n")
    print("```")
    with open(env) as fh:
        print(fh.read().rstrip())
    print("```\n")


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
    """Compile-phase medians across repetitions, with the control used as the significance test.

    TOTAL_DRIVER_COMPILE_TIME is source -> device binary, done by the GPU driver, which cannot know
    whether the source came from JVMCI or reflection. It therefore MUST NOT move. However much it
    moves anyway is this machine's noise floor, and any compile delta smaller than that is not a
    result -- it is the same noise showing up in a different column.
    """
    keys = ["TOTAL_GRAAL_COMPILE_TIME", "TOTAL_CODE_GENERATION_TIME", "TOTAL_DRIVER_COMPILE_TIME"]
    runs = defaultdict(lambda: defaultdict(list))   # (label, jdk) -> key -> [ms per run]
    for name in sorted(os.listdir(results)):
        m = re.match(r"(.+)__on-jdk(\d+)\.profiler(?:\.\d+)?\.json$", name)
        if not m:
            continue
        try:
            text = open(os.path.join(results, name)).read()
        except OSError:
            continue
        for k in keys:
            tot = sum(int(x.group(1)) for x in
                      re.finditer(r'"' + k + r'"\s*:\s*"?(\d+)"?', text))
            if tot:
                runs[(m.group(1), int(m.group(2)))][k].append(tot / 1e6)
    if not runs:
        return
    print("## Compile phase, profiler ON (median ms per run; 20 kernels per run)\n")
    print("| SDK | run JDK | n | " + " | ".join(k.replace("TOTAL_", "").replace("_TIME", "") for k in keys) + " |")
    print("|---" * (len(keys) + 3) + "|")
    for (label, jdk), d in sorted(runs.items()):
        n = max((len(v) for v in d.values()), default=0)
        print(f"| {label} | {jdk} | {n} | " + " | ".join(fmt(med(d.get(k, []))) for k in keys) + " |")
    print()

    for base, new_, jdk in (("jvmci-jdk21", "reflect-jdk21", 21),
                            ("jvmci-jdk25", "reflect-jdk22plus", 25)):
        b, n_ = runs.get((base, jdk)), runs.get((new_, jdk))
        if not b or not n_:
            continue
        ctrl_b, ctrl_n = med(b.get(keys[2], [])), med(n_.get(keys[2], []))
        floor = abs(ctrl_n - ctrl_b) / ctrl_b * 100 if ctrl_b else float("nan")
        print(f"### Compile-phase delta vs the noise floor — JDK {jdk}\n")
        print(f"Control (`DRIVER_COMPILE`) moved **{floor:+.1f}%** between the two SDKs. It should be 0%, "
              f"so treat {floor:.1f}% as this machine's noise floor: a compile delta smaller than that "
              f"is not a measurement.\n")
        print("| phase | JVMCI | reflection | delta % | above the floor? |")
        print("|---|---|---|---|---|")
        for k in keys[:2]:
            mb, mn = med(b.get(k, [])), med(n_.get(k, []))
            if mb != mb or mn != mn or not mb:
                continue
            pct = (mn - mb) / mb * 100
            verdict = "**yes**" if abs(pct) > floor else "no — within noise"
            print(f"| {k.replace('TOTAL_', '').replace('_TIME', '')} | {fmt(mb)} | {fmt(mn)} | {pct:+.1f}% | {verdict} |")
        print()


def scaling(tiered, base, new, jdk):
    """Per-kernel compile cost by tier -- the slope that says whether the reflection overhead grows
    with graph size. A flat delta means a fixed per-kernel cost; a rising one means it scales with
    the number of metadata calls, which is what decides the cost for a real application."""
    b, n = tiered.get((base, jdk)), tiered.get((new, jdk))
    if not b or not n:
        return
    print(f"## Per-kernel compile cost by kernel complexity — JDK {jdk} (warm-new, median ms)\n")
    print("| tier | kernel | JVMCI | reflection | delta | delta % | n |")
    print("|---|---|---|---|---|---|---|")
    for tier in TIERS:
        xb, xn = b.get(("warm-new", tier), []), n.get(("warm-new", tier), [])
        if not xb or not xn:
            continue
        mb, mn = med(xb), med(xn)
        d = mn - mb
        pct = (d / mb * 100) if mb else float("nan")
        print(f"| {tier} | {TIER_DESC.get(tier, '')} | {fmt(mb)} | {fmt(mn)} | {d:+.2f} | {pct:+.1f}% | {len(xn)} |")
    print()


def main():
    if not os.path.isdir(RESULTS):
        sys.exit(f"no results directory: {RESULTS}")
    samples, tiered = load()
    if not samples:
        sys.exit(f"no .tsv results in {RESULTS}")
    print(f"# JVMCI-removal compile cost — {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    environment()
    table(samples)
    ab(samples, "jvmci-jdk21", "reflect-jdk21", 21, "MECHANISM",
       "`upstream/develop` is a clean ancestor of the removal branch, so this delta **is** the cost "
       "of the JVMCI removal. This is the number to defend.")
    ab(samples, "jvmci-jdk25", "reflect-jdk22plus", 25, "PRODUCT",
       "`upstream/jdk25` has diverged from the removal branch by ~730 commits plus 437 of its own, "
       "so this delta includes unrelated change. Read it as *what a JDK 25 user sees*, **not** as "
       "the cost of the removal.")
    scaling(tiered, "jvmci-jdk21", "reflect-jdk21", 21)
    scaling(tiered, "jvmci-jdk25", "reflect-jdk22plus", 25)
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

    os.makedirs(os.path.dirname(REPORT) or ".", exist_ok=True)
    with open(REPORT, "w") as fh:
        fh.write("\n".join(_buf) + "\n")
    latest = os.path.join(os.path.dirname(REPORT) or ".", "report-latest.md")
    try:
        if os.path.islink(latest) or os.path.exists(latest):
            os.remove(latest)
        os.symlink(os.path.basename(REPORT), latest)
    except OSError:
        pass
    import builtins
    builtins.print(f"\n[report written to {REPORT}]", file=sys.stderr)


if __name__ == "__main__":
    main()

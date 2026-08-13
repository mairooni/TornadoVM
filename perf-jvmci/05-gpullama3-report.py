#!/usr/bin/env python3
"""
Aggregate the GPULlama3 runs into the report tables.

Two metrics, and it matters which is which:

  tokens/second   steady-state decode. Nothing about metadata access can reach a kernel that is
                  already compiled and resident on the device, so this is the CONTROL: it must not
                  move. If it does, the measurement is contaminated (thermal, clocks, other load)
                  and the startup numbers next to it cannot be trusted either.
  startup         wall clock minus decode time: JVM start, model load, and every kernel compile.
                  This is the half where a per-kernel compile regression accumulates, and also
                  where the removal branch's smaller module path wins time back.

Usage:  ./05-gpullama3-report.py [results-dir] [-o report.md]
"""
import os
import re
import statistics
import sys
import datetime

_args = sys.argv[1:]
_out = None
if "-o" in _args:
    i = _args.index("-o")
    _out = _args[i + 1]
    del _args[i:i + 2]

PERF_WORK = os.path.expanduser(os.environ.get("PERF_WORK", "~/perf-jvmci-work"))
RESULTS = _args[0] if _args else os.path.join(PERF_WORK, "gpullama3-results")
REPORT = _out or os.path.join(
    PERF_WORK, "gpullama3-report-" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S") + ".md")

SDKS = ["jvmci-jdk21", "reflect-jdk21"]
_buf = []


def emit(line=""):
    _buf.append(line)
    print(line)


def num(s):
    """The launcher prints 'achieved tok/s: 156.10.' -- the sentence's full stop ends up in the
    captured field, so strip any trailing punctuation before parsing."""
    return float(s.rstrip("."))


def load():
    """data[model][sdk] -> dict of metric -> [values]"""
    data = {}
    for name in sorted(os.listdir(RESULTS)):
        m = re.match(r"(.+)__(.+)\.tsv$", name)
        if not m:
            continue
        model, sdk = m.group(1), m.group(2)
        rows = data.setdefault(model, {}).setdefault(sdk, {"toks": [], "decode": [], "wall": [], "startup": []})
        for line in open(os.path.join(RESULTS, name)):
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            try:
                rows["toks"].append(num(p[1]))
                rows["decode"].append(num(p[3]))
                rows["wall"].append(num(p[4]))
                rows["startup"].append(num(p[5]))
            except ValueError:
                continue
    return data


def med(xs):
    return statistics.median(xs) if xs else float("nan")


def f(x, nd=2):
    return "—" if x != x else f"{x:.{nd}f}"


def main():
    if not os.path.isdir(RESULTS):
        sys.exit(f"no results directory: {RESULTS}")
    data = load()
    if not data:
        sys.exit(f"no .tsv results in {RESULTS}")

    emit(f"# GPULlama3 — JVMCI vs reflection, JDK 21 — {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}\n")

    cfg = os.path.join(RESULTS, "config.txt")
    if os.path.isfile(cfg):
        emit("## Configuration\n")
        emit("```")
        emit(open(cfg).read().rstrip())
        emit("```\n")

    emit("## Tokens/second — the control\n")
    emit("Steady-state decode. Metadata access cannot reach an already-compiled, resident kernel, "
         "so this **must not move**. A difference here means the measurement is contaminated.\n")
    emit("| model | JVMCI tok/s | reflection tok/s | delta | delta % |")
    emit("|---|---|---|---|---|")
    for model in sorted(data):
        a, b = data[model].get(SDKS[0]), data[model].get(SDKS[1])
        if not a or not b:
            continue
        ma, mb = med(a["toks"]), med(b["toks"])
        emit(f"| {model} | {f(ma)} | {f(mb)} | {mb - ma:+.2f} | {(mb - ma) / ma * 100:+.1f}% |")
    emit()

    emit("## Startup — JVM start + model load + every kernel compile\n")
    emit("Wall clock minus decode. This is where a per-kernel compile cost accumulates, and where "
         "the removal branch's smaller module path wins time back.\n")
    emit("| model | JVMCI s | reflection s | delta | delta % |")
    emit("|---|---|---|---|---|")
    for model in sorted(data):
        a, b = data[model].get(SDKS[0]), data[model].get(SDKS[1])
        if not a or not b:
            continue
        ma, mb = med(a["startup"]), med(b["startup"])
        emit(f"| {model} | {f(ma)} | {f(mb)} | {mb - ma:+.2f} | {(mb - ma) / ma * 100:+.1f}% |")
    emit()

    emit("## End-to-end wall clock\n")
    emit("| model | JVMCI s | reflection s | delta | delta % |")
    emit("|---|---|---|---|---|")
    for model in sorted(data):
        a, b = data[model].get(SDKS[0]), data[model].get(SDKS[1])
        if not a or not b:
            continue
        ma, mb = med(a["wall"]), med(b["wall"])
        emit(f"| {model} | {f(ma)} | {f(mb)} | {mb - ma:+.2f} | {(mb - ma) / ma * 100:+.1f}% |")
    emit()

    emit("## Raw medians and spread\n")
    emit("| model | SDK | n | tok/s (min–max) | decode s | startup s | wall s |")
    emit("|---|---|---|---|---|---|---|")
    for model in sorted(data):
        for sdk in SDKS:
            d = data[model].get(sdk)
            if not d:
                continue
            t = sorted(d["toks"])
            emit(f"| {model} | {sdk} | {len(t)} | {f(med(t))} ({f(t[0])}–{f(t[-1])}) "
                 f"| {f(med(d['decode']))} | {f(med(d['startup']))} | {f(med(d['wall']))} |")
    emit()

    os.makedirs(os.path.dirname(REPORT) or ".", exist_ok=True)
    with open(REPORT, "w") as fh:
        fh.write("\n".join(_buf) + "\n")
    print(f"\n[report written to {REPORT}]", file=sys.stderr)


if __name__ == "__main__":
    main()

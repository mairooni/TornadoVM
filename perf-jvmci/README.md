# Measuring the compile / code-gen cost of the JVMCI removal

The JVMCI removal replaced HotSpot's native metadata access with Java reflection plus classfile
parsing. This harness prices that: **how much slower is kernel compilation, cold and warm, on each
JDK.**

## The one thing to get right

There are two comparisons and they answer different questions. Never merge them into one table.

| | comparison | what it means |
|---|---|---|
| **MECHANISM** | `upstream/develop` vs removal branch, **both on JDK 21** | `develop` is a clean ancestor of the removal branch (0 behind / 167 ahead), so the delta **is** the JVMCI removal. **This is the number to defend.** |
| **PRODUCT** | `upstream/jdk25` vs removal branch, **both on JDK 25** | `upstream/jdk25` is the maintained JDK 25 line but has diverged (~730 commits behind `develop`, 437 of its own). The delta includes unrelated change — quote it as *what a JDK 25 user sees*, never as the cost of the removal. |

There is deliberately **no** baseline for JDK 22/23/24 (those branches are 1,800–2,700 commits
behind — a comparison would attribute years of unrelated change to reflection), none for JDK 26 (no
such branch), and none possible for JDK 27 (JVMCI is gone from the JDK). For those, the harness
reports **cross-JDK variation of the new path**, which needs no baseline.

## What is measured

Three states, because "cold vs warm" hides the one that matters:

| state | how | what it prices |
|---|---|---|
| **cold** | first kernel in a fresh JVM | classfile read + parse, cache fill, cold JIT |
| **warm-new** | kernels 2..8, each a **new holder class**, warm JIT | the per-new-kernel cost — **what a large application pays repeatedly** |
| **warm-same** | re-executing an already-compiled graph | everything cached; should be ~free. If it isn't, a cache is missing. |

The workload is 20 kernels in four complexity tiers, taken from `tornado-benchmarks/ComputeKernels`
so the graph shapes are ones real applications compile:

| tier | kernel | why |
|---|---|---|
| S | elementwise add | the floor: smallest possible graph |
| M | math-heavy elementwise | loop plus `TornadoMath` calls |
| L | DFT | nested loop, trig, many locals |
| XL | n-body | triple-nested, array allocation, heaviest |

Tiers exist so the result is a **slope, not a point**: metadata cost should grow with graph size
while the one-off start-up difference does not, and only several sizes can separate those.

### Which metric to trust

**Use the profiler's compile-phase timers, not the wall clock, for per-kernel cost.** A
first-execution wall clock bundles compilation with GPU execution, and on the reference box the
execution variance (IQRs of 3–5 ms) completely swamped a sub-millisecond compile delta — all four
tiers came out statistically indistinguishable that way. `TOTAL_GRAAL_COMPILE_TIME` excludes
execution and resolves far better.

Wall clock remains the right metric for **cold**, where start-up is the thing being measured.

### The control decides what counts as a result

`TOTAL_DRIVER_COMPILE_TIME` is source → device binary, performed by the GPU driver, which cannot
know whether the source came from JVMCI or reflection. **It must not move.** However much it moves
anyway is the machine's noise floor, and `03-aggregate.py` prints it as such: any compile delta
smaller than the control's own drift is not a measurement. On a noisy developer laptop that floor
was ~5%, which is precisely why the real run belongs on a quiet machine.

## Prerequisites

- JDKs 21 and 25 (mandatory), plus any of 22/23/24/26/27 you want in the variation sweep.
  Discovered automatically under `/opt/jenkins/jdks` and `~/.sdkman/candidates/java`; each
  candidate is **verified** with `java -version`, so a misleading directory name cannot select the
  wrong JDK.
- Whatever the normal TornadoVM build needs (`CMAKE_ROOT`, `CUDA_PATH`, …). On cyclone, source the
  same environment the CI job uses.
- ~4 GB free for four SDKs plus worktrees.

## Running it

```bash
export PERF_WORK=$HOME/perf-jvmci-work      # SDKs, logs, results (nothing lands in git)
export TORNADO_REPO=/path/to/TornadoVM
export BACKEND=cuda                          # one backend per campaign

./01-build-matrix.sh                         # ~4 sequential builds, 15-60 min
./02-run-probe.sh                            # the measurements
./03-aggregate.py                            # markdown tables
```

Knobs: `REPS` (fresh-JVM repetitions, default 7), `WARM_REPS` (default 20),
`REFLECT_RUN_JDKS` (default `22 23 24 25 26 27`), `JDK_ROOTS`.

To rebuild or re-measure one configuration only:

```bash
./01-build-matrix.sh reflect-jdk22plus
./02-run-probe.sh    reflect-jdk22plus
```

## Why it is built this way

Three decisions exist because of failures we have already hit on this branch:

- **Each build gets its own git worktree.** `graalJars/` is per-checkout state, and the removal
  branch *deletes* the raw `compiler-23.1.0.jar` after relocating it — while `develop` and
  `jdk25` still need that jar for `--upgrade-module-path`. One shared checkout means whichever
  builds second fails, or silently links the wrong Graal.
- **Probes compile against the SDK's own jars, never `~/.m2`.** `develop` and the removal branch's
  jdk21 SDK publish the **same coordinate** (`5.2.1-jdk21-dev`) with different content, so a
  repository lookup would silently pick whichever was installed last.
- **`javac` flags are derived from the SDK**, by reading the class-file version of its
  `TaskGraph.class` (major → `--release`, minor `65535` → `--enable-preview`). Hardcoding these is
  exactly how the sample apps ended up compiling at `release 21 --enable-preview` on JDK 22 and
  failing outright.

Plus one discarded warm-up run per configuration: the first run after a build also faults the SDK
jars in from disk, which measured ~4× higher warm-new times. That is page cache, not compilation.

## Reading the output

`03-aggregate.py` prints wall-clock medians, the two A/B tables (each labelled with its caveat),
the cross-JDK variation table, and phase attribution from the profiler.

**The control row.** `TOTAL_DRIVER_COMPILE_TIME` is source → device binary, done by the GPU driver,
which cannot know how the source was produced. **It must not move between JVMCI and reflection.**
If it does, the measurement is wrong — stop and find out why before believing anything else.

**Profiler numbers are attribution only.** A profiler-on run is a different workload (TornadoVM's
profiler perturbs short task graphs). Use it to say *where* time goes; use profiler-off wall clock
for *how much*.

## Large workload: GPULlama3

The probe prices one kernel at a time. GPULlama3 prices a real application, where the cold cost is
paid across many kernels at once and shows up as **time to first token**.

It needs [PR #146](https://github.com/beehive-lab/GPULlama3.java/pull/146) (same `[21,22)` +
`[22,)` `jdk22plus` consolidation as the ray tracer and kfusion):

```bash
git clone https://github.com/beehive-lab/GPULlama3.java.git && cd GPULlama3.java
git fetch origin pull/146/head:jdk27 && git checkout jdk27

# Build against each SDK in turn. The pom auto-selects its profile from JAVA_HOME, and
# -Dtornadovm.version pins the artifact to the SDK under test.
export TORNADOVM_HOME=$PERF_WORK/sdks/reflect-jdk22plus
export JAVA_HOME=$(ls -d /opt/jenkins/jdks/*25*/ | head -1)
mvn clean package -Dtornadovm.version=5.2.1-jdk22plus-dev

# Time to first token, 5 fresh runs. Compare the same model + same prompt across SDKs.
for i in $(seq 1 5); do
  /usr/bin/time -f "run$i %e s" ./llama-tornado --gpu --model <model>.gguf \
      --prompt "hello" --max-tokens 1 2>&1 | tail -1
done
```

Two things to hold fixed or the numbers mean nothing: the **same model file** (kernel shapes depend
on it) and the **same `--max-tokens`**. Use `--max-tokens 1` to isolate startup + compile, and a
larger value to check that steady-state decode is unaffected — decode is GPU-bound, so it is the
natural **control**: it must not move.

## What invalidates a run

- Mixing backends, or mixing GPUs. One machine, one backend, one GPU = one data point. Say which.
- Comparing across machines: on the CI runners, `thunder` has CUDA 11.5 and a Quadro GP100 while
  `cyclone` has CUDA 12.0 and an RTX 5070 Ti. Not comparable.
- Forgetting `-Dtornado.recover.bailout=False` (the harness always sets it). Without it a failed
  kernel silently runs sequential Java and you are timing the CPU.
- Reusing a worktree between builds — always let `01-build-matrix.sh` recreate it.
- Other load on the box. Check nothing else is using the GPU.

## Reporting

Per the campaign guidance: medians, the configuration (GPU, driver, CUDA, JDK, backend), a control
row that does not move, and the negatives — where the cost is *not* significant, and what stayed
unmeasured. A report that only lists wins reads as unmeasured.

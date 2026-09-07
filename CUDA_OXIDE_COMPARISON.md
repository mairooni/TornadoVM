# cuda-oxide (NVIDIA Labs) vs TornadoVM

An in-depth comparison of two **source-to-GPU compilers**, an assessment of what TornadoVM can adopt
from cuda-oxide, and a critique of the two projects' websites with concrete suggestions for
[tornadovm.org](https://www.tornadovm.org/).

Scope note. This document compares **NVIDIA Labs' cuda-oxide**([nvlabs.github.io/cuda-oxide](https://nvlabs.github.io/cuda-oxide/)) — an experimental**Rust → CUDA PTX** compiler — with TornadoVM. It is *not* about the unrelated Protryon/cuda-oxidecrate (a thin safe-Rust CUDA driver wrapper). The NVIDIA-Labs project is a peer to TornadoVMbecause both compile a high-level language directly to GPU code.

---

## TL;DR

Two compilers that let you write GPU kernels in a mainstream, memory-safe language instead of CUDA
C++:


**cuda-oxide** — a custom `rustc` codegen backend that compiles **idiomatic (SIMT) Rust to PTX**.
  "Not a DSL." CUDA-only, `v0.1.0` alpha, NVIDIA-backed. Its edge is **language-level safety**
  (ownership/borrow-checking, type-checked kernel arguments) and it will track **NVIDIA hardware
  features** (tensor cores, and by lineage fp8/async-copy) quickly.

**TornadoVM** — a Graal **JIT** that compiles **Java bytecode to five backends**
  (CUDA, OpenCL, PTX, Metal, SPIR-V), with automatic parallelization (`@Parallel`), an explicit SIMT
  API (`KernelContext`), task-graph scheduling, tensor-core MMA, a hybrid native-library API, and
  multi-GPU. `v5.0.1`, production, real workloads (GPULlama3 LLM inference).

**Bottom line:** TornadoVM leads decisively on **breadth, maturity, and ecosystem**; cuda-oxide leads
on **compile-time kernel safety** and (as an NVIDIA-Labs project) **NVIDIA feature velocity**. The
three ideas most worth borrowing are: stronger static verification of kernels, a faster path to
native NVIDIA features (fp8 MMA, `cp.async`/TMA), and cleaner async composition ergonomics.

---

## 1. What each is

| | **cuda-oxide (NVIDIA Labs)** | **TornadoVM** |
|---|---|---|
| Kind | Rust → PTX compiler (a `rustc` codegen backend) | Java bytecode → GPU JIT compiler + runtime |
| Source language | Idiomatic, safe(ish) Rust | Idiomatic Java (`@Parallel` + `KernelContext`) |
| Output | PTX (NVIDIA only) | OpenCL C, PTX, CUDA C, Metal, SPIR-V |
| Kernel model | Explicit SIMT (thread/block/warp) | Auto-parallel loops **and** explicit SIMT |
| Positioning | "Not a DSL. A custom `rustc` codegen backend." | "Java on GPUs. And more." |
| Maturity | `v0.1.0` alpha | `v5.0.1`, production |
| Backing | NVIDIA Labs | APT Group, University of Manchester (beehive-lab) |

Both reject the "write a separate C++ kernel + glue" model. cuda-oxide keeps you in Rust and its
crate ecosystem; TornadoVM keeps you in Java and Maven with no JNI.

---

## 2. Compilation model


**cuda-oxide: ahead-of-time.** A custom codegen backend for `rustc` lowers Rust MIR to PTX at
  build time. You get static typing, borrow-checking, and monomorphized specialization "for free"
  from the Rust front-end. Trade-off: PTX/NVIDIA-only, and specialization is fixed at compile time.

**TornadoVM: runtime JIT (Graal).** Java bytecode is sketched and compiled to the target backend
  on first execution, per device. Trade-offs: first-run compile latency (mitigated by CUDA graphs
  and caching), but runtime specialization to the actual device and **portability across five
  backends** from one source. TornadoVM also has an interpreter path for the task-graph bytecode.

The two occupy the same conceptual slot ("high-level source → GPU"), differing on when/where codegen
happens and how many targets it serves.

---

## 3. Feature matrix

Legend: **Y** yes · **~** partial · **N** no. TornadoVM classes are under
`tornado-drivers/cuda/src/main/java/uk/ac/manchester/tornado/drivers/cuda/`.

| Capability | cuda-oxide | TornadoVM | Notes |
|---|:--:|:--:|---|
| Source language | Rust | Java | Both memory-managed, mainstream |
| GPU backends | CUDA/PTX only | 5 (CUDA, OpenCL, PTX, Metal, SPIR-V) | TornadoVM: write once, run on NVIDIA/AMD/Intel/Apple |
| Automatic parallelization | N | Y | TornadoVM `@Parallel`; cuda-oxide is explicit-SIMT only |
| Explicit SIMT (shared mem, barriers, warp shuffle) | Y | Y | TornadoVM `KernelContext` (`allocate*LocalArray`, `localBarrier`, `simdShuffleDown`) |
| Tensor cores / MMA | Y | Y | TornadoVM `KernelContext` MMA (`mma.sync`) + Metal simdgroup |
| FP8 / low precision | ~ (NVIDIA lineage) | ~ | TornadoVM FP8 storage + cuBLASLt FP8 GEMM landed (PR #928); software decode is compute-bound |
| Async op graphs | Y (`DeviceOperation` graphs) | Y (`TaskGraph` + CUDA graphs) | TornadoVM CUDA graphs via `CUDADeviceContextInterface` |
| Streams | Y | Y | `CUDACommandQueue`, `CUDAStreamType` |
| Events | Y | Y | `CUDAEvent`, `CUDAEventPool` |
| Unified memory | ~ | Y | TornadoVM UM in `CUDADevice` (`withCudaUM`) |
| Stream-ordered memory pools (`cuMemAllocAsync`) | ~ | N | Neither strong; a real TornadoVM gap (grow-only buffers) |
| Multi-GPU | ~ | Y | TornadoVM NCCL collectives (PR #920) |
| NVRTC / runtime kernel compile | Y | Y | TornadoVM `CUDACodeCache` (NVRTC path + diagnostics) |
| Native-library interop | N | Y | TornadoVM hybrid API: cuBLAS/cuBLASLt/cuDNN/cuFFT/cuSPARSE/cuTENSOR as task-graph library tasks |
| Compile-time kernel type/memory safety | **Y (strong)** | ~ (weak, runtime) | cuda-oxide's clearest edge — see §4 |
| Maturity | alpha | production (`v5.0.1`) | |
| Ecosystem | Rust crates | JVM / Maven | |

---

## 4. Safety model contrast — cuda-oxide's clearest edge

cuda-oxide inherits Rust's **ownership and borrow-checking** and adds **type-checked kernel
arguments**, so a large class of host:left_right_arrow:device mistakes (mismatched argument types, misuse of buffers)
is caught **at compile time**. The GPU-side story is "safe(ish)" — it does not make arbitrary kernel
code memory-safe, but the surface it exposes is far more guarded than raw CUDA C++.

TornadoVM's safety is **CPU-side and largely runtime**: the JVM plus Panama `MemorySegment` off-heap
buffers give managed host memory, but kernel correctness surfaces at **codegen / NVRTC / runtime**,
not at compile time. Two concrete CUDA-backend codegen mis-lowerings found while bringing up FP8 this
cycle illustrate the gap — they compiled to *wrong* PTX rather than being rejected:


Early `return`s in a helper inlined into a `@Parallel` loop mis-lowered (guard captured only the
  first statement) — filed as **beehive-lab/TornadoVM#931**.

Two loops in separate `if/else` branches merged incorrectly in the phi lowering — filed as
  **beehive-lab/TornadoVM#930**.

A Rust-style front-end would have far less room for this class of silent mis-compilation. This is the
single most valuable idea to borrow (see §5.1).

---

## 5. What TornadoVM can adopt / improve — ranked

### 5.1 Stronger compile-time kernel verification — HIGH impact / LARGE effort
Add static checks over `KernelContext` kernels (argument types, obvious bounds/aliasing, local-memory
sizing) so the class of bugs that today only fail at PTX/NVRTC time fails **at build time**, closer
to cuda-oxide's Rust guarantees. Directly motivated by #930/#931. Filed as an enhancement issue.

### 5.2 NVIDIA hardware-feature velocity — HIGH / LARGE — largely CLOSED by #940
FP8 *storage* + cuBLASLt FP8 GEMM landed (#928); the follow-up (#940) closed most of the
in-kernel gap on the CUDA backend:

- **Native FP8 conversion**: `FP8.e4m3ToFloat`/`e5m2ToFloat` now compile to
  `__nv_cvt_fp8_to_halfraw` (cuda_fp8.h) inside CUDA kernels — the dequant benchmark went from
  0.68 ms/iter (software decode, compute-bound) to 0.22 ms/iter (memory-bound, 2.3x faster than
  the fp16 path). Encode stays software on all backends (rounding-mode parity).
- **Native FP8 MMA**: `mma.sync.aligned.m16n8k32` with e4m3/e5m2 operands and f32 accumulator
  (`KernelContext.mmaFP8E4M3/mmaFP8E5M2`), gated at sm_89 by the tensor-core compiler phase.
- **cp.async**: `KernelContext.asyncCopyToLocal` + commit/wait emit
  `cp.async.ca.shared.global` for the MMA GEMM tile loads — 1.30x on a 2048^3 fp16 GEMM,
  1.36x on the Llama-3.2-1B FFN gate-up shape.
- **BF16 as a first-class MMA input**: `mmaBF16` (m16n8k16 `.bf16.bf16.f32`, sm_80+) plus a
  kernel-safe `BFloat16` codec with native `__bfloat162float` decode on CUDA. Throughput
  parity with fp16 MMA (6161 vs 5995 GFLOP/s on a 2048^3 GEMM) at ~7x the fp16 rounding
  error — a range/robustness tradeoff, not a speed win.

Remaining gap: **TMA + warp-specialization (Hopper sm_90+)** — deliberately not implemented,
since no Hopper-class hardware was available to validate against (an unvalidated Hopper path
is worse than an absent one). Still where a NVIDIA-Labs Rust→PTX compiler would move fastest.

### 5.3 Async / stream / graph ergonomics — MED
cuda-oxide's `DeviceOperation` async graphs are a clean composition model. TornadoVM already has the
machinery (`TaskGraph`, CUDA graphs, streams, events) but could tighten the *author-facing* async API
so overlap/dependencies are expressed as ergonomically.

### 5.4 Idiomatic-source positioning — LOW
cuda-oxide's identity is razor-sharp: "Not a DSL. A custom `rustc` codegen backend." TornadoVM is
"pure Java" but leans on `@Parallel`/`KernelContext`, which read as a light DSL. Worth sharpening the
message (and, over time, reducing DSL-feel).

---

## 6. Where TornadoVM already leads

The features cuda-oxide lacks are exactly TornadoVM's differentiators, and they validate its design:


**Multi-backend** from one source (CUDA/OpenCL/PTX/Metal/SPIR-V) — cuda-oxide is NVIDIA-only.
**Automatic parallelization** (`@Parallel`) in addition to explicit SIMT.
**Production maturity** (`v5.0.1`) vs `v0.1.0` alpha.
**Native-library hybrid API** — cuBLAS/cuBLASLt/cuDNN/cuFFT/cuSPARSE/cuTENSOR as first-class
  task-graph library tasks; cuda-oxide has no native-library interop story.

**Multi-GPU** (NCCL collectives) and a **real end-to-end workload** (GPULlama3 LLM inference).


A source-to-PTX compiler is a strong *floor*; TornadoVM's runtime, portability, and library ecosystem
are the *ceiling*.

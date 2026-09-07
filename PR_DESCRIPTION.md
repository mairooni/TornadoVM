#### Description

Implements #940 (NVIDIA feature parity for the CUDA backend): native FP8 conversion, FP8
tensor-core MMA, cp.async global→shared copies, and BF16 as a first-class MMA input type.
Five commits (all benchmarks: RTX 4070 Laptop, sm_89):

**1. Native FP8 conversion + FP8 tensor-core MMA (58b86b9b8)**

- `FP8.e4m3ToFloat`/`e5m2ToFloat` now compile to `__nv_cvt_fp8_to_halfraw` + `__half2float`
  inside CUDA kernels, with `cuda_fp8.h` injected under the same source-scan gating as
  `cuda_fp16.h`. Decode only: the hardware encoder rounds ties to even while the software
  codec rounds half away from zero, so encode stays software on every backend, and all other
  backends keep inlining the Java decoder. The interception is also conditional on the
  toolkit (cae742c98): `cuda_fp8.h` only exists since CUDA 11.8, so the plugins consult a
  one-time NVRTC probe and fall back to inlining the Java software decoder when the header
  cannot compile — an old or mismatched toolkit degrades FP8 decode performance, not
  correctness (caught by CI on a pre-11.8 toolkit host).
- FP8 (OCP E4M3/E5M2) is wired as a third MMA input format: `mma.sync.aligned.m16n8k32` with
  an f32 accumulator (`KernelContext.mmaFP8E4M3`/`mmaFP8E5M2`). The A/B tile layout is
  byte-identical to the int8 m16n8k32 path, so the FP8 loads reuse the int8 load nodes; only
  the compute op carries the element type (`MMAComputeStmt` no longer infers it from the tile
  shape, which m16n8k32 no longer identifies uniquely).
- Gating: sm_89 (Ada/Hopper — FP8 tensor cores don't exist at the sm_80 floor), enforced
  through the same compiler phase (`TornadoDeviceMMANotSupported`) that gates fp16/int8 MMA
  at sm_80 — not the dispatch-time probe used by the cuBLASLt library path in #928.
- Benchmark (acceptance criterion from #940): the FP8 dequant benchmark drops from
  0.522 ms/iter (software decode, compute-bound) to **0.222 ms/iter** — now **2.29x faster
  than the fp16 hardware-convert path** in the same run, as a memory-bound kernel should be.
  The in-kernel decode gap flagged in #928 is closed.

**2. cp.async global→shared copies (aadac999f)**

- `KernelContext.asyncCopyToLocal` + `asyncCopyCommit`/`asyncCopyWaitGroup` expose Ampere's
  LDGSTS: one 4-byte slot (two fp16/bf16 halves or four fp8/int8 bytes) is copied straight
  from a global Tornado array into a shared int tile, bypassing the register round-trip of
  the previous load-pack-store loops. `cp.async.commit_group`/`cp.async.wait_group` order the
  copies before the barrier that publishes the tile.
- Route taken: inline PTX asm, consistent with the existing dp4a/ldmatrix/mma statements.
  `cuda::pipeline` was rejected because it needs libcudacxx headers under NVRTC, which the
  toolkit include probing does not guarantee on every install.
- Gating: no new gate. cp.async shares mma.sync's sm_80 floor, so its nodes are gated by the
  existing tensor-core support phase — stated explicitly here rather than silently assumed.
- Benchmark: 2048³ fp16 GEMM: 6082 → 7931 GFLOP/s (**1.30x**) over the synchronous-load MMA
  kernel; Llama-3.2-1B FFN gate-up shape (M=128 N=8192 K=2048): 0.944 → 0.696 ms (**1.36x**).

**3. BF16 as a first-class MMA input type (55ff301c1)**

- New `types.BFloat16` codec (kernel-safe software encode/decode, same style as `FP8.java`);
  on CUDA the decoder lowers to a single `__int_as_float(bits << 16)` bit reinterpretation
  (cae742c98) — bf16 is the high half of the f32 pattern, and `__int_as_float` is a core
  CUDA builtin, so the decode needs no header and no minimum toolkit version (the first
  cut used `cuda_bf16.h` helpers that older 11.x header revisions lack, which CI caught).
  Encode stays software everywhere (ties-to-even vs half-away-from-zero, as with FP8).
- `KernelContext.mmaBF16` emits `mma.sync.aligned.m16n8k16.row.col.f32.bf16.bf16.f32`.
  Confirmed against the PTX ISA: bf16 shares fp16's m16n8k16 tile shape and per-lane fragment
  layout, so the existing `mmaLoadA`/`mmaLoadB` fragment loads work unchanged on raw bf16 bit
  pairs; only the mma.sync element type differs. Same sm_80 gate as fp16/int8 MMA.
- HalfFloat-replacement audit: every backend (including CUDA) runs a
  `TornadoHalfFloatReplacement` phase that keys off `HalfFloat`/`HalfFloatArray` node types.
  BF16 values deliberately travel as raw bits in `ShortArray`, so the phase never sees them
  on any backend; and the `HalfFloat[]` MMA fragment handles are the same opaque register
  tuples the existing fp16 MMA path already runs through that phase safely — they are never
  element-accessed, so no fp16 decode is ever applied to bf16 bits. No fp16-only assumption
  is reachable from this path.
- Numerical delta: M=N=K=64, U[−1,1] inputs, same data quantized to both formats and compared
  against the exact f32 GEMM: max |err| = **0.0207 (bf16) vs 0.0029 (fp16), 7.2x** — tracking
  the 3-bit mantissa gap (2³ = 8x expected). Asserted in
  `TestMatrixMultiplicationMMABF16#testBF16VersusFP16Delta`, printed in the test log.
- Benchmark: 2048³ GEMM: **6161 GFLOP/s (bf16) vs 5995 GFLOP/s (fp16)** — throughput parity
  (1.03x, within run-to-run noise), as expected for identical tile shapes and bandwidth.
  BF16 is a range/robustness tradeoff, not a speed win.

**4. Non-CUDA backends reject the new intrinsics (15975fcb0)**

- The OpenCL, PTX, SPIRV and Metal graph-builder plugins register `unimplemented()` stubs for
  all ten new `KernelContext` methods (FP8 loads/computes, `mmaBF16`, the three
  `asyncCopyToLocal` overloads, `asyncCopyCommit`/`asyncCopyWaitGroup`) — the same pattern
  #867 used for the original MMA intrinsics — so a kernel using them fails the sketch with a
  clear "only supported for the CUDA backend" diagnostic instead of inlining the Java
  fallback bodies and silently computing wrong results. PTX gets the stubs too: it implements
  the fp16/int8 MMA ops but not these extensions. (The `BFloat16`/`FP8` codecs themselves
  stay usable on every backend — they're plain kernel-safe Java; only the tensor-core and
  cp.async intrinsics are gated.)

**Skipped: TMA + warp-specialization (Hopper, sm_90+)** — no Hopper-class hardware was
available to validate against, and an unvalidated Hopper path is worse than an absent one.

**Notes on #930/#931 (latent CUDA codegen bugs):** neither bug recurred in the new codegen
(the FP8 convert node, cp.async statements and BF16 convert node are all straight-line
single-statement emissions). The one place the patterns could have appeared is the software
`BFloat16` codec, which must compile on all backends: it is written in the FP8.java
kernel-safe style (single return through a result variable to avoid #931; sequential
power-of-two loops on one accumulator instead of branched loops to avoid #930), with comments
citing both issues at each spot. This is the documented workaround, not a root fix — the
phi-lowering (#930) and inliner (#931) fixes remain open follow-ups, deliberately out of this
PR's scope.

**Reviewer note — NVRTC/toolkit header-mismatch hazard (out of scope, flagging like
#930/#931):** while testing, kernels that `#include <cuda_fp8.h>` failed to compile when
`CUDA_PATH` was unset on a machine with libnvrtc 12.5 loaded but CUDA 13.3 headers in
`/usr/local/cuda/include`: the JNI's include probing falls back to `/usr/local/cuda/include`,
and 13.3's `cuda_fp8.hpp` uses macros (`__NV_SILENCE_DEPRECATION_BEGIN`) that NVRTC 12.5's
builtin headers don't define. The probing could verify header/NVRTC version agreement (e.g.
compare `CUDA_VERSION` in the probed `cuda.h` against `nvrtcVersion()`) before adding an
include dir. Filing as a follow-up rather than fixing here.

#### Problem description

Not a bug fix — this closes the NVIDIA hardware-feature gaps tracked in #940. The concrete
problem it removes: in-kernel FP8 decode was the software arithmetic path from `FP8.java`
(2–10x slower than a hardware convert per the #928 benchmark), MMA GEMM tile loads
round-tripped through registers instead of using Ampere's async copies, and BF16 had no MMA
path at all.

#### Backend/s tested

Mark the backends affected by this PR.

- [x] OpenCL
- [x] PTX
- [x] CUDA
- [x] SPIRV
- [x] Metal

CUDA is the functional target: all new/modified test classes plus `make fast-tests` are
green on it. The other four backends are affected only by the `unimplemented()` stubs;
OpenCL was additionally tested at runtime (kernels invoking `asyncCopyToLocal`/`mmaBF16`
abort compilation with the expected diagnostic), while PTX/SPIRV/Metal are compile-verified
(identical stub code; no hardware/OS available here to execute them).

#### OS tested

Mark the OS where this PR is tested.

- [x] Linux
- [ ] OSx
- [ ] Windows

#### Did you check on FPGAs?

If it is applicable, check your changes on FPGAs.

- [ ] Yes
- [x] No

Not applicable — the new intrinsics are CUDA-only; FPGA (OpenCL) backends only see the new
`unimplemented()` stubs.

#### How to test the new patch?

Unit tests (CUDA backend, sm_80+; the FP8 MMA tests need sm_89+):

```bash
tornado-test -V uk.ac.manchester.tornado.unittests.arrays.TestFP8
tornado-test -V uk.ac.manchester.tornado.unittests.arrays.TestBFloat16
tornado-test -V uk.ac.manchester.tornado.unittests.kernelcontext.matrices.TestMatrixMultiplicationMMAFP8
tornado-test -V uk.ac.manchester.tornado.unittests.kernelcontext.matrices.TestMatrixMultiplicationMMACpAsync
tornado-test -V uk.ac.manchester.tornado.unittests.kernelcontext.matrices.TestMatrixMultiplicationMMABF16
tornado-test -V uk.ac.manchester.tornado.unittests.kernelcontext.matrices.TestMatrixMultiplicationMMA
```

All test classes are also registered in `tornado-test`'s `__TEST_THE_WORLD__`, so a plain
`make tests` / `make fast-tests` covers them.

Benchmarks (reproduce the numbers above):

```bash
# FP8 dequant: hardware vs software decode vs fp16
tornado -m tornado.examples/uk.ac.manchester.tornado.examples.arrays.FP8Benchmark

# GEMM variants: baseline / MMA / MMA+cp.async / MMA bf16 / swizzled
tornado -m tornado.examples/uk.ac.manchester.tornado.examples.compute.MatrixMultiplicationMMA
```

To see the unsupported-backend behaviour, run any kernel using the new intrinsics on a
non-CUDA device (e.g. force an OpenCL device with `-D<taskgraph>.<task>.device=X:Y`): the
sketch aborts with `unimplemented: ... only supported for the CUDA backend.`

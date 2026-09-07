# Merge decision log: upstream/develop → feat/remove-ptx-backend (PR #872)

Merge base: `b4e1e1ca3`. Develop tip: `604e3c598` (361 commits).
My branch: single commit `f743780c2` (removes the PTX backend; adds no non-PTX
functionality — so for every conflict the correct result is "develop's version
minus the PTX parts").

Policy: PTX removal authoritative; keep all non-PTX develop changes; strip PTX
from mixed files; keep PTX-only files deleted; flag genuine ambiguity.

Merge left in progress (NOT committed). Review with `git diff --cached`, `git diff`.

---

## PTX-only from develop (dropped)

**A. PTX files my branch deleted that develop modified** (32, resolved via `git rm`):
all under `tornado-drivers/ptx/` and `tornado-drivers/ptx-jni/` (PTXContext,
PTXStream, PTXDevice, PTXAssembler, PTXKind, mm/*, ptx-jni cpp/h, poms, etc.).
Policy: file deleted on my branch + modified on develop → keep deleted.

**B. NEW PTX-only files develop added after my base** (15, resolved via `git rm -f`):
these auto-merged in because my deletion set predated them. All under
`tornado-drivers/ptx/`:
  - `PTXEventRegistry.java`, `PTXExecutionStreamSet.java`, `PTXStreamPool.java`,
    `PTXStreamType.java`
  - `graal/nodes/`: `MMAComputeNode`, `MMAFragmentNode`, `MMALoadANode`,
    `MMALoadAInt8Node`, `MMALoadBNode`, `MMALoadBInt8Node`, `MMALoadBSwizzledNode`,
    `MMAStoreNode`, `MMAStoreBSwizzledNode`, `PTXConvertHalfBitsToIntNode`
  - `graal/phases/PTXTensorCoreSupportPhase.java`
  Note: the PTX module was already fully unwired (no pom/settings/service/enum
  references), so these were orphan sources; removed to honor the deletion.

## Non-PTX from develop (accepted — deletions honored)

- `docs/source/cuda-backend.rst`, `docs/source/spirv-backend.rst` — develop DELETED
  both in `fc34766ea [docs] Removed duplicated or old files` (also profiler.rst).
  My branch had only PTX-stripped them. Honored develop's deletion (`git rm`).
- `Jenkinsfile` — develop DELETED it in `cfc5c5698 [cleanup] remove deprecated
  files`. My branch had only swapped ptx→cuda. Honored develop's deletion.
- `docs/source/profiler.rst` → develop renamed/replaced it with `dev-tools.rst`;
  handled under dev-tools.rst below.

## Mixed files (surgically resolved — kept develop's non-PTX, stripped PTX)

### Build / CI / release
- `.github/workflows/build-test-jdk21.yml`, `build-test-jdk25.yml` — develop added a
  new ray-tracer regression step to `build-opencl` AND kept a `build-ptx` job. Kept
  the ray-tracer step (belongs to build-opencl); dropped the entire `build-ptx` job.
  (The two PTX unit-test steps auto-merged out cleanly.)
- `.github/workflows/build-test-with-graalvm-and-polyglot.yml` — develop restructured
  3 jobs into one `build-test-graalvm-polyglot` (non-PTX) and left `BACKEND=ptx,opencl,
  spirv`. Took develop's single-job structure; set `BACKEND=cuda,opencl,spirv`. The two
  PTX test steps auto-merged out.
- `.github/workflows/publish-archives-to-sdkman-jdk21.yml`, `-jdk25.yml` — develop added
  both a `ptx` and a NEW `cuda` SDKMAN matrix entry. Kept cuda, dropped ptx.
- `.github/workflows/build-release-sdks.yml` — job display names listed
  `opencl+ptx+spirv+cuda`; dropped `ptx` (actual builds come from build-release-sdks.py).
- `.github/workflows/1-prepare-release-jdk21.yml`, `-jdk25.yml` — sed regex alternation
  `(opencl|cuda|ptx|spirv|metal|full)` → dropped `ptx|`.
- `.github/pull_request_template.md` — develop added a `- [ ] PTX` backend checkbox; dropped it.
- `Makefile` — MIXED: develop added `cuda-backend` AND `metal-backend` to the clean
  target while keeping `ptx-backend`. Result: `opencl-backend,spirv-backend,cuda-backend,
  metal-backend` (kept metal, dropped ptx). Comment line: took HEAD's `{opencl,cuda,spirv,metal}`.
- `bin/compile` — took HEAD: "full" = `{opencl, cuda, spirv}` (develop had `{opencl,ptx,spirv,cuda}`).
- `scripts/build-release-sdks.py` (new on develop) — removed the `ptx` build variant from
  linux+windows BUILDS; changed "full" combo and label logic from `{opencl,ptx,spirv,cuda}`
  to `{opencl,spirv,cuda}` (matches bin/compile); updated docstring/comments.
- `scripts/HOW_TO_BUILD_SDK.md` (new on develop) — dropped `ptx` rows; set "full" =
  `opencl,spirv,cuda`. NOTE: the develop table omitted `cuda` rows entirely (stale even on
  develop); I added `cuda` rows so the doc matches build-release-sdks.py. (Slightly beyond
  pure strip — flag if undesired.)
- `docs/source/developer-guidelines.rst` — stripped ptx from installer option list,
  `--backend` set, and `make BACKEND=opencl,cuda,spirv` line.
- `tornado-assembly/src/bin/tornado.py` — develop added `validate_ptx_backend()` +
  `validate_cuda_backend()` functions and call sites. Kept the CUDA validation, dropped
  all PTX validation.
- `tornado-assembly/src/bin/tornado` (shell launcher) — removed a vestigial `*cuda*)` case
  calling the undefined `check_ptx_backend` (dead PTX-backend validation branch).

### Docs (took develop's version, reworded PTX-as-backend out; kept legit "NVRTC→PTX" IR mentions)
- `docs/source/index.rst` — "five backends (…NVIDIA PTX…)" → "four backends"; dropped the
  NVIDIA PTX bullet; reworded the CUDA C bullet.
- `docs/source/introduction.rst`, `docs/source/faq.rst` — dropped "PTX" from backend
  enumerations; kept "(compiled to PTX via NVRTC)".
- `docs/source/installation.rst` — took develop's rewritten (SDKMAN-based) installation
  guide; removed `ptx` from the two prebuilt-SDK lists.
- `docs/source/flags.rst` — stripped ptx from `--printKernel` text and the
  `{opencl,…}.priority` row (+ default list); DELETED the entire "PTX Backend Specific
  (CU_JIT Flags)" section (PTX-only, configures `-Dtornado.ptx.compiler.flags`); reworded
  the "CUDA C Backend Specific" intro that had described PTX + CUDA C as two backends.
- `docs/source/dev-tools.rst` (develop's rename of my profiler.rst) — reworded 8
  PTX-as-backend mentions to CUDA (NVML section heading, profiler value descriptions,
  BACKEND example JSON, Nsight bullet).

### API / runtime / drivers (Java)
- `tornado-drivers/cuda/.../plugins/CUDAGraphBuilderPlugins.java` — develop IMPLEMENTED
  swizzled FP16 MMA load/store (`CUDASwizzledLoad/StoreFP16Stride32Node`, new node classes
  present in tree); HEAD had `unimplemented(...)` stubs. Took develop's implementation (new
  non-PTX CUDA feature).
- `tornado-drivers/{opencl,metal,spirv}/.../plugins/*GraphBuilderPlugins.java` — develop
  added many MMA stub methods each throwing `unimplemented("MMA instructions only supported
  for the PTX backend.")`. Kept the stubs; changed the message to "…for the CUDA backend."
  (15 each × 3 files).
- `tornado-runtime/.../sketcher/TornadoDataflowAnalysis.java` — removed a string reference
  to the deleted `drivers.ptx.graal.nodes.MMAStoreNode` from an `||` chain.
- `tornado-runtime/.../common/TornadoXPUDevice.java` — two javadocs said the PTX backend
  overrides `setIntraPlanConcurrency`/`setStagedTransfers`; both are now overridden by the
  CUDA backend → reworded to "CUDA backend".
- `tornado-runtime/.../common/TornadoOptions.java` — "(PTX and CUDA backends)" → "(CUDA backend)".
- `tornado-runtime/.../interpreter/TornadoVMInterpreter.java` — comments referenced deleted
  `PTXEventPool`/`PTXEventRegistry`; reworded to `CUDAEventPool` (exists) / generic
  "event-registry" (no CUDAEventRegistry class exists).
- `tornado-drivers/cuda/.../CUDADeviceContext.java` — dropped a "Mirrors the PTX backend,
  which frees its ring in PTXStream.cuDestroyStream()" sentence (referenced deleted class).
- `tornado-api/.../KernelContext.java` ("Lowered by the PTX backend" → CUDA),
  `TornadoExecutionPlan.java` (×2 "PTX and CUDA backends" → CUDA),
  `enums/MMAShape.java` ("PTX and CUDA backends" → CUDA),
  `tornado-drivers/cuda/.../CUDAConvertFP8ToFloat.java` (backend list OpenCL/PTX/… → OpenCL/…).

### Tests (Java) — several were compile-breaking
- **COMPILE-BREAKING (fixed):** my branch removed the `TornadoVMBackendType.PTX` enum
  constant, but develop-added tests referenced it:
  - `MatrixMultiplicationMMA.java` (example): `isPTXorCUDABackend()` → `isCUDABackend()`,
    dropped `== PTX`, updated message.
  - 10 unittests with `switch(backendType){ case OPENCL, PTX, SPIRV, METAL -> … }` →
    dropped `PTX` from the case label: TestBFloat16, TestFP8, TestHalfFloatInlineWrite,
    TestCuBlas, TestCuBlasLt, TestCuDnn, TestCuFft, TestCusparse, TestCutlass, TestNvtx.
  - 3 MMA tests with `assertNotBackend(TornadoVMBackendType.PTX);` (BF16, CpAsync, FP8) →
    removed the line (they already exclude every non-CUDA backend).
- `TestCUDAStreams.java`, `TestStreamsPerformance.java` — javadoc called these CUDA-stream
  tests "PTX backend" tests → reworded to CUDA.
- `TestSimdgroupMatrix.java`, `TestSimdgroupTiledMatrix.java` — dropped "PTX" from a
  "no equivalent in OpenCL, PTX, SPIR-V or CUDA" comment.
- `HYBRID_API_GUIDE.md` — "UNSUPPORTED on OpenCL/PTX/SPIR-V/Metal" → dropped PTX.

## Intentionally KEPT (legitimate — NOT PTX backend)
These describe the PTX **instruction set / NVRTC intermediate representation** that the
CUDA backend genuinely emits, so they are correct and were left as-is:
- "compiled to PTX via NVRTC" / "emits CUDA C … to PTX" (index, intro, faq, flags,
  cuda-backend narrative, INSTALL_FROM_SOURCE.md:93).
- `KernelContext.java` "PTX equivalent: …", "PTX ISA m16n8k16", "shared memory in PTX".
- `MMAShape.java` "PTX instruction emitted:".
- CUDA driver files: "inline-PTX asm", "PTX ISA 8.4", "ptxas", "the CUDA-C counterpart of
  the PTX dp4a", "mirrors the PTX address arithmetic", `CUDAProgram.java` PTX-ISA gating.
- `tornado-cutlass/README.md`, `HYBRID_API_GUIDE.md:386` "compute_80 PTX" (SASS/PTX target).
- `docs/source/CHANGELOG.rst` — historical release notes mentioning `BACKEND=ptx,opencl`
  (historical record; not rewritten).
- `.claude/skills/tornadovm/SKILL.md`, `references/codegen-map.md` — internal Claude-skill
  tooling docs that already treat PTX as out-of-scope / distinct from the CUDA module.
  (Non-shipped; left as-is. Mildly stale now that PTX is gone.)

## Flagged for human review
1. **CUDA-backend "mirrors the PTX backend" design-lineage comments (11).** All in the
   kept CUDA backend; they explain that the CUDA implementation copied the (now-removed)
   PTX backend's design, so the references are dangling but non-functional. Rewording each
   to preserve the rationale (vs. deleting the clause) is an authoring judgment I did not
   want to guess. Locations:
   - `CUDACommandQueue.java:94`
   - `CUDADeviceContext.java:104`, `:379`
   - `enums/CUDADeviceInfo.java:133`
   - `mm/CUDAMemorySegmentWrapper.java:219`
   - `mm/CUDAPinnedMemoryRegistry.java:43`, `:58`, `:62`, `:139` (these discuss a
     now-impossible "PTX backend in the same process also pinned this memory" scenario —
     may warrant more than a rename)
   - `scheduler/CUDAKernelScheduler.java:148`, `:174`
2. **README.md** — develop's README heavily markets the NVIDIA path and blends the PTX
   backend with the CUDA backend's real PTX codegen. I removed the two unambiguous items
   (the installable `5.2.0-ptx` SDKMAN row; a "four production backends (OpenCL, PTX,
   SPIR-V, Metal)" list that omitted CUDA → set to CUDA). I LEFT the blended NVIDIA
   marketing prose, which needs the author's wording decision:
   - `:8` "NVIDIA CUDA PTX … via CUDA/PTX … beyond generating PTX"
   - `:24` "emitted as CUDA PTX and compiled through NVRTC"
   - `:121` "more than a PTX code generator"
   - `:125` table row "**CUDA PTX backend** | … Graal IR → CUDA PTX → NVRTC → cubin"
   - `:129` "real PTX/CUDA generated from Java"
   - `:244` CUDA SDKMAN row "PTX codegen"
   - `:305` "JIT-compiles … to CUDA PTX"
   - `:323` "driven either through the **PTX/CUDA** backend … or … OpenCL"
   (Most are arguably accurate — the CUDA backend does emit PTX — but the "PTX backend" /
   "PTX/CUDA backend" phrasings read as a separate backend.)
3. **scripts/HOW_TO_BUILD_SDK.md** — I added `cuda` rows to match the build script (the
   develop doc listed only opencl/ptx/spirv). Confirm the added cuda rows are desired.

## Build

Command: `./mvnw -Pjdk21,opencl-backend,cuda-backend,spirv-backend -DskipTests -fae package`
(env setup: JDK 21; had to clean stale JDK-25 jars out of the gitignored `graalJars/`
so only the JDK-21 Graal 23.1.0 jars remained — otherwise "duplicate module on upgrade
module path". `graalJars/` is a build artifact, not part of the merge.)

Result: **ALL Java modules compile**, including every module touched by this merge —
tornado-api, tornado-runtime, tornado-drivers-common, tornado-drivers-opencl,
tornado-drivers-spirv, tornado-drivers-cuda, tornado-examples, tornado-benchmarks,
the cublas/cufft/cudnn/cusparse/cutlass Java wrappers, and **tornado-unittests**
(confirms the `TornadoVMBackendType.PTX` enum-removal ripple is fully fixed).

Failing modules — **environmental only, NOT merge-related** (merge touched no native code):
the 5 native JNI modules fail in the cmake/C++ phase on missing NVIDIA dev headers:
  - `tornado-drivers-cuda-jni` — `fatal error: nvrtc.h: No such file or directory`
  - `tornado-drivers-cublas-jni` — `cublas_v2.h` not found
  - `tornado-drivers-cufft-jni` — `cufft.h` not found
  - `tornado-drivers-cudnn-jni` — `cudnn.h` not found
  - `tornado-drivers-cutlass-jni` — CMake compiler-detection error (CUTLASS)
  (`tornado-drivers-opencl-jni` and `tornado-drivers-cusparse-jni` native builds passed.)
These need the cuBLAS/cuFFT/cuDNN/CUTLASS + NVRTC dev packages installed; install them
(or build with `BACKEND=opencl,spirv` only) to get a green native build. Not caused by
the PTX removal.

Metal backend not exercised (macOS-only); the sole Metal edit was a one-word string in
`MetalGraphBuilderPlugins.java` (PTX→CUDA in an error message).

# Merge decision log: upstream/develop → feat/purge-spirv-backend (PR #951)

Merge base: `287420465` (Bump version 5.1.1). Develop tip: `a5c1681e1` (104 commits;
**includes the merged PTX-removal PR #872**, so this merge also brings PTX removal in).
My branch: single commit `c4e94c832` "Remove the SPIR-V backend" (removes SPIR-V only;
adds no new functionality). My branch KEEPS `beehive-spirv-toolkit` (standalone SPIR-V lib).

Key structural fact driving nearly every conflict: **HEAD still had PTX** (my branch was cut
before develop removed it) while **develop still had SPIR-V**. The correct final state has
BOTH gone, leaving **OpenCL / CUDA / Metal**. So most conflicts resolved to "develop's side
(PTX already removed) minus SPIR-V".

Merge left in progress (NOT committed). 479 files, +11,414 / −38,145.

---

## SPIR-V-only from develop (dropped)
- `tornado-drivers/spirv/pom.xml`, `graal/compiler/SPIRVCompilerConfiguration.java`,
  `graal/compiler/plugins/SPIRVGraphBuilderPlugins.java` — SPIR-V backend files my branch
  deleted, develop modified. Kept deleted (`git rm`).
- No new SPIR-V driver files auto-merged in (verified 0 under `drivers/spirv/`).

## Non-SPIR-V from develop (accepted)
- `tornado-drivers/ptx/.../runtime/PTXTornadoDevice.java` — my branch modified it
  (SPIR-V-stripping); develop DELETED it as part of its own PTX-backend removal (#872).
  Honored develop's deletion (`git rm`). The rest of develop's PTX removal came in cleanly.

## Mixed files (surgically resolved — kept develop's non-SPIR-V changes, dropped BOTH PTX & SPIR-V)

### Backend registry / enums / drivers (Java)
- `runtime/.../enums/TornadoBackends.java` — enum body already {OpenCL,Metal,CUDA};
  removed both the PTX and SPIR-V priority imports (neither constant exists post-merge).
- `drivers-common/.../logging/Logger.java` — dropped the PTX/SPIRV enum entry → {OpenCL, Metal}.
- `drivers-common/.../TornadoDeviceQuery.java` — kept `CUDA→GREEN`, dropped PTX & SPIRV colour rows.
- `api/.../TornadoBackend.java`, `TornadoDeviceMap.java`, `TaskGraph.java` — javadoc backend
  lists → dropped PTX & SPIR-V, kept CUDA/OpenCL.
- `api/.../TornadoExecutionPlan.java` (×2) — "realised on CUDA backend; no-op for OPENCL,
  METAL" (dropped SPIRV from the no-op list; PTX already gone).
- `opencl/CUDACompiler/MetalCompiler` — the "FIXME Remove the inheritance (See PTX/SPIRV)"
  comment referenced now-removed backends → genericized to "Remove the inheritance".
- `runtime/.../TornadoOptions.java` (×3 comments), `TornadoKernelContextReplacement.java`,
  `interpreter/InterpreterUtilities.java` (if-else: kept `CUDA→GREEN`, dropped PTX & SPIRV),
  `interpreter/TornadoVMInterpreter.java` (comment; and `updateMeta` now sets only OPENCL
  compiler flags — develop set OPENCL+SPIRV, no CUDA line existed to keep. **See flag #3**).

### Build / CI / release
- `pom.xml` — dropped the 3 ptx/spirv `exportLists` `<arg>` lines (both export files removed;
  common/opencl/cuda/metal-exports remain).
- `tornado-drivers/pom.xml`, `tornado-assembly/pom.xml` — the conflicted profile was the
  `ptx-backend`(HEAD)/`spirv-backend`(develop) profile → removed the whole profile
  (opencl/metal/cuda-backend profiles remain).
- `Makefile`, `Makefile.mak` — dropped the `ptx:`/`spirv:` targets; clean target →
  `opencl,cuda[,metal]` (dropped both ptx-backend and spirv-backend); comment list → {opencl,cuda,metal}.
- `bin/compile` — clean = `opencl-backend,cuda-backend`; "full" set = `{opencl, cuda}`;
  `--backend` help = `opencl,cuda,metal`.
- `scripts/build-release-sdks.py`, `scripts/HOW_TO_BUILD_SDK.md` — dropped ptx & spirv build
  variants; "full" = `opencl,cuda` (matches bin/compile); label logic `{opencl,cuda}`.
- `.github/workflows/`: dropped the develop-added `build-spirv` job (build-test-jdk21);
  `BACKEND=…` lines → `cuda,opencl` and reverted "for SPIR-V with OpenCL runtime" step names
  (build-test-with-jdks21, graalvm-polyglot, heavy-memory-tests); release-sdk job display
  names → `opencl+cuda[+full]`; `1-prepare-release-jdk21/25.yml` sed alternation dropped `spirv|`;
  `ISSUE_TEMPLATE/bug_report.md` dropped the "Level Zero & SPIR-V Versions" line.
- `bin/intellij_build.py`, `bin/post_installation.py`, `bin/update_paths.py`,
  `bin/tornadovm-installer`, `tornado-assembly/src/bin/{gen-tornado-argfile-template,
  idea_xml_utils,tornado-benchmarks}.py`, `tornado.bat` — normalized backend
  lists/sets to `{opencl, cuda, metal}` (and "full" bundle = opencl,cuda).
- `tornado-assembly/src/bin/test-native.cmd` — removed the PTX and SPIR-V native-test blocks
  (they shared a closing paren); kept the OpenCL native-test block.
- `tornado-assembly/src/bin/tornado` (launcher) — removed a vestigial duplicate `*cuda*)`→
  `check_ptx_backend` case and the develop `*spirv*)`→`check_spirv_backend` case; the real
  `*cuda*)`→`check_cuda_backend` case remains.
- `tornado-assembly/src/bin/tornado.py` — took develop's version (PTX already removed) and
  stripped SPIR-V: removed `validate_spirv_backend()`, `__SPIRV_EXPORTS__`, the spirv/levelzero
  lib-existence checks (.so/.dylib), the `spirv-backend` export/module-flag blocks, and the
  spirv validation call site; kept validate_opencl/cuda/metal. Python syntax-checked OK.

### Docs
- `README.md` — intro backend list → NVIDIA CUDA/OpenCL C/Metal; SDKMAN table dropped the
  `spirv` row (took develop's 5.2.0 table); "four production backends (…SPIR-V…)" → "three
  production backends (OpenCL, CUDA, Metal)".
- `docs/source/{index,introduction,faq,flags,dev-tools,installation,developer-guidelines,
  ide-integration,multi-device,simple-start}.rst`, `INSTALL_FROM_SOURCE.md`, `CONTRIBUTING.md`,
  `HYBRID_API_GUIDE.md`, `tornado-assembly/src/etc/README.md` — took develop's prose, dropped
  SPIR-V (and Level Zero SYSMAN, the SPIR-V bullet in index, the SPIR-V/Level-Zero SDK rows),
  kept legit "compiled to PTX via NVRTC" (CUDA IR) mentions.

### Tests (Java) — several were compile-breaking
- **COMPILE-BREAKING (fixed):** develop-added tests referenced `TornadoVMBackendType.SPIRV`,
  which my branch removed from the enum. Fixed:
  - `TestBFloat16`, `TestHalfFloatInlineWrite`: `case OPENCL, SPIRV, METAL ->` → `OPENCL, METAL`.
  - `TestMatrixMultiplicationMMABF16/CpAsync/FP8`, `TestHalf2Packed` (×7): dropped the
    `assertNotBackend(TornadoVMBackendType.SPIRV)` guard lines.
- Conflict-resolved tests: `TestAtomics` (19 guards), `TestHello` (2), `TestInheritedFields`(4),
  `TestSignedComparisonsCodegen`, `TestVirtualDevice{FeatureExtraction,Kernel}`, the CUDA-lib
  tests `TestCuBlas/CuBlasLt/CuDnn/CuFft/Cusparse/Cutlass/FP8/Nvtx` (`case OPENCL, …, METAL`),
  `TestSimdgroupMatrix/TiledMatrix/MatrixPrimitives`, `TestConcurrentBackends` — dropped both
  PTX and SPIR-V guards/refs, keeping OpenCL/CUDA/Metal.
- `TornadoTestBase.java` — dropped develop's SPIR-V-only helpers (`getSPIRVSupportedDevice`,
  `assertNotBackendOptimization` using `SPIRVOptNotSupported`, `isSPIRVSupported`) — my branch
  had removed them and their supporting symbols.
- `TornadoHelper.java` — removed both the `TornadoVMPTXNotSupported` and `TornadoVMSPIRVNotSupported`
  imports (both exception classes deleted; their reporting blocks auto-merged out).
- `TestCompilerFlagsAPI.java` — removed the whole `testPTX`/`testSPIRV` method (tested a removed
  backend's compiler flags).
- `PrebuiltTests.java` — removed BOTH HEAD's PTX `testPrebuiltMutiBackend` and develop's three
  SPIR-V prebuilt tests (incl. `testPrebuilt04SPIRVThroughOpenCLRuntime`). **See flag #2.**
- `TestProfiler.java` — the `isBackendPTXOrSPIRV` dispatch-timer guard: renamed to
  `isBackendWithoutDispatchTimers` returning false (no remaining backend lacks dispatch timers;
  behaviour identical). **See flag #4** if you'd rather delete the now-vestigial guard.

## SPIR-V-as-IR retained per policy (NOT the removed backend) — CONFIRMED KEPT (PR #951)
My branch never touched these; they are the OpenCL/CUDA/Metal backends' ability to **load
prebuilt SPIR-V binaries** as an intermediate format (`clCreateProgramWithIL`,
`SPIRV_MAGIC_NUMBER = 119734787`, `isInputSourceSPIRVBinary`, `isSPIRVBinary`,
`createProgramWithIL(spirvBinary, …)`). Per the policy ("SPIR-V as an intermediate format is
not banned everywhere") and per PR #951 ("OpenCL/CUDA/Metal ability to load prebuilt .spv IL
binaries … backend-agnostic runtime capability, not the SPIR-V backend — deliberately kept"),
these were KEPT as-is. **Confirmed by author.** Files:
- OpenCL: `OCLCodeCache.java`, `OCLContext.java`, `OCLDeviceContext.java`, `graal/OCLInstalledCode.java`
- CUDA:   `CUDACodeCache.java`, `CUDAContext.java`, `CUDADeviceContext.java`, `graal/CUDAInstalledCode.java`
- Metal:  `MetalContext.java`, `MetalDeviceContext.java`, `graal/MetalInstalledCode.java`

## Flagged for human review
1. **SPIR-V-as-IR ingestion (above)** — CONFIRMED kept per PR #951 (deliberate). No action.
2. **`PrebuiltTests.testPrebuilt04SPIRVThroughOpenCLRuntime` dropped** — CONFIRMED correct: the
   test was removed in PR #951, so it stays dropped (no re-implementation needed).
3. **`TornadoVMInterpreter.updateMeta`** now propagates compiler flags for OPENCL only (develop
   propagated OPENCL + SPIRV; there was no CUDA/Metal line on either side). I removed the SPIRV
   line rather than guess a CUDA one — confirm CUDA/Metal compiler-flag propagation isn't needed
   here (pre-existing gap, not introduced by this merge).
4. **`TestProfiler` guard — DONE.** Removed the vestigial `isBackendWithoutDispatchTimers`
   method, unwrapped the 3 `if (...) {}` guards so the dispatch-timer asserts run
   unconditionally, and dropped the now-unused `driverIndex` locals. Recompiled: tornado-unittests
   BUILD SUCCESS.
5. **Backend-list normalization — REVISED to a pure strip (per author).** Reverted the earlier
   metal-normalization. Each backend list/set/enumeration is now **develop's value with only
   `spirv` (and any leftover `ptx`) removed** — metal is kept ONLY where develop already had it
   (`valid_backends`, `known_backends`, `__SUPPORTED_BACKENDS__`, and the metal-in-develop docs)
   and NOT added anywhere develop lacked it (`all_backends`, `detected_backends`, docstrings/
   examples, and the prose backend lists in flags/dev-tools/simple-start/multi-device/
   ide-integration/CONTRIBUTING/etc.). "full" bundle stays `opencl,cuda` (develop's `opencl,cuda,
   spirv` minus spirv), matching `bin/compile`. Also cleaned a develop-pre-existing duplicate
   `elif '-cuda'` branch in `bin/update_paths.py` (the vestige of the removed `-ptx` branch).

## Build
`./mvnw -Pjdk21,opencl-backend,cuda-backend -DskipTests -fae package` → **BUILD SUCCESS**.
All 27 modules compiled, including every module touched by this merge and, notably,
**tornado-unittests** (confirms the `TornadoVMBackendType.SPIRV` enum-removal ripple is fully
fixed) and all native JNI modules (cuda/cublas/cufft/cudnn/cusparse/cutlass-jni). 0 failures.
Metal backend not exercised (macOS-only); its edits were string/comment-level only.

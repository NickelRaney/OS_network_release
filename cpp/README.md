# Standalone C++17 / AVX2 Accelerator

This folder contains the self-contained C++17 / AVX2 implementation of the V1
layer 4C$\alpha$ network simulator (`v1_simulator.cpp`, `v1_4_neuron_kernel.hpp`,
`v1_l6.cpp`, `v1_rng.hpp`, `mat5.cpp`). It provides two ways to run the exact
same C++ engine:

1. **MATLAB MEX interface (`V1_4_fast_cpp_native`)**: Call the C++ engine
   directly from MATLAB (or via `cfg.backend = 'cpp_native'` in
   `run_simulation`), with no intermediate files or external processes.
2. **Standalone CLI executable (`v1_4_standalone`)**: Run directly from the
   command line on a MAT-v5 input file with no MATLAB runtime dependency.

On a modern x86-64 desktop CPU, a 2.5 s full-network simulation takes **~6.5 s**
in C++ versus **~110 s** in pure MATLAB (**~17x speedup**).

## Important Note on Random Number Generation (RNG)

This C++ engine uses its own self-contained native random stream (`v1_rng.hpp`):

- **Generator**: `xoshiro256** 1.0`, seeded via four `SplitMix64` outputs
- **Uniform draws**: 53-bit double-precision values in $[0, 1)$
- **Normal draws**: NumPy 2.3.3 double-precision Ziggurat transform

**It does NOT call MATLAB's `rand`/`randn` (`twister`), and the legacy
MATLAB-RNG MEX variant is not included in this repository.** As a consequence:

- For the **same seed**, individual spike trajectories from `cpp_native` /
  `v1_4_standalone` will **differ** from pure MATLAB (`V1_4_fast.m`).
- All model equations, double precision, update order, and population-level
  statistics match the MATLAB simulator.
- Repeated runs of the C++ engine with the same seed and build are
  deterministically reproducible (both CLI and MEX share the exact same engine
  and produce bit-identical outputs for the same inputs and seed).

## Requirements

- **CPU**: x86-64 processor with **AVX2** support
- **Compiler**: C++17 compiler (tested with MSVC 2022 on Windows x64;
  `CMakeLists.txt` and `build_v1_native_mex.m` also configure GCC / Clang with
  `-O2 -mavx2 -fno-fast-math -ffp-contract=off`)

Pre-compiled binaries are intentionally excluded from Git; build once locally.

## Option 1: Build and Use from MATLAB (`cpp_native` backend)

From MATLAB in the repository root:

```matlab
setup
build_v1_native_mex
```

Then either run through `run_simulation`:

```matlab
cfg = model_config('disordered_lgn');   % or 'low_lgn'
cfg.backend = 'cpp_native';
cfg.condition = 'onset';
cfg.duration_s = 2.5;
cfg.onset_s = 1.5;
cfg.analysis_start_s = 1.5;
folder = run_simulation(cfg);
```

or call `V1_4_fast_cpp_native` directly with MATLAB structs (`param_v1`,
`init_v1`, `ex`):

```matlab
[res, stats] = V1_4_fast_cpp_native(param_v1, init_v1, ex, uint64(2026112302));
```

`param_v1.f1_ee_scale` and `param_v1.f1_ei_scale` are respected automatically
(defaulting to `2.25` and `1.5` if omitted, matching `V1_4_fast.m`).

## Option 2: Build and Run the Standalone CLI (`v1_4_standalone`)

### Build with CMake (Windows / Linux / macOS x86-64)

```sh
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build --config Release
```

### Or build on Windows with MSVC PowerShell script

```powershell
.\cpp\build.ps1
```

### Run CLI

The CLI reads a MAT-v5 file containing `param_v1`, `init_v1`, and an external
LGN input struct (default variable name `ex_vert_long`, overrideable with
`--external`):

```sh
v1_4_standalone \
  --input path/to/input_case.mat \
  --external ex \
  --f1-ee 2.25 \
  --f1-ei 1.5 \
  --duration 2.5 \
  --seed 20261005 \
  --output result.mat \
  --summary summary.json
```

Use `--f1-ee 1.3 --f1-ei 1.2` for the Low-LGN model, or `--help` for all
options. Note that the lightweight MAT-v5 reader/writer supports standard
MAT-v5 files (not HDF5 / `-v7.3`).

## Third-Party Components and Licenses

- **miniz 3.0.2** (`third_party/miniz/`): MIT License (used for zlib
  decompression in the standalone MAT-v5 reader; license in
  [third_party/miniz/LICENSE](third_party/miniz/LICENSE)).
- **xoshiro256\*\* 1.0 & SplitMix64** (`vendor/`): Public domain (Blackman &
  Vigna; reference files in `vendor/xoshiro256starstar_reference.c` and
  `vendor/splitmix64_reference.c`).
- **NumPy 2.3.3 Ziggurat normal tables** (`vendor/ziggurat_constants.h`):
  BSD-3-Clause License (complete notice in
  [vendor/NUMPY_LICENSE.txt](vendor/NUMPY_LICENSE.txt)).

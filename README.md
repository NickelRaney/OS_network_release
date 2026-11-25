# OS Network

MATLAB implementation (plus an optional standalone C++17/AVX2 accelerator) of a
large-scale spiking network model of macaque primary visual cortex (V1 layer
4C$\alpha$) driven by magnocellular LGN inputs, with short-term synaptic
depression, GABA-B inhibition, and Layer 6 feedback.

The repository includes two fixed 4-hypercolumn network models (13,689
excitatory neurons, 2,401 inhibitory neurons, and 110 LGN afferents),
simulation entry points for background activity, drifting gratings, and
stimulus onset, steady-state and onset latency analysis tools, standalone
network-construction scripts, and optional C++ source code for ~17x faster
simulations.

## Included Models

| Preset | Alias | Network file | Complex-cell LGN inputs | F1 scales (EE / EI) |
| --- | --- | --- | --- | --- |
| `disordered_lgn` (default) | `high_cx` | `data/disordered_lgn.mat` | 2--5 (low), 7--8 (mid), 10--11 (high) | `2.25` / `1.5` |
| `low_lgn` | -- | `data/low_lgn.mat` | 0, 1, or 2 | `1.3` / `1.2` |

Both models use pre-generated fixed network realizations stored under `data/`
(no external downloads required). Standalone scripts illustrating how the LGN
templates and L4 torus connectivity are generated are available in
[scripts/](scripts/README.md).

## Requirements

- **MATLAB** (tested on **R2025b Update 2**, Windows x64; GNU Octave is not
  tested)
- **Statistics and Machine Learning Toolbox** (required for `exprnd` in `NLIF`)
- **Memory & runtime**: Allow at least **4--8 GB RAM**. On a desktop CPU, pure
  MATLAB takes ~15--20 s for the 0.4 s demo and ~1.5--2 min for a 2.5 s onset
  trial. With the optional C++ backend (`cpp/`), a 2.5 s simulation takes
  **~6.5--7 s** (requires an x86-64 AVX2 CPU and C++17 compiler).

## Quick Start

Open MATLAB in the repository root and run:

```matlab
setup
run_tests                               % Verify data, LGN gains, analysis, and scripts
folder = run_demo('disordered_lgn');    % 0.4 s grating demo (Disordered-LGN model)
folder_low = run_demo('low_lgn');       % 0.4 s grating demo (Low-LGN model)
```

Or run non-interactively from a terminal in the repository root:

```sh
matlab -batch "setup; run_tests; run_demo('disordered_lgn');"
```

Each simulation writes a new folder under `outputs/` containing:
- `result.mat`: full-network spike trains, central-ROI voltage/conductance
  traces, LGN inputs, parameters, and before/after `rng` states
- `run.json`: configuration, population firing-rate/ISI summaries, runtime, and
  SHA-256 checksums of code and data
- `raster.png`: spike raster of all central-ROI neurons (841 E in red, 145 I in blue)

## Running Custom Conditions

Configure stimulus conditions explicitly with `model_config` and
`run_simulation`:

```matlab
setup
cfg = model_config('disordered_lgn');   % or 'low_lgn'
cfg.condition = 'onset';                % 'bg', 'grating', or 'onset'
cfg.direction = [1 0];                  % [1 0] / [-1 0]: optimal (vertical); [0 1]: orthogonal
cfg.duration_s = 2.5;                   % 0.4 to 10 s, in multiples of 0.1 s
cfg.onset_s = 1.5;                      % must be 0 for 'bg' and 'grating'
cfg.analysis_start_s = 1.5;
cfg.lgn_seed = 2026112311;
cfg.v1_seed = 2026112312;
folder = run_simulation(cfg);
```

### Stimulus Direction Convention

The central analysis ROI (`inde(1,:)`) prefers **vertical** gratings.
`cfg.direction` specifies the grating drift direction perpendicular to the bars:
- **Optimal orientation (vertical)**: `[1 0]` (rightward) or `[-1 0]` (leftward)
- **Orthogonal orientation (horizontal)**: `[0 1]` (upward) or `[0 -1]` (downward)

## Optional C++ Accelerator (`cpp/`)

For multi-trial campaigns (such as steady-state or onset protocols), `cpp/`
provides the standalone C++17 / AVX2 implementation of `V1_4_fast` (~17x faster
than pure MATLAB, ~6.5 s for a 2.5 s trial), usable either as a MATLAB MEX
backend or as a standalone CLI executable (`v1_4_standalone`).

**Important difference in RNG**: This C++ engine uses its own self-contained
native RNG (`xoshiro256** 1.0` + `SplitMix64` + NumPy 2.3.3 double-precision
Ziggurat normal transform) and does **not** use MATLAB's `twister` RNG (the
legacy MATLAB-RNG MEX variant is not included). Therefore:
- Model equations, double precision, and population-level statistics match
  pure MATLAB (`V1_4_fast.m`).
- **Individual spike trajectories for the same seed will differ** between
  `cfg.backend = 'matlab'` and `cfg.backend = 'cpp_native'`.
- Repeated runs of the C++ engine with the same seed and build are
  deterministically reproducible.

To build and enable the C++ backend in MATLAB (requires an AVX2 x86-64 CPU and
C++17 compiler configured via `mex -setup C++`):

```matlab
setup
build_v1_native_mex                     % Build once into cpp/

cfg = model_config('disordered_lgn');
cfg.backend = 'cpp_native';             % Use C++ engine instead of 'matlab'
cfg.condition = 'onset';
cfg.duration_s = 2.5;
cfg.onset_s = 1.5;
cfg.analysis_start_s = 1.5;
folder = run_simulation(cfg);
```

To build the standalone MATLAB-free executable, see [cpp/README.md](cpp/README.md).

## Examples and Analysis Tools

- `examples/run_conditions('disordered_lgn')`: Runs single 2.5 s trials for
  `bg`, optimal `grating`, and `onset` in sequence.
- `examples/analyze_steady_state('disordered_lgn', 3)`: Runs BG, optimal
  forward (`[1 0]`), optimal reverse (`[-1 0]`), and orthogonal (`[0 1]`)
  gratings and computes firing rates (`compute_firing_rates`), modulation ratio
  $F_1/F_0$ (`compute_modulation_ratio`), directional selectivity (`DS`), and
  orthogonal/preferred (`OP`) ratio (`compute_selectivity`).
- `examples/run_onset_protocol('disordered_lgn', 'simple', 2, 4)`: Multi-phase
  onset workflow calling `compute_simple_onset` (sliding-window periodic
  detector) or `compute_complex_onset` (Gaussian-smoothed rate threshold
  detector).
- `examples/reproduce_historical_onset`: Evaluates `historical_t50` on the
  archived complex-cell histograms in `data/historical_onset_reference.mat`
  (741 exact phase-by-cell reference checks).

## Repository Layout

| Directory | Contents |
| --- | --- |
| `src/` | Core LGN/V1 simulators (`LGN3`, `LGN4`, `NLIF`, `ex_gen2`, `ex_gen_onset_test`, `V1_4_fast`, `V1_4_fast_cpp_native`), runner (`run_simulation`), and analysis functions |
| `cpp/` | Optional standalone C++17 / AVX2 simulator, MATLAB MEX bridge (`build_v1_native_mex.m`), CMake/MSVC build scripts, and third-party RNG/miniz notices |
| `config/` | Model presets and simulation defaults (`model_config.m`) |
| `data/` | Fixed network snapshots (`disordered_lgn.mat`, `low_lgn.mat`) and complex-onset histogram fixture |
| `examples/` | End-to-end scripts for conditions, steady-state metrics, onset protocols, and historical T50 check |
| `scripts/` | Standalone Section 2 network-construction utilities (LGN grid, templates, LGN/L4 wiring) |
| `tests/` | Automated unit and integration tests (`run_tests.m`) |
| `docs/` | Detailed parameter, output, and analysis documentation |
| `outputs/` | Local simulation outputs (ignored by Git) |

## Documentation

- [Model parameters and dynamics](docs/PARAMETERS.md)
- [Output files and data conventions](docs/OUTPUTS.md)
- [Analysis functions and stimulus conventions](docs/ANALYSIS.md)
- [Standalone C++ accelerator and RNG details](cpp/README.md)
- [Network construction scripts](scripts/README.md)

## License

Project code and data are released under the [MIT License](LICENSE). Third-party
components under `cpp/third_party/` and `cpp/vendor/` retain their respective
upstream licenses (`miniz`: MIT; `xoshiro256**` / `SplitMix64`: public domain;
NumPy Ziggurat tables: BSD-3-Clause).

# MATLAB V1 Disordered-LGN and Low-LGN models

An isolated research-code release candidate for the historical MATLAB V1
cortical-network project: 13,689 excitatory cells, 2,401 inhibitory cells and
110 LGN afferents. It includes the fixed network, visual/LGN input generation,
the MATLAB network simulator, two fixed model parameter sets, short demonstrations
and a high-cx frozen onset-analysis reference. No MEX compiler or C++ simulator is used.

The high-cx model itself contains simple cells and low/medium/high LGN-input
complex-cell groups. Those within-model analysis labels differ from the
separate Low-LGN model, now included as `low_lgn` with its own fixed network.

**Status: review candidate, not an archival reproduction of every paper
figure.** The original working directory mixed model variants in hard-coded
constants. The default `high_cx` preset uses the author-confirmed F1 scales
**2.25/1.5**, and the fixed Disordered-LGN/New-model network. The `low_lgn` preset
loads the original Low-LGN snapshot and uses its documented F1 **1.3/1.2**.
The `historical_disk` option is only a source-comparison
control, using the same high-cx connectivity with the old disk constants.
Authorship, citation and licensing still require confirmation. See [provenance](docs/PROVENANCE.md) and the
[publication checklist](docs/PUBLICATION_CHECKLIST.md).

Repository: <https://github.com/NickelRaney/OS_network_release>.

The [high-cx paper-alignment audit](docs/PAPER_ALIGNMENT.md) checks the
recent paper figures and their saved data, including a complete fresh onset
campaign. It finds unresolved Figure 5 input provenance and different archived
simple-latency arrays, plus measured differences in fresh firing rates.
**Implementation agreement does not certify agreement with every paper figure.**
Low-LGN source selection and validation are documented [separately](docs/LOW_LGN.md).

## Quick start

Open MATLAB in this repository's root directory:

```matlab
setup
report = run_tests;                 % Input checks + 741 historical t50 checks
folder = run_demo;                 % 0.4 s, complete network, fresh grating trial
low_folder = run_demo('low_lgn');   % Same entry point, separate Low-LGN network
S = load(fullfile(folder,'result.mat'));
S.summary
```

From a terminal whose current directory is the repository root:

```sh
matlab -batch "setup; run_tests; run_demo;"
```

Tested with MATLAB R2025b Update 2 on Windows x64. **Statistics and Machine
Learning Toolbox is required** for `exprnd` in the original LGN neuron model.
Parallel Computing Toolbox is not needed. MATLAB's bundled JVM is used for
SHA-256 run records. GNU Octave and other MATLAB releases have not been tested.
All required fixed input data are included under `data/`; no download or access
to the original experiment directory is needed for normal use.

The 0.4 s demo takes about 17 s for V1 simulation on the validation machine,
plus MATLAB startup, input generation, saving and plotting. Runtime depends on
hardware. It runs the full network and retains its dense connectivity and
temporary trace arrays; allow several GB of available RAM. Longer runs use
substantially more memory. The demo includes startup transients and is not a
steady-state firing-rate estimate.

## Choose an explicit model and condition

- `model_config('high_cx')`: `data/disordered_lgn.mat`; F1 2.25/1.5;
  BG/grating initial LGN gains 1.15/1.20.
- `model_config('low_lgn')`: `data/low_lgn.mat`; F1 1.3/1.2;
  BG/grating initial LGN gains 1.11/1.15.

Both use the historical onset helper with initial gains 1.11/1.15. All
connectivity, L6/F2 parameters and initial states come from the selected MAT
file; choosing `low_lgn` does more than change two gain coefficients.

```matlab
setup
cfg = model_config('high_cx');
cfg.condition = 'onset';
cfg.duration_s = 2.5;
cfg.onset_s = 1.5;
cfg.analysis_start_s = 1.5;
cfg.lgn_seed = 2026112311;
cfg.v1_seed = 2026112312;
folder = run_simulation(cfg);
```

`run_conditions('high_cx')` or `run_conditions('low_lgn')` runs one BG, one
direct grating and one onset trial of the selected model. This is a
three-condition example, not the paper's many-trial onset
protocol. Inspect [model_config.m](config/model_config.m) for the code-defined
gains as well as the fixed network. BG/direct grating use initial LGN gains
1.15/1.20 from the New-model notes. Onset retains its separately audited
20-spike curves with initial gains 1.11/1.15; these condition-specific rules
are not silently unified.

The audited visual setting is SF 2.5 cycles/degree, speed 4 degrees/second
(TF 10 Hz), a 0.1 ms integration step, and the saved 16/25 ms OFF/ON kernel
padding. The runner accepts normalized motion directions and a phase in units
of pi. It deliberately keeps SF, speed and dt fixed: the original LGN code has
assumptions about sampling and padding that need further validation before
these are generalized. Use durations from 0.4 to 10 s in multiples of 0.1 s.
Setting `kernel_extra_delay_ms=4` delays both kernels; no analysis offset is
then silently added.

## Outputs and reproducibility

Every run creates a fresh directory under `outputs/` with `result.mat`,
`run.json`, `raster.png` and copies of the source/configuration used. The MAT
file stores the full-network spikes, ROI traces, LGN spikes, visual drive,
actual parameters, seeds and before/after random-generator states. JSON records
the configuration, category rates, within-cell ISI summaries, MATLAB version,
simulation time and source/data SHA-256 hashes. See [output definitions](docs/OUTPUTS.md).

LGN and V1 use separate explicit `rng(seed,'twister')` resets, and trials run
serially. A seed is for replay of this release and environment. Historical
simulation seeds were not preserved, and bitwise agreement across MATLAB
versions or platforms is not promised. No historical MFE rejection, trial
resampling or cell filtering occurs in the new example runs.

For an analysis-only historical benchmark:

```matlab
setup
reference = reproduce_historical_onset;
```

This high-cx-only benchmark reuses compact archived histograms and reproduces 741 complex-cell
phase-specific latency values, including the historical binning/threshold
conventions. It does **not** regenerate those histograms from new simulation
spikes. [Validation](docs/VALIDATION.md) separates these two types of evidence.

## Repository contents

- `src/`: the MATLAB simulator, original input generators, extracted local
  helpers and small run/analysis utilities.
- `config/`: visible model constants and simulation defaults.
- `data/`: separate minimal high-cx and Low-LGN fixed networks, plus the
  high-cx frozen onset fixture.
- `examples/`: three-condition simulation and historical-analysis entry points.
- `tests/`: executable checks against independently evaluated input kernels and
  archived onset values.
- `tools/`: optional maintainer extraction and original-source comparison tools.
- `validation/`: machine-readable records from the local release validation.
- `docs/`: provenance, parameters, outputs, validation and publication decisions.

The dated experiment directories, interactive workspace scripts, manuscript
drafts, private correspondence, temporary files, generated results, native
binaries and later C++ experiments are outside this release. Network wiring is
distributed as a fixed realization; its original generation seeds and manual
editing history are not reconstructed here.

## Citation and reuse

Author names, paper title/identifier, software version and
license have not been established for this release. Do not cite an invented
paper or infer authors from a local account name. See
[CITATION.cff.in](CITATION.cff.in) for the metadata template and
[LICENSE_STATUS.md](LICENSE_STATUS.md) for the unresolved code/data permissions.

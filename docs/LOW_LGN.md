# Low-LGN model

The Low-LGN model is packaged alongside the original high-cx / Disordered-LGN
model at the author's request. Select it explicitly:

```matlab
setup
folder = run_demo('low_lgn');
S = load(fullfile(folder,'result.mat'));
S.cfg
S.summary
```

For a background or onset example:

```matlab
cfg = model_config('low_lgn');
cfg.condition = 'onset';
cfg.duration_s = 2.5;
cfg.onset_s = 1.5;
cfg.analysis_start_s = 1.5;
cfg.lgn_seed = 2026112301;
cfg.v1_seed = 2026112302;
folder = run_simulation(cfg);
```

`run_conditions('low_lgn')` runs BG, grating and onset examples in sequence.
The same portable MATLAB entry points serve both models; no original
experiment directory or MEX compiler is needed.

## Parameter source and distinctions

`data/low_lgn.mat` contains exactly these 11 variables from
`OS_network/25_9_28/lowlgn/low lgn full param.mat`: `param`, `param_v1`,
`init_v1`, `kernel3`, `grid`, `Cat`, `Cat_i`, `LVwire`, `LVwire_i`, `inde`,
and `indi`. Every variable was checked with MATLAB `isequaln` after export,
including numeric classes. The original file is unchanged. The paper's
`revised/fig5_icell_generate_data.m` identifies this same snapshot for its
Low-LGN model. Hashes and comparison results are in
`validation/low_lgn_parameters.json` and `docs/source_inventory.json`.

The original `25_9_28/parameters.txt` specifies:

- F1 EE/EI scales: **1.3/1.2**, giving initial multipliers about 1.041333/1.07.
- BG/direct-grating LGN initial E/I gains: **1.11/1.15**.
- Onset uses the shared historical `ex_gen_onset_test` helper, with initial
  gains **1.11/1.15** and its original 20-spike curves.

The 13,689 E / 2,401 I network retains its own E-to-E connectivity, LGN
wiring and saved L6/F2 parameters. The six base weights equal those of the
high-cx snapshot, but `cmee`, `Cat`, `LVwire`, `l6cof`, `seetfm`, `sei_list`
and `sincre_e` differ. Remaining fields are retained exactly, including
historical fields that the selected simulator does not consume.

The original kernels are shared: OFF/ON padding 16/25 ms, with no added
delay. The central analysis ROI has 841 E and 145 I cells. Of the 841 E cells,
594 are simple and 247 complex; complex cells have 0, 1 or 2 LGN inputs,
with counts **74, 81 and 92**. The Low-LGN Figure 4 nLGN and presynaptic E
count arrays match this extracted network exactly.

`historical_disk` remains a diagnostic preset using the high-cx connectivity
with the disk constants. Use `low_lgn` to load the actual Low-LGN network.

## Output groups

Low-LGN summaries contain All E, Simple, Complex, Cx 0 LGN, Cx 1 LGN,
Cx 2 LGN and I. The Complex aggregate includes all three LGN subgroups.
High-cx retains its existing summary labels and definitions. Read
`summary.labels` and `summary.cell_count`; do not assume identical category
positions across models.

## Validation boundary

The release checks exact parameter extraction, Low-LGN network identity,
LGN connection totals, zero-LGN cell accounting, paper Figure 4 network
descriptors, and original-source simulation agreement. The fresh examples
exercise BG and grating for 0.4 s and onset for 2.5 s with a 1.5 s onset.
The short direct-condition runs are dependency checks, not steady-state
firing-rate estimates.

The original-source comparison runs the unchanged working-copy MATLAB
simulator and input generators with the Low-LGN snapshot and explicit seeds.
The 0.4 s comparison passed with all eight output fields and final RNG state
exactly equal. The three fresh examples completed and passed spike-data checks.
For a maintainer with access to the original archive:

```matlab
setup
addpath('tools')
report = verify_against_archive(archive_root, 'low_lgn');
```

Machine-readable evidence is in `validation/low_lgn_archive_comparison.json`,
`dual_model_checks.json`, `high_cx_summary_regression.json`, and
`low_lgn_grating_run.json`, `low_lgn_bg_run.json`, `low_lgn_onset_run.json`.
The original high-cx summary was also compared against the pre-extension
implementation using an existing full-network trial.

A Git-index export into a separate path containing spaces also passed the
tests and ran `run_demo('low_lgn')` in a fresh MATLAB process. Full-network
spikes, LGN spikes and drive, RNG states, actual parameters, configuration and
summaries matched the original demo exactly. See
`validation/low_lgn_portable_checkout.json`.

The Low-LGN grating/BG/onset runs took approximately 14.88/11.31/92.54 s
for V1 simulation on the validation machine. Their E/I mean rates were
10.6926/49.3448, 2.0482/15.4310 and 25.9917/77.2966 Hz in the explicitly
recorded windows. These single-trial values are execution records, not a fit
to the paper or a steady-state estimate from the short trials.

This extension does not certify full Low-LGN paper-figure reproduction.
The packaged historical onset fixture and the earlier 220-slot campaign
belong to high-cx only. They must not be relabeled as Low-LGN validation.
Original Low-LGN raw trial collections and figure-reconstruction scripts
remain outside the public payload until their provenance and protocols are
separately reconciled. Author, citation, license and data-permission fields
remain deferred as requested.

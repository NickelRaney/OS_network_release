# Output Files and Data Conventions

Every call to `run_simulation(cfg)` creates a fresh directory under `outputs/`
containing:

- `result.mat`: Full simulation output (MATLAB `-v7.3` format)
- `run.json`: Human- and machine-readable configuration, summary rates, ISI
  statistics, runtime, and SHA-256 hashes of code and data
- `raster.png`: Spike raster of all central-ROI neurons (E in red, I in blue)
- `src/` and `config/`: Exact copy of the code and configuration used for the run

Generated folders under `outputs/` are excluded from Git.

## Units and Indexing

- **Time units**:
  - `cfg.duration_s`, `cfg.onset_s`, `cfg.analysis_start_s`, V1 spike times
    (`spike`), `trace_time_s`, and `l6_block_start_s` are in **seconds**.
  - `dt` (`0.1 ms`), refractory periods, LGN spike times (`lgn`), and membrane /
    synaptic time constants are in **milliseconds**.
- **V1 spikes (`spike`)**: One column per neuron: E cells `1:13689`, followed
  by I cells `13690:16090`. Row 1 stores the spike count `n`; rows `2:(n+1)`
  store spike times in seconds; remaining rows are zero.
- **LGN spikes (`lgn`)**: `110 x N` matrix, one row per afferent. Column 1
  stores the spike count; subsequent entries store spike times in milliseconds.
  `out` stores the filtered 1 ms visual drive.
- **Central ROI**: `roi_e` (841 global E indices) and `roi_i` (145 global I
  indices, offset by `13689`).

## Variables in `result.mat`

| Variable | Description |
| --- | --- |
| `cfg` | Run configuration struct (preset, condition, direction, phase, duration, onset, seeds, gains) |
| `param`, `param_v1` | Effective LGN and V1 parameter structs used by the run |
| `kernel3` | Effective LGN spatial and temporal kernels |
| `lgn`, `out` | LGN spike times (ms) and filtered visual drive (1 ms sampling) |
| `random_state` | MATLAB `rng` states before and after LGN (`lgn_before`, `lgn_after`) and V1 (`v1_before`, `v1_after`) |
| `spike` | Compact full-network spike matrix (trimmed to `max_count + 1` rows) |
| `roi_e`, `roi_i` | Central ROI E and I global indices |
| `trace_time_s` | 1 ms sample timestamps for ROI membrane/conductance traces |
| `l6_block_start_s` | Block start times for 1 ms L6 event counts |
| `traces` | ROI traces: `VE`, `VI` (normalized voltage), `HE`, `HI` (summed synaptic conductances), `L6` (proposed L6 event counts per 1 ms block), `ambient`, and `LGN` |
| `summary` | Per-population mean firing rates and within-cell ISI statistics in `(analysis_start_s, duration_s]` |
| `meta` | MATLAB version, platform, V1 wall-clock time, and SHA-256 manifest of code and data |

## Population Labels in `summary`

- For `disordered_lgn` (and alias `high_cx`):
  `{'All E', 'Simple', 'Cx low', 'Cx medium', 'Cx high', 'I'}` with ROI counts
  `[841, 594, 58, 86, 103, 145]`.
- For `low_lgn`:
  `{'All E', 'Simple', 'Complex', 'Cx 0 LGN', 'Cx 1 LGN', 'Cx 2 LGN', 'I'}` with
  ROI counts `[841, 594, 247, 74, 81, 92, 145]`.

Always index `summary` using `summary.labels` rather than hard-coded positions.
ISIs are computed strictly within individual neurons inside the analysis window;
groups with no intervals report `NaN` (`null` in JSON).

# Output definitions and numerical conventions

## Units and ordering

- `param.duration`, `param_v1.duration`, `cfg.duration_s`, V1 spike times and
  analysis windows are in seconds. The integrator's `dt`, refractory periods,
  LGN spike times and membrane time constants use milliseconds.
- `spike` has one column per cell: E IDs 1:13689 followed by I IDs 13690:16090.
  Row 1 is the count, rows 2:(count+1) contain times in seconds; remaining rows
  are zero. Always respect the count row.
- `lgn` has one row per afferent. Column 1 is its spike count; subsequent values
  are times in milliseconds. `out` is the filtered visual drive at 1 ms
  sampling and includes the original helper's padding beyond the simulated
  interval.
- The central ROI is the snapshot's `inde(1,:)` / `indi(1,:)`: 841 E and 145 I
  for both packaged networks. I global IDs include the E offset. This is not the
  entire cortical orientation domain.
- Connectivity is presynaptic row / postsynaptic column:
  `cmee`: E->E, `cmie`: E->I, `cmei`: I->E, `cmii`: I->I.

## Saved variables

`cfg` records preset, condition, gains, duration, direction, phase, delays,
analysis start and independent LGN/V1 seeds. `param` / `param_v1` contain the
actual run parameters, including the fixed connectivity. `kernel3` is the
actual kernel used. The input profile is determined exactly by cfg and the
saved source; the LGN visual drive and spikes are also saved.

`random_state` holds `lgn_before`, `lgn_after`, `v1_before`, `v1_after`, each a
MATLAB `rng` state. Save these alongside seeds when implementing continuations;
this runner always starts a fresh simulation and does not implement restart
from final network state.

`traces.VE` / `VI` contain normalized membrane voltage in ROI order.
`HE` / `HI` contain the model's summed excitatory / inhibitory synaptic
conductance terms for E followed by I. `ambient` and `LGN` are the E external
input components from the simulator. Do not label voltage in mV or conductance
in physical units without the model's normalization conventions.

`trace_time_s` corresponds to integrator steps 10,20,... (1 ms apart at
dt=0.1 ms). The original allocated final unused row is removed when saving ROI
traces. `traces.L6` instead contains the count of proposed L6 events in each
ten-step block, as returned by the original simulator; its block start times
are `l6_block_start_s` (0.2,1.2,... ms). These counts are not a conductance or a
sample of instantaneous L6 rate, and later facilitation can weight delivery.

`summary` uses raw spikes in `(analysis_start_s, duration_s]`. Each group's
mean rate is its total count divided by cell count and elapsed seconds.
For `high_cx` and the diagnostic `historical_disk` preset, classes follow the
historical central-ROI definitions: simple = `Cat(:,1)~=2`;
complex low = `Cat(:,1)==2` and 2--5 LGN afferents; medium = 7--9; high = >9.
The main snapshot has 594/58/86/103 E cells in those groups. ISIs are differences
between adjacent spikes of the same cell with both spikes inside the window.
For `low_lgn`, `summary.labels` instead contains All E, Simple, Complex,
Cx 0 LGN, Cx 1 LGN, Cx 2 LGN, I, with counts 841/594/247/74/81/92/145.
The Complex row is the aggregate; its three LGN-count subgroups partition it.
Zero-LGN complex cells are included, rather than omitted by high-cx grouping.
Consumers should read the labels instead of assuming fixed column positions.

ISIs are never formed across cells or trials; a group with no intervals has
NaN median (JSON null). The pooled ISI median weights intervals, not cells.

`raster.png` shows every raw ROI spike in the full simulated interval. E cells
are red and retain their original order; I cells are blue. The title specifies
the separate window used for rate summaries. No raster binning, subsampling or
category reordering changes the stored spike data.

`meta` records the version, platform, V1 wall-clock time and hashes of code and
fixed data. The runtime source/config copies are retained under the result
directory. The top-level entry scripts are identified by hashes and supplied
in the repository. Generated output directories are excluded from Git.

## Scope

BG means zero visual stimulus with active stochastic LGN, ambient and L6 input.
`grating` denotes the fixed direction/phase trial; it is not automatically
declared the empirically optimal direction. There is no scan selecting "opt"
or "ortho" in this small release example. `onset` is a global grating onset
after a background interval. One such trial cannot estimate the paper's
many-trial simple or complex onset distributions. Spectral analysis and the
later delay/adaptive-threshold experiments are not part of this release.

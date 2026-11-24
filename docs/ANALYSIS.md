# Analysis Functions and Stimulus Conventions

All analysis functions live in `src/` and are added to the MATLAB path by
`setup`.

## Stimulus Direction Convention

The central analysis ROI (`inde(1,:)`, 841 E cells) prefers **vertical**
gratings. In `cfg.direction`, the 2D unit vector specifies the **drift
direction** perpendicular to the grating bars:

- **Optimal orientation (vertical grating)**: `cfg.direction = [1 0]` (rightward)
  and `cfg.direction = [-1 0]` (leftward).
- **Orthogonal orientation (horizontal grating)**: `cfg.direction = [0 1]`
  (upward) or `[0 -1]` (downward).
- **Background (`bg`)**: Zero visual contrast with active spontaneous LGN,
  ambient, and L6 activity.

## Steady-State Metrics

### Firing Rates (`compute_firing_rates`)

```matlab
rates = compute_firing_rates(spike, [1.0, 3.0]);
```

Returns per-cell mean firing rates (Hz) in `(window_s(1), window_s(2)]` across
all 16,090 neurons (`rates.all`, `rates.e`, `rates.i`).

### Modulation Ratio (`compute_modulation_ratio`)

```matlab
mr_fwd = compute_modulation_ratio(spike_fwd, [1.0, 6.0], 10);
% Or evaluate both opposite optimal directions and select each cell's preferred direction:
mr = compute_modulation_ratio(spike_fwd, [1.0, 6.0], 10, spike_rev);
```

Bins spikes at 10 ms resolution over whole seconds, folds across 1 s epochs,
takes the FFT, and computes $\mathrm{MR} = F_1 / F_0 = 2 |F(f_{\text{TF}})| / |F(0)|$
at the grating temporal frequency $f_{\text{TF}} = 10\text{ Hz}$. When both
opposite directions are supplied, `mr.mr_pref` selects the modulation ratio from
the direction yielding the higher mean firing rate (`F0`).

### Directional Selectivity and OP Ratio (`compute_selectivity`)

```matlab
sel = compute_selectivity(fr_fwd, fr_rev, fr_ortho, Cat, roi_e, 'disordered_lgn');
```

- **Directional Selectivity (DS)**:
  $\mathrm{DS} = \max(r_{\text{fwd}}, r_{\text{rev}}) / \min(r_{\text{fwd}}, r_{\text{rev}})$
  (uncapped raw values in `sel.ds_raw`; capped at 9 when $>8$ or `NaN` in
  `sel.ds_capped`, matching the paper's histogram binning).
- **Orthogonal-to-Preferred (OP) ratio**:
  $\mathrm{OP} = r_{\text{ortho}} / \max(r_{\text{fwd}}, r_{\text{rev}})$.

Returns per-cell arrays and per-group medians for Simple, Complex, and Complex
LGN subgroups.

See [examples/analyze_steady_state.m](../examples/analyze_steady_state.m) for a
complete runnable workflow.

## Onset Latency Estimators

### Simple Cell Onset Latency (`compute_simple_onset`)

```matlab
out = compute_simple_onset(spike_trials, onset_s, period_ms);
```

Implements the periodic sliding-window algorithm for simple cells:
1. For each stimulus phase (4 trials per phase), folds post-onset spikes modulo
   `period_ms` (100 ms for TF = 10 Hz) into 5 ms bins and finds the 35 ms
   (7-bin) window with the highest total spike count.
2. Steps forward cycle by cycle (`+100 ms`) from stimulus onset until the
   earliest cycle in which at least 3 of 4 trials fire within that 35 ms window.
3. Averages the spike times inside that earliest window to obtain the phase
   latency, and takes the minimum across stimulus phases for each cell.

### Complex Cell Onset Latency (`compute_complex_onset`)

```matlab
out = compute_complex_onset(spike_phase_groups, onset_s, 'ThresholdFraction', 0.5);
```

Implements the instantaneous-rate threshold estimator for complex cells:
1. Aggregates spikes across trials within each phase group into 1 ms bins
   spanning `[onset_s - 1.0, onset_s + 1.0]` (2,000 bins) and convolves with a
   21-point normalized Gaussian kernel ($\sigma = 5\text{ ms}$).
2. Computes background rate $f_{\text{bg}}$ (first 1,000 ms) and steady-state
   rate $f_{\text{ss}}$ (last 801 ms, bins 1200--2000).
3. Finds the earliest post-onset time $t \in [1, 300]\text{ ms}$ such that the
   minimum instantaneous rate over the subsequent 25 ms strictly exceeds
   $f_{\text{bg}} + a (f_{\text{ss}} - f_{\text{bg}})$ (default $a = 0.50$;
   `301` indicates no crossing). Takes the minimum across phase groups.

### Frozen Complex-Onset Benchmark (`historical_t50` / `reproduce_historical_onset`)

`reproduce_historical_onset` runs `src/historical_t50.m` on the compact
archived complex-cell histograms in `data/historical_onset_reference.mat` to
verify exact agreement on 741 phase-by-cell reference values.

See [examples/run_onset_protocol.m](../examples/run_onset_protocol.m) for a
fresh multi-trial onset workflow.

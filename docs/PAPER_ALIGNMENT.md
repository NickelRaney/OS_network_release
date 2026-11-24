# Alignment with the recent paper figures

**Status: not certified as a complete reproduction of the paper figures.**
This audit concerns only the high-cx / Disordered-LGN network with the
author-confirmed F1 scales 2.25/1.5. It does not change model equations or
choose parameters to improve agreement.

The later Low-LGN packaging extension is documented in [LOW_LGN.md](LOW_LGN.md).
The high-cx results below must not be attributed to that separate model.

The dated baseline covers figures modified from August 23 through November 23,
2026: the September 3 Figure 4 in `fig4/`, revised Figures 5 and 6, and the
September supplementary figures. The older MAT files loaded by these figures
remain necessary dependencies. File modification dates do not establish when
a stochastic simulation was run.

## What is reproduced from the saved data

- Figure 4's LGN counts and presynaptic E counts match the packaged network.
- Reanalysis of the original spike data reproduces all 274 simple latencies
  to floating-point roundoff (maximum error 2.3e-13 ms), all 247 complex
  latencies exactly, and all 247 corresponding complex steady rates exactly.
- The six saved Figure 3 raster examples match the Figure 5 onset spikes
  exactly when the original cell-selection rule is applied.
- The supplementary complex-rise plot uses unshifted high/medium/low medians
  42/51.5/65 ms. Figure 4 uses the same values with its display-only +4 ms
  offset, giving 46/55.5/69 ms.

These are historical-data analysis checks, not fresh reproduction of the
original stochastic trajectories. Figure 1 includes anatomical/reference
illustrations. Figure 3's F1/F0 panel has not been independently certified.

## Accelerated fresh runs

A separate local C++ backend uses MATLAB's random-number stream and preserves
the release's operation order. All eight output fields and final RNG state
match the release exactly for both a 0.4 s grating run and a 2.5 s onset run.
The C++ code is an audit accelerator outside this MATLAB-only public payload.
See `validation/paper_backend_short.json` and `paper_backend_onset.json`.

The audit runs four 3 s trials each for BG, preferred grating and orthogonal
grating, with rates counted strictly between 1 and 3 s. All-E means in the
saved Figure 6 versus fresh trial means are:

- BG: 5.0089 versus 5.4759 Hz (+9.32%).
- Opt: 28.2515 versus 28.8343 Hz (+2.06%).
- Ortho: 12.2883 versus 11.4261 Hz (-7.02%).

The broad rate pattern is similar. Four new trials, without the historical
RNG states or complete historical run metadata, do not establish statistical
equivalence. Population-specific results and trial ranges are included in
`validation/paper_alignment.json`.

The complete fresh onset campaign uses 180 complex-cell trials (3 phase
groups x 60) and 40 simple-cell trials (10 phases x 4), with the original
phase jitter and fixed 274 preferred-direction simple cells. It retains
first attempts and all retries needed by the historical pre-onset MFE rule.
Both first-attempt and accepted-trial statistics are reported in
`validation/paper_onset.json`; the +4 ms display offset is explicitly applied
when comparing with Figure 4. This offset is not added to the input kernels.

There were 241 onset attempts in total. First attempts failed the historical
MFE criterion in 18/180 complex slots and 2/40 simple slots. All attempts were
retained. After the historical rule, simple / Cx high / Cx medium / Cx low
median latencies are **53.85 / 47 / 62.5 / 77.5 ms**, versus the paper's
**52.66 / 46 / 55.5 / 69 ms**. Without rejection they are
52.98 / 48 / 64 / 78 ms. Thus the low- and medium-LGN-input complex groups in
this high-cx model remain slower than the archived reference in this campaign.

## Unresolved Figure 5 provenance

Figure 5 panel E stores a different simple-cell latency array: its displayed
median is 76.4625 ms, while Figure 4's is 52.66 ms. Adding the Figure 4 offset
to Figure 5's data does not resolve the difference (maximum per-cell
difference 168.4667 ms). The intended relationship between these archived
arrays must be established before treating them as one reference result.

The packaged parameter snapshot's saved `lgn_vert` and `lgn_vert_long` do not
reconstruct Figure 5's frozen LGN conductance traces with the release input
generators. In the selected cell (ROI 662, global E index 7667), maximum
errors within the checked plot windows are 0.05383 and 0.07126 in model
conductance units for steady and onset input, respectively. This check occurs
before stochastic V1 simulation, so a missing V1 seed cannot explain it.
It indicates an incomplete input provenance chain, not a demonstrated cause
in a particular historical simulator implementation.

The alternative `DS_all.mat` snapshot also fails this input reconstruction;
its effective E LGN wiring differs from the packaged snapshot. It is not
silently substituted into the release.

## Review boundary

No original figures or experimental files were overwritten. Every fresh
trial and rejected attempt, its LGN spikes and RNG states, actual parameters,
source hashes, audit code and comparison plots are retained in a separate
local audit workspace. The portable payload contains compact evidence only.
Authorship, citation, licensing and data-publication decisions remain
deferred. This audit preceded the subsequently requested GitHub upload;
repository availability does not change its scientific validation limits.

Before certifying all paper figures, resolve the Figure 5 input sources and
the two simple-latency arrays, then decide whether the observed fresh-run
differences are compatible with the intended scientific reproducibility
criterion. Exact source equivalence alone is insufficient for that claim.

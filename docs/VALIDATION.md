# Validation of this release candidate

Validation was performed locally on 2026-11-23 with MATLAB
25.2.0.3055257 (R2025b) Update 2, Windows x64. The machine-readable evidence is
in `validation/`. These are checks of the MATLAB package and selected archived
analysis, not a claim that every historical paper result has been regenerated.

## Fixed data and dependency checks

The selected 11 variables were exported from the original parameter snapshot
and loaded again; every retained variable passed `isequaln`, including its
classes, connectivity, weights, initial voltages and kernels. The compact
network file is about 6.5 MB, versus the original 69 MB workspace archive.
No saved `ex_vert_long` is used as simulation input.

`run_tests` checks network dimensions, ROI counts, saved OFF/ON padding, LGN
connection totals, independent sums of AMPA impulse responses and all archived
complex-cell onset values. It separately exercises the disk and high-cx
steady gains and verifies that high-cx steady settings do not alter the
onset curves. Full measured errors are in `input_and_onset_checks.json`.

`dependencies.json` records MATLAB's static dependency closure and required
products. The examples were also executed using only this release's paths.
The Python payload checker rejects absolute machine paths, accidental native
binaries/experiment spreadsheets and files larger than its review threshold
of 50 MB. This threshold is a project packaging check, not a statement about
hosting-service policy. `release_manifest.json` hashes the published payload.

A Git-index export into a separate directory whose name contains spaces also
passed `run_tests` and `run_demo`. Its complete saved LGN arrays, full-network
spikes and before/after RNG states were exactly equal to the original high-cx
demo. Runtime function resolution pointed into that exported copy; see
`portable_checkout.json`. `source_integrity.json` confirms that the 18 original
source/evidence files still match their recorded hashes.

## Original-source equivalence

`tools/verify_against_archive.m` works with a read-only path to the original
archive. It copies the original MATLAB functions into a local scratch directory
and changes only their public names to avoid path shadowing. The original
files remain untouched. It compares the two input generators exactly and
runs the full network for 0.4 s with LGN seed 2026112301 and V1 seed 2026112302.

Using the explicit disk-control F1 values, all eight simulator fields (`VE`,
`VI`, `HE`, `HI`, `L6`, `amb`, `lgn`, `spike`) and final MATLAB RNG state are
**exactly equal**. The final original/release calls took approximately 17.2/17.1 s.
See `archive_comparison.json`. This checks packaging changes against the
available source; it does not endorse the disk-control constants as high-cx.

`tools/verify_high_cx.m` constructs a second reference by changing only the
two original F1 literals to the author's confirmed 2.25/1.5. The release and
this corrected reference again have **exactly equal** eight output fields and
final RNG state over 0.4 s, using high-cx LGN input. See
`high_cx_comparison.json`. This is a reference implementation check of the
authorized correction, not a new fit to archived paper spikes.

## Fresh high-cx examples

`demo_run.json`, `bg_run.json` and `onset_run.json` record the actual high-cx
configuration, seeds, source/data hashes, rates and runtime for:

- 0.4 s direct grating, full network, window (0,0.4] s.
- 0.4 s background, full network, window (0,0.4] s.
- 2.5 s onset trial with stimulus at 1.5 s, window (1.5,2.5] s.

The saved full-network spikes pass finite-value, count, ordering, duration and
0.1 ms grid checks. Complete results and copied sources are local generated
outputs; the public payload retains the small JSON records and a demo raster.
Rates from short trials are not steady-state estimates, and no trial was
rejected or rerun for a favorable scientific outcome.

Measured V1 times were 17.78 s (grating), 14.31 s (BG) and 105.01 s (onset).
The corresponding ROI E/I means were 17.241379/68.689655 Hz,
3.406659/22.172414 Hz and 27.605232/84.448276 Hz in the stated windows.
These are deterministic replay checks for the documented seeds, not population
estimates from repeated trials. The demo raster was visually inspected.

Source files are normalized to LF line endings with trailing whitespace
removed, and `.gitattributes` preserves that format across checkouts. The
format-only transformation is recorded in `source_formatting.json`; the
original-source comparisons and a portable demo replay were run on the final
formatted code. Earlier BG/onset records retain their actual run-time hashes,
which can be mapped to the formatted files through that record.

## Archived analysis benchmark

`reproduce_historical_onset` reproduces all **741/741** archived complex-cell
phase-specific T50 values with **0 ms maximum error**. The 247 cells belong to
the same high-cx network: 103 high-, 86 medium- and 58 low-LGN-input complex
cells. After the original display-only +4 ms shift, medians are 46, 55.5 and
69 ms. These group labels are not separate network models.

The fixture consists of frozen histograms, so this validates the estimator,
not historical LGN/network simulation. Historical seeds and the exact runtime
source were not archived alongside those spike results. The previous native
reproduction audit explicitly reported unresolved differences; this release
does not erase that limitation.

## Low-LGN extension

The separate Low-LGN snapshot and its runtime configuration are now packaged
and validated independently. See [Low-LGN validation](LOW_LGN.md),
`low_lgn_parameters.json`, `low_lgn_archive_comparison.json`,
`dual_model_checks.json` and the three `low_lgn_*_run.json` records. The
high-cx paper campaign above is not a Low-LGN reproduction claim.

## Limits and execution notes

- A subsequent separate audit ran the full 220-slot onset protocol and
  twelve fresh BG/opt/orthogonal trials; see [paper alignment](PAPER_ALIGNMENT.md).
  It also validates the historical simple-cell analysis. It does not certify
  every paper figure. The public examples remain lightweight and do not ship
  the full historical raw trial collection or the separate audit accelerator.
- No direction-selection scan, spectral analysis or new network-generation
  pipeline is supplied. Low-LGN is included as a fixed network; a complete
  Low-LGN paper-figure reproduction campaign has not been performed.
- Only the documented dt/SF/speed configuration was exercised. Other MATLAB
  releases, operating systems and Octave have not been validated.
- Two initial analysis-only attempts exposed integer types introduced by
  SciPy when exporting fixture metadata. The export now uses double reference
  values and the reader casts numeric metadata explicitly. The complete checks
  were rerun after the fix; scientific simulation code was not altered by it.
- Bitwise comparisons are local to the tested environment. They should not
  be used as a cross-platform guarantee or as proof of original paper seeds.
- Licensing and bibliographic completeness are author decisions still listed
  in PUBLICATION_CHECKLIST.md; numerical validation does not resolve them.

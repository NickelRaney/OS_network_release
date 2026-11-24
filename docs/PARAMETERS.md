# Model parameters and their sources

The starting points are two fixed historical parameter snapshots, not newly
generated networks. `data/disordered_lgn.mat` and `data/low_lgn.mat` each retain
the 11 named runtime variables documented in PROVENANCE.md. Use `load` with an output structure;
loading the whole original workspace into the base workspace is unnecessary.

## Network and base weights

The fixed network has 13,689 E cells on a 117 x 117 grid, 2,401 I cells and
110 LGN afferents. `Cat`, `Cat_i`, `LVwire` and `LVwire_i` retain the cell classes
and LGN connections. `grid` determines LGN spatial phase. `inde` / `indi`
select the historical analysis blocks. The runner uses their first block.

Both snapshots' base weights in `param_v1` include `see=0.0245`, `sie=0.0176`, `sei=0.051`,
`sii=0.09`, `selgn=0.0529`, and `silgn=0.09881`. They are not overwritten by
F1 or LGN gain configuration. `L6tfm`, `l6cof`, `seetfm`, `sei_list`, the
refractory vectors, `mask`, `seifac`, `segabab` and the remaining fields are
loaded from the same snapshot. Consult the actual simulator for the subset
of retained historical fields it consumes; unused fields are not promoted to
active model parameters merely because they are present in the MAT structure.

The Low-LGN snapshot has its own `cmee`, `Cat`, `LVwire`, `l6cof`, `seetfm`,
`sei_list` and `sincre_e`. Its 247 central complex cells have 0/1/2 LGN
inputs (74/81/92 cells). The high-cx complex groups have 2--5, 7--9 or >9
inputs. Do not transplant high-cx connectivity or lookup tables into a Low-LGN
run merely because the scalar base weights agree.

`init_v1` contains saved E/I voltages. The original simulator still generates
initial synaptic states, ambient activity and L6 events from its V1 random
stream. This is not a saved full network-state continuation.

## Code-defined dynamics

At the audited `dt=0.1 ms`, the original MATLAB simulator uses E/I membrane
leak times 20/15 ms, base normalized spike threshold 1, reset voltage 0 and
E adaptive-threshold increment 0.5 with 5.05 ms decay. Excitatory/inhibitory
reversal values in the voltage equations are 14/3 and -2/3. Time constants
for rise/decay are AMPA 0.5/2 ms, NMDA 2/80 ms, GABA-A 0.5/5 ms and
GABA-B 3/100 ms. These are in `src/V1_4_fast.m`, not inferred from every
similarly named snapshot field.

The high-cx F1 curves are evaluated literally as

```matlab
e_decay = ([1.2 1.2 1.2 1.20 1.19 1.16 1.12 1.07 1 0.96 0.94 0.92]-0.92)*2.25/3+0.92;
i_decay = ([1.1 1.02 0.98 0.96 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95]-0.95)*1.5/1.5+0.95;
```

Their initial multipliers are 1.13 and 1.10. The recurrence, recovery rules,
F2 lookup behavior and interaction with base weights are inherited from
`V1_4_fast`. There are no new recurrent synaptic delays or later native-model
threshold modifications in this package.

For `low_lgn`, the same two expressions use 1.3 instead of 2.25 and 1.2
instead of 1.5. Initial F1 multipliers are approximately 1.041333 and 1.07,
as documented in the original Low-LGN parameter notes. These are separate
from the base weights and initial LGN gains.

## Visual/LGN inputs

The saved kernel and spatial grid are reused without regeneration. `param.c`
is 0.001, `IB` 0.1, `tau_leak` 0.1, `lambda` 1, `dt` 0.1 ms,
`s_ex` 0.075, `v_th` 1.15, and speed 0.004 degrees/ms. The original helper
upsamples each 1 ms visual-drive sample tenfold before calling `NLIF`.

BG/grating `ex_gen2` uses high-cx initial gains E=1.15, I=1.20, or Low-LGN
initial gains E=1.11, I=1.15. Its original
history transitions are retained exactly: E switches at its 12-spike/240 ms
and 20-spike/400 ms tests to 1.12 and 1; I switches at its 7-spike/140 ms
and 12-spike/240 ms tests to 1.05 and 1. Read the code for exact indexing,
strict inequalities and one-way flags; these are not continuous decay laws.

Onset uses the separate `ex_gen_onset_test` rule. Before onset E/I gains are
1.11/1.15; each afferent's first spike strictly after onset starts its original
20-spike curve, and subsequent spikes after its 20th have gain 1. Background
history does not advance the onset curve. These settings match the historical
onset helper and numerical audit; they are deliberately exposed independently
of the steady-input initial gains in `model_config`.

The numerical and experimental literature sources for the model require an
author-approved bibliography. The release identifies the concrete local
parameter/source evidence without inventing missing paper references.

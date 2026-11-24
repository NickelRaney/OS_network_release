# Model Parameters and Network Snapshots

The repository ships two fixed network realizations representing the two V1
layer 4C$\alpha$ models:

| Preset | Alias | Data file | Complex LGN inputs | F1 scales (EE / EI) | BG/Grating initial LGN gain (E / I) | Onset initial LGN gain (E / I) |
| --- | --- | --- | --- | --- | --- | --- |
| `disordered_lgn` (default) | `high_cx` | `data/disordered_lgn.mat` | 2--5 (low), 7--8 (mid), 10--11 (high) | `2.25` / `1.5` | `1.15` / `1.20` | `1.11` / `1.15` |
| `low_lgn` | -- | `data/low_lgn.mat` | 0, 1, or 2 | `1.3` / `1.2` | `1.11` / `1.15` | `1.11` / `1.15` |

Each `.mat` file contains 11 runtime variables: `param`, `param_v1`, `init_v1`,
`kernel3`, `grid`, `Cat`, `Cat_i`, `LVwire`, `LVwire_i`, `inde`, and `indi`.
Use `S = load('data/disordered_lgn.mat')` to inspect them without polluting the
base workspace.

## Network Geometry and Connectivity

- **Cortical patch**: 4 hypercolumns covering $0.5^\circ \times 0.5^\circ$ of
  visual space, with periodic (torus) boundary conditions.
- **Populations**: 13,689 excitatory (E) neurons on a $117 \times 117$ grid,
  2,401 inhibitory (I) neurons on a $49 \times 49$ grid, and 110 LGN afferents
  (55 ON, 55 OFF) on a $10 \times 11$ grid (`grid`).
- **Central analysis ROI**: `inde(1,:)` selects 841 E neurons (594 simple, 247
  complex) and `indi(1,:)` selects 145 I neurons in the central hypercolumn,
  which prefers **vertical** gratings.
  - In `disordered_lgn`, the 247 central complex cells partition into 58 low
    (2--5 LGN), 86 medium (7--9 LGN), and 103 high (>9 LGN).
  - In `low_lgn`, the 247 central complex cells partition into 74 (0 LGN), 81
    (1 LGN), and 92 (2 LGN).
- **Connectivity convention**: Presynaptic row to postsynaptic column
  (`cmee`: E->E, `cmie`: E->I, `cmei`: I->E, `cmii`: I->I).
- **Network generation scripts**: Standalone scripts illustrating Section 2
  template and connectivity construction are provided under `scripts/`. By
  default, simulations load the fixed `.mat` snapshots above so runs are
  directly comparable across users.

## Base Synaptic Weights and Cortical Dynamics

Both snapshots store base synaptic weights in `param_v1`:

- `see = 0.0245`, `sie = 0.0176`, `sei = 0.051`, `sii = 0.09`
- `selgn = 0.0529`, `silgn = 0.09881`
- `sel6 = 0.0072`, `sil6 = 0.0042`, `segabab = 1.5`

Model-specific connectivity (`cmee`, `Cat`, `LVwire`) and lookup tables
(`l6cof`, `seetfm`, `sei_list`, `sincre_e`) are loaded from the corresponding
snapshot.

In `src/V1_4_fast.m` (integration step `dt = 0.1 ms`):

- **Membrane leak**: $\tau_{\text{leak}}^E = 20\text{ ms}$,
  $\tau_{\text{leak}}^I = 15\text{ ms}$.
- **Reversal potentials & threshold**: $V^E = 14/3$, $V^I = -2/3$, base spike
  threshold $V_{\text{thresh}} = 1$, reset $V_{\text{rest}} = 0$. After an E
  spike, the threshold increment is set to `0.5` and decays exponentially with
  $\tau_{\text{thresh}} = 5.05\text{ ms}$.
- **Synaptic time constants (rise / decay)**:
  - AMPA: `0.5 / 2 ms` (LGN, L4, L6; 80% of L4/L6 EPSC)
  - NMDA: `2 / 80 ms` (20% of L4/L6 EPSC)
  - GABA-A: `0.5 / 5 ms`
  - GABA-B: `3 / 100 ms` (triggered when an E cell receives >265 I spikes in 30 ms)
- **Short-term depression ($F_1$ curves)**:

```matlab
e_decay = ([1.2 1.2 1.2 1.20 1.19 1.16 1.12 1.07 1 0.96 0.94 0.92] - 0.92) * f1_ee_scale / 3 + 0.92;
i_decay = ([1.1 1.02 0.98 0.96 0.95 0.95 0.95 0.95 0.95 0.95 0.95 0.95] - 0.95) * f1_ei_scale / 1.5 + 0.95;
```

  With `disordered_lgn` (`2.25 / 1.5`), the initial multipliers are `1.13` and
  `1.10`. With `low_lgn` (`1.3 / 1.2`), they are `~1.0413` and `1.07`.
- **Ambient conductance**: Ornstein-Uhlenbeck process updated each step.
- **Layer 6 feedback**: Computed from a $21 \times 21$ Gaussian-weighted local
  E firing rate (89% simple, 11% complex) mapped through `L6tfm`, with
  activity-dependent facilitation.

## Visual Stimuli and LGN Inputs

- **Spatial & temporal kernels**: Loaded from `kernel3` (1 ms temporal sampling;
  16 ms OFF leading zeros and 25 ms ON leading zeros). Optional extra uniform
  delay can be added via `cfg.kernel_extra_delay_ms`.
- **LGN LIF parameters**: `c = 0.001`, `IB = 0.1`, `tau_leak = 0.1` (rate in
  $\text{ms}^{-1}$, i.e., 10 ms time constant), `lambda = 1`, `dt = 0.1 ms`,
  `s_ex = 0.075`, `v_th = 1.15`.
- **Drifting grating geometry**: Fixed at SF = 2.5 cycles/degree and speed = 4
  degrees/second (`param.speed = 0.004` deg/ms, TF = 10 Hz).
  - `cfg.direction = [1 0]` or `[-1 0]`: **vertical** grating drifting right or
    left (**optimal** orientation for the central ROI).
  - `cfg.direction = [0 1]` or `[0 -1]`: **horizontal** grating drifting up or
    down (**orthogonal** orientation for the central ROI).
- **LGN synaptic gain rules**:
  - `bg` and `grating` use `src/ex_gen2.m` (history-dependent step transitions
    from `lgn_gain_e_initial / lgn_gain_i_initial` to `1.12 / 1.05` and `1.0`).
  - `onset` uses `src/ex_gen_onset_test.m` (initial gains `1.11 / 1.15` prior
    to `onset_s`, followed by the 20-spike post-onset depression curve per
    afferent and gain `1.0` thereafter).

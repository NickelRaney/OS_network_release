# Network Construction Scripts

These standalone MATLAB functions document and regenerate the Section 2 network
layout and wiring rules:

- `build_lgn_grid.m`: Generates the $10 \times 11$ LGN lattice (`110` ON/OFF
  afferents) and cortical $117 \times 117$ E and $49 \times 49$ I grids over
  the $0.5^\circ \times 0.5^\circ$ patch (4 hypercolumns).
- `build_lgn_templates.m`: Constructs the oriented ON/OFF simple-cell LGN
  template library for $0^\circ$, $45^\circ$, $90^\circ$, and $135^\circ$
  orientation domains.
- `assign_lgn_connections.m`: Assigns LGN-to-L4 connections (`LVwire`, `Cat`,
  `LVwire_i`, `Cat_i`) for either `'disordered_lgn'` or `'low_lgn'`.
- `build_l4_connectivity.m`: Builds periodic-boundary (torus) truncated
  Gaussian L4 connectivity matrices (`cmee`, `cmie`, `cmei`, `cmii`) and the
  nearest-E mapping (`map_ei`).

## Example

```matlab
setup
rng(42, 'twister');
[grid, x, y] = build_lgn_grid();
templates = build_lgn_templates();
[LVwire, Cat, LVwire_i, Cat_i] = assign_lgn_connections(x, y, templates, 'disordered_lgn');
conn = build_l4_connectivity(x, y, Cat, 'disordered_lgn');
```

> **Note on reproducibility**: Standard simulations (`run_demo`,
> `run_simulation`, `examples/`) load the pre-generated fixed network snapshots
> in `data/disordered_lgn.mat` and `data/low_lgn.mat` by default so that all
> runs share the exact same cortical topology and calibrated lookup tables.

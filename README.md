# OS Network

MATLAB code for simulating macaque primary visual cortex (V1) responses to
inputs from the lateral geniculate nucleus (LGN).

The repository includes the **Disordered-LGN** (`high_cx`) and **Low-LGN**
(`low_lgn`) models, their parameter files, and examples for background
activity, drifting gratings, and stimulus onset.

## Requirements

MATLAB with **Statistics and Machine Learning Toolbox**. The data needed
to run the examples are included in this repository.

## Quick start

Download or clone this repository, then open MATLAB in its root folder.
Run a short demo of the Disordered-LGN model:

```matlab
setup
folder = run_demo('high_cx');
```

For the Low-LGN model, use:

```matlab
folder = run_demo('low_lgn');
```

Each run saves its results in a new folder under `outputs/`, including
simulation data (`result.mat`), a run summary (`run.json`), and a spike raster
(`raster.png`). The returned `folder` gives the location.

To check the installation, run `run_tests` after `setup`.

## Further information

- [Model configuration](config/model_config.m) and [parameter documentation](docs/PARAMETERS.md)
- [Output definitions](docs/OUTPUTS.md)
- [Validation](docs/VALIDATION.md) and [paper reproduction status](docs/PAPER_ALIGNMENT.md)
- [Code and data sources](docs/PROVENANCE.md)

Citation and licensing information is pending; see [LICENSE_STATUS.md](LICENSE_STATUS.md).

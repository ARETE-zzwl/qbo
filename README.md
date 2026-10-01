# QBO background experiments

[中文](README.zh-CN.md) · [Architecture](docs/architecture.md) · [Data and workflows](docs/workflows.md)

Input preparation and forcing diagnostics for QBO experiments in CESM/WACCM. The experiment holds the westerly and easterly QBO targets fixed across two SST and sea-ice backgrounds to study how the atmospheric response changes with the background climate.

## Experiment design

| | Early background, 1981–1994 | Late background, 1995–2024 |
|---|---|---|
| Westerly QBO | W × early | W × late |
| Easterly QBO | E × early | E × late |

QBO targets are built from MERRA-2 equatorial monthly zonal wind. Years from 1981–2024 are ranked by the March 70-hPa wind anomaly; the highest and lowest eight years form the W and E composites. Each target contains seven pressure levels and thirteen monthly nodes, from September before the selected March through September of that year. Both backgrounds use the same target file for each phase.

SST and sea-ice backgrounds are monthly HadISST climatologies for the two periods. They are mapped to the baseline FV grid by nearest-neighbour interpolation with periodic longitude. Missing source values are filled from the baseline.

## Installation

Python 3.10 or later:

```bash
git clone https://github.com/ARETE-zzwl/qbo.git
cd qbo
python -m venv .venv
```

Activate with `source .venv/bin/activate` on Linux/macOS or `.venv\Scripts\Activate.ps1` in PowerShell, then install from the repository root:

```bash
python -m pip install -e ".[test]"
```

## Try the workflow

```bash
python -m pip install -e ".[test,demo]"
qbo-demo --output outputs/demo
```

The demo generates synthetic wind, SST, sea-ice and tropopause data, runs both input builders, and writes target QC, boundary weights and an overview figure in PNG, SVG and PDF. Open `outputs/demo/README.md` to inspect the results. Each run uses a new output directory.

To diagnose your own snapshot:

```bash
qbo-diagnose --input data/snapshot.nc --output outputs/diagnostics
```

The [snapshot format](docs/workflows.md#diagnose-a-snapshot) uses explicit pressure fields on a complete longitude grid. The [research notes](docs/research-directions.md) describe the next model checks and scientific questions.

## Build inputs

Place the MERRA-2 monthly wind table, HadISST archives and CESM SST baseline under `data/`. File formats and processing steps are described in [Data and workflows](docs/workflows.md).

```bash
qbo-build-targets --input data/QBO_MERRA2-Uvals_00N_GSFC.txt --output outputs/targets
qbo-validate-targets outputs/targets --output outputs/targets/qc.json
qbo-build-backgrounds --sst-gz data/HadISST_sst.nc.gz --ice-gz data/HadISST_ice.nc.gz --baseline data/sst_baseline.nc --output outputs/backgrounds
```

Outputs are NetCDF files and JSON metadata. QBO files store `qbo(time, level)` to match WACCM's Fortran reader.

## Boundary weights

The column weight tapers smoothly to zero as the layer bottom approaches the most restrictive of three tropopause estimates. Taking the minimum across longitude gives a uniform weight for the latitude ring.

```python
import numpy as np
from qbo.boundary import safe_pressure, column_taper
from qbo.zonal import zonal_envelope

# Three tropopause diagnostics, four longitudes; pressure is in hPa.
pressure = np.array([[95., 100., 105., 110.]] * 3)
safe = safe_pressure(pressure, np.ones_like(pressure))
column_weights = column_taper(80., safe)
ring_weight = zonal_envelope(column_weights)
```

## Tests

```bash
python -m pytest
```

Tests cover input generation, NetCDF layout, longitude mapping, boundary weights, forcing diagnostics and CAM source preparation. Test inputs are generated locally.

[GitHub Actions](https://github.com/ARETE-zzwl/qbo/actions/workflows/tests.yml) runs the Python tests and demo on Linux and Windows with Python 3.10 and 3.13, plus the native Fortran/MPI checks on Linux. Reports and compiler logs are saved with each run.

With GNU Fortran and MPI installed:

```bash
python scripts/check_native.py --output outputs/native-check
```

This compares the Fortran kernel with the Python reference and checks zonal reductions and MPI column coverage. Use a new output directory for each run.

## Repository layout

```text
src/qbo/          Input builders, numerical operators and diagnostics
scripts/          CAM source preparation and archived-run analysis
native/           Fortran boundary, zonal, NetCDF and MPI contract tests
tests/            Python regression tests
docs/             Architecture, workflows and source provenance
vendor/cam/       Pinned CAM QBO source and upstream license
```

The v5 CAM driver integration is at the design stage; see [CAM integration](docs/architecture.md#cam-integration). To analyze existing model runs, point `QBO_WORKSPACE` at the experiment archive as described in [Data and workflows](docs/workflows.md#analyze-an-existing-experiment-archive).

`vendor/cam/` contains the CAM QBO source used by the adapter, with its [upstream license](vendor/cam/LICENSE.txt).

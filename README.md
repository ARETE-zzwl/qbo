# QBO experiments

[中文](README.zh-CN.md) · [Architecture](docs/architecture.md) · [Data and workflows](docs/workflows.md)

Tools for studying how the atmospheric response to the quasi-biennial oscillation (QBO) changes between early and late climate backgrounds in CESM/WACCM. The experiment pairs the same westerly and easterly QBO targets with two SST and sea-ice climatologies.

The repository contains target and background builders, tropopause-aware forcing operators, diagnostics, and Fortran tests for CAM integration.

## Experiment design

| | Early background, 1981–1994 | Late background, 1995–2024 |
|---|---|---|
| Westerly QBO | W × early | W × late |
| Easterly QBO | E × early | E × late |

Targets use MERRA-2 monthly zonal wind at the equator. Years are ranked by the March 70-hPa anomaly over 1981–2024; the upper and lower eight years define the W and E composites. Each target has seven pressure levels and thirteen monthly nodes, from the preceding September to the following September. The target file for a given phase is identical across backgrounds.

HadISST supplies the monthly SST and sea-ice climatologies. A periodic nearest-neighbour mapping places them on the baseline FV grid, with baseline values filling missing source cells.

## Install and test

Python 3.10 or later:

```bash
python -m venv .venv
```

Activate with `source .venv/bin/activate` on Linux/macOS or `.venv\Scripts\Activate.ps1` in PowerShell, then install from the repository root:

```bash
python -m pip install -e ".[test]"
python -m pytest
```

The Python tests build small arrays and NetCDF files locally. They cover input layout, periodic longitude mapping, boundary weights, forcing diagnostics, operator decomposition, and CAM source preparation.

## Build inputs

Place the MERRA-2 monthly wind table, HadISST archives and CESM SST baseline under `data/`. File formats and processing steps are described in [Data and workflows](docs/workflows.md).

```bash
qbo-build-targets --input data/QBO_MERRA2-Uvals_00N_GSFC.txt --output outputs/targets
qbo-validate-targets outputs/targets --output outputs/targets/qc.json
qbo-build-backgrounds --sst-gz data/HadISST_sst.nc.gz --ice-gz data/HadISST_ice.nc.gz --baseline data/sst_baseline.nc --output outputs/backgrounds
```

The builders write NetCDF files and JSON metadata with input and output hashes. QBO files store `qbo(time, level)`, matching the WACCM Fortran reader's `u_inp(level, time)` array.

## Use the operators

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

The column operator applies a smooth taper above the most restrictive valid tropopause. The zonal operator takes the minimum over longitude, giving a uniform weight that satisfies every column's boundary constraint.

## Repository layout

```text
src/qbo/          Input builders, numerical operators and diagnostics
scripts/          CAM source preparation and archived-run analysis
native/           Fortran boundary, zonal, NetCDF and MPI contract tests
tests/            Python regression tests
docs/             Architecture, workflows and source provenance
vendor/cam/       Pinned CAM QBO source and upstream license
```

Fortran checks use GNU Fortran and an MPI toolchain:

```bash
python scripts/check_native.py --output outputs/native-check
```

The command checks Python/Fortran boundary agreement, zonal reductions across MPI layouts, and twelve epoch/column-coverage cases.

## Development status

Input preparation, offline operators, the instrumented v4 CAM adapter, and isolated v5 MPI contracts have test coverage. The v5 prepare/reduce/resume integration remains a design described in the [architecture notes](docs/architecture.md#cam-integration). Full-model background-response experiments are the next stage.

Historical analysis scripts read a separate experiment archive through `QBO_WORKSPACE`. Large input datasets, model histories, job logs and machine-specific connection settings stay outside Git.

The pinned CAM source retains its [upstream license](vendor/cam/LICENSE.txt).

# QBO background experiments

[中文](README.zh-CN.md) · [Architecture](docs/architecture.md) · [Data formats](docs/workflows.md)

Input builders and forcing diagnostics for QBO experiments in CESM/WACCM.

The code builds westerly and easterly QBO targets from MERRA-2 equatorial winds, and early (1981–1994) and late (1995–2024) SST/sea-ice backgrounds from HadISST. Both backgrounds use the same pair of QBO targets. Python and Fortran routines calculate the column boundary weights and their minimum around each latitude ring.

## Install

Requires Python 3.10 or later.

```bash
git clone https://github.com/ARETE-zzwl/qbo.git
cd qbo
python -m venv .venv
```

Activate the environment with `source .venv/bin/activate` on Linux/macOS or `.venv\Scripts\Activate.ps1` in PowerShell, then install:

```bash
python -m pip install -e ".[test,demo]"
```

## Run the demo

```bash
qbo-demo --output outputs/demo
```

The demo generates synthetic inputs and runs the builders and boundary diagnostics. Results are collected in `outputs/demo/README.md`, with NetCDF data, CSV/JSON summaries and a figure in PNG, SVG and PDF. Use a new output directory for each run.

## Build experiment inputs

Place the source files under `data/` and run:

```bash
qbo-build-targets --input data/QBO_MERRA2-Uvals_00N_GSFC.txt --output outputs/targets
qbo-validate-targets outputs/targets --output outputs/targets/qc.json
qbo-build-backgrounds --sst-gz data/HadISST_sst.nc.gz --ice-gz data/HadISST_ice.nc.gz --baseline data/sst_baseline.nc --output outputs/backgrounds
```

QBO years are ranked by the March 70-hPa wind anomaly over 1981–2024. The highest and lowest eight years form the W/E composites, each with seven pressure levels and thirteen monthly nodes from the preceding September to September of the selected year. Target files use `qbo(time, level)` for the WACCM Fortran reader.

HadISST monthly climatologies are mapped to the baseline FV grid with periodic nearest-neighbour interpolation. Missing cells use the baseline values. Required variables and output files are listed in [Data formats](docs/workflows.md).

## Diagnose a snapshot

```bash
qbo-diagnose --input data/snapshot.nc --output outputs/diagnostics
```

This writes column and zonal weights to `weights.nc`, with CSV and JSON summaries. See the [snapshot format](docs/workflows.md#diagnose-a-snapshot) for pressure variables and grid requirements.

## Test

```bash
python -m pytest
```

With GNU Fortran and MPI installed:

```bash
python scripts/check_native.py --output outputs/native-check
```

[CI](https://github.com/ARETE-zzwl/qbo/actions/workflows/tests.yml) runs Python tests and the demo on Linux and Windows, plus native Fortran/MPI checks on Linux.

## Layout

```text
src/qbo/       Input builders, operators and diagnostics
scripts/       CAM source preparation and archive analysis
native/        Fortran and MPI kernels and tests
tests/         Python tests
docs/          Architecture and usage
vendor/cam/    CAM source and upstream license
```

CAM adapter setup and driver integration are documented in [Architecture](docs/architecture.md#cam-integration). The vendored CAM code retains its [upstream license](vendor/cam/LICENSE.txt).

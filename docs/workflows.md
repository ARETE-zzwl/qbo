# Data and workflows

[中文](workflows.zh-CN.md) · [README](../README.md)

Run commands from the repository root after installing `.[test]`.

## Input files

| Input | Required content |
|---|---|
| `QBO_MERRA2-Uvals_00N_GSFC.txt` | ASCII monthly table with a `P (hPa):` header and `YYYYMM` rows; selected levels are 10, 20, 30, 40, 50, 70 and 100 hPa |
| `HadISST_sst.nc.gz` | `sst(time, latitude, longitude)` and a CF time coordinate |
| `HadISST_ice.nc.gz` | `sic(time, latitude, longitude)` on the same grid, with sea-ice fractions in 0–1 |
| `sst_baseline.nc` | `lat`, `lon`, `time`, `date`, `datesec`, `SST_cpl`, `ice_cov`; the two fields need `long_name` and `units` attributes |

The QBO selection period is 1981–2024. Supply the preceding September–December as well, since composites span the September before each selected March. HadISST needs every month from January 1981 through December 2024. The SST baseline supplies twelve monthly records on the target grid.

The wind parser rejects duplicate/invalid months, non-finite winds and missing or duplicate target pressure levels. HadISST records are sorted by calendar year/month; each requested month must occur once. NetCDF missing-value masks are retained when forming monthly means. Cells missing throughout the period are filled from the baseline after regridding.

```bash
qbo-build-targets --input data/QBO_MERRA2-Uvals_00N_GSFC.txt --output outputs/targets
qbo-validate-targets outputs/targets --output outputs/targets/qc.json
qbo-build-backgrounds --sst-gz data/HadISST_sst.nc.gz --ice-gz data/HadISST_ice.nc.gz --baseline data/sst_baseline.nc --output outputs/backgrounds
```

The target builder writes `qbo_target_W.nc`, `qbo_target_E.nc`, `selection.csv` and `target_metadata.json`. NetCDF dates are month-start nodes in model years 1 and 2, with September repeated as the endpoint. A March node therefore represents an interpolation node rather than a calendar-month integral.

The background builder writes one early and one late climatology plus `background_metadata.json`. It clips ice fractions to 0–1 and writes equal initial values to the main and `_prediddle` fields. The CESM monthly-mean-preserving adjustment is a subsequent model-input processing step; it is not performed by this builder.

## Diagnose a snapshot

`qbo-diagnose` reads one snapshot using these NetCDF names and dimensions:

| Variable | Dimensions | Meaning |
|---|---|---|
| `lat`, `lon` | `(lat)`, `(lon)` | Degrees; equally spaced latitude centers excluding the poles, and a complete equally spaced longitude ring without a duplicate seam |
| `tropopause_hpa` | `(diagnostic, lat, lon)` | TROP, TROPP and TROPF pressures, in that order; diagnostic dimension length 3 |
| `found` | `(diagnostic, lat, lon)` | Corresponding found flags; values greater than 0.5 mean found |
| `layer_bottom_hpa` | `(level)` or `(level, lat, lon)` | Positive layer-bottom pressures in hPa |

The optional global attribute `input_kind` is copied into outputs; the demo uses `synthetic`. For CAM data, select one time, reconstruct layer-bottom pressures from the hybrid interfaces and surface pressure, convert Pa to hPa, and retain every longitude in each selected latitude ring. Regridding before taking the minimum changes the operator being evaluated.

```bash
qbo-diagnose --input data/snapshot.nc --output outputs/diagnostics
```

The new output directory contains `weights.nc`, `summary.csv` and `report.json`. Summaries use cosine latitude weights. Invalid tropopause columns receive zero weight. These are geometric support statistics; the actual relaxation rate also depends on the model's latitude taper, relaxation strength and other forcing conditions.

`qbo-demo --output outputs/demo` creates a small example of this format, input-builder outputs and a three-panel figure. Install `.[demo]` to draw the figure. The CLI itself uses the base package dependencies.

## Native tests

On Linux with `gfortran`, `mpifort` and `mpiexec` on `PATH`:

```bash
python scripts/check_native.py --output outputs/native-check
```

Use a new output directory for each run. The runner saves compiler/run logs, generated vectors, numerical comparisons and a JSON summary. It checks:

1. The Fortran column kernel against 10,047 Python vectors, including boundary and invalid-input cases.
2. A small zonal fixture with one/two ranks and cyclic, block and single-owner partitions, plus an incomplete-coverage case.
3. The twelve two-rank contract cases defined in `native/contract/cases.json`.

`native/netcdf/read_qbo_native.f90` is a separate reader for environments with the NetCDF Fortran library. The Python tests also verify the target layout through the NetCDF C API.

## CAM source preparation

```bash
python scripts/prepare_cam_adapter.py outputs/cam-adapter
```

This creates a new directory containing `SourceMods/src.cam/`, native harness sources, `manifest.json` and `SHA256SUMS`. The source baseline is `vendor/cam/qbo.F90`; the adapter targets the archived CAM `cam_cesm2_1_rel_60` source. Preparing files is local. Building and running a CESM case uses the site's CESM toolchain and configuration.

## Analyze an existing experiment archive

The analysis scripts retain the September 2026 archive layout and its input hashes. Set `QBO_WORKSPACE` to the directory containing the original `outputs/`, `scripts/`, protocols and `source_paths.json`:

```bash
export QBO_WORKSPACE=/path/to/qbo-archive
python scripts/analysis/evaluate_boundary_v4.py --output outputs/column-assessment.json
python scripts/analysis/audit_qbo_operator_scope.py --output outputs/operator-assessment.json
python scripts/analysis/assess_zonal_v5.py --output outputs/zonal-assessment.json
```

In PowerShell, set `$env:QBO_WORKSPACE = 'D:/path/to/qbo-archive'`. The local `source_paths.json` identifies the original input experiment and HadISST directory. It belongs to the archive and is not needed for installation, input builders or regression tests.

| Script | Archive content it reads |
|---|---|
| `audit_inputs.py` | Original/corrected targets and backgrounds; pinned CAM source archive |
| `audit_remote_qbo.py` | Native reader output, vertical grid and job checkpoint |
| `audit_qbo_technical.py` | Short-run history package and coupler termination log |
| `audit_v3_preconditions.py` | Case configuration, namelist and build evidence |
| `screen_qbo_boundaries.py` | History packages for hard-boundary comparisons |
| `evaluate_boundary_v4.py` | Archived states and the v3 assessment |
| `audit_qbo_operator_scope.py` | v4 assessment and source fingerprints |
| `assess_zonal_v5.py` | v4 assessment and the same six archived states |
| `summarize_remote_increment.py` | Runtime reports and Slurm accounting |

Scripts with an `--output` argument support a new destination. The input, native-read, technical-run and summary audits write reports at their established archive locations. Use a copy of the archive when retaining those reports unchanged matters.

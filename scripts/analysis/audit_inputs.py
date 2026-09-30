"""Independently compare v1/v2 inputs and reconstruct QBO composites."""
from __future__ import annotations

import hashlib
import json
import os
import tarfile
import xml.etree.ElementTree as ET
from pathlib import Path

import netCDF4
import numpy as np


ROOT = Path(os.environ.get("QBO_WORKSPACE", Path(__file__).resolve().parents[2]))
CONFIG = json.loads((ROOT / "source_paths.json").read_text())
SOURCE = Path(CONFIG["experiment"])
OUT = ROOT / "outputs/20260913_input_audit"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_arrays(path):
    with netCDF4.Dataset("memory", memory=path.read_bytes()) as dataset:
        arrays = {name: np.asarray(var[:]) for name, var in dataset.variables.items()}
        dimensions = {name: var.dimensions for name, var in dataset.variables.items()}
    return arrays, dimensions


def audit_qbo():
    raw_path = SOURCE / "qbo_source/QBO_MERRA2-Uvals_00N_GSFC.txt"
    header = raw_path.read_text(encoding="ascii").splitlines()
    pressure = np.array([float(x) for x in header[13].split(":")[1].split()])
    records = np.loadtxt(raw_path, skiprows=15)
    stamps = records[:, 0].astype(int)
    lookup = {stamp: i for i, stamp in enumerate(stamps)}
    if len(lookup) != len(stamps):
        raise ValueError("Duplicate source month")
    levels = [10, 20, 30, 40, 50, 70, 100]
    columns = [int(np.flatnonzero(pressure == p)[0]) + 1 for p in levels]
    years = np.arange(1981, 2025)
    march70 = np.array([records[lookup[y * 100 + 3], columns[5]] for y in years])
    order = np.argsort(march70, kind="stable")
    selections = {"W": years[order[-8:][::-1]], "E": years[order[:8]]}
    old_metadata = json.loads((SOURCE / "qbo_targets_v1/target_metadata.json").read_text())
    new_metadata = json.loads((OUT / "qbo_targets_v2/target_metadata.json").read_text())
    assert old_metadata["input_sha256"] == new_metadata["input_sha256"] == sha256(raw_path)
    assert old_metadata["source_years"] == new_metadata["source_years"]
    results = {}
    for phase, selected in selections.items():
        profiles = []
        for year in selected:
            # Independent year/month arithmetic, not the builder's date_sequence.
            month_indices = (int(year) - 1) * 12 + 8 + np.arange(13)
            dates = (month_indices // 12) * 100 + month_indices % 12 + 1
            profiles.append(records[[lookup[int(d)] for d in dates]][:, columns])
        expected = np.mean(profiles, axis=0)
        filename = f"qbo_target_{phase}.nc"
        old_path = SOURCE / "qbo_targets_v1" / filename
        new_path = OUT / "qbo_targets_v2" / filename
        old, old_dims = read_arrays(old_path)
        new, new_dims = read_arrays(new_path)
        np.testing.assert_array_equal(old["qbo"].T, new["qbo"])
        np.testing.assert_allclose(new["qbo"], expected, rtol=0, atol=1e-12)
        for coordinate in ["date", "secs", "level", "time"]:
            np.testing.assert_array_equal(old[coordinate], new[coordinate])
        assert new_dims["qbo"] == ("time", "level")
        results[phase] = {
            "source_years": selected.tolist(),
            "v1_dimensions": old_dims["qbo"],
            "v2_dimensions": new_dims["qbo"],
            "max_value_change_mps": float(np.max(np.abs(old["qbo"].T - new["qbo"]))),
            "independent_composite_max_error_mps": float(np.max(np.abs(expected - new["qbo"]))),
            "march_1_70hpa_mps": float(new["qbo"][6, 5]),
            "march_continuous_linear_mean_70hpa_mps": float(new["qbo"][6:8, 5].mean()),
            "v1_sha256": sha256(old_path),
            "v2_sha256": sha256(new_path),
        }
    results["source_sha256"] = sha256(raw_path)
    results["same_qbo_model_fidelity_verified"] = False
    results["source_definition"] = "MERRA-2 00N absolute monthly wind, not ERA5 5S-5N"
    results["date_convention"] = "Month-start nodes; March node is not March mean"
    return results


def audit_backgrounds():
    baseline, _ = read_arrays(SOURCE / "background_source/sst_baseline.nc")
    metadata = json.loads((OUT / "backgrounds_v2/background_metadata.json").read_text())
    old_metadata = json.loads((SOURCE / "backgrounds_v1/background_metadata.json").read_text())
    assert metadata["source_hashes"] == old_metadata["source_hashes"]
    results = {}
    for period in ["early", "late"]:
        filename = f"sst_HadISST_{period}_0.9x1.25_clim.nc"
        old_path = SOURCE / "backgrounds_v1" / filename
        new_path = OUT / "backgrounds_v2" / filename
        old, _ = read_arrays(old_path)
        new, _ = read_arrays(new_path)
        for coordinate in ["lat", "lon", "time", "date", "datesec"]:
            np.testing.assert_array_equal(old[coordinate], new[coordinate])
        seam = np.flatnonzero(new["lon"] == 0)
        assert seam.size == 1
        nonseam = new["lon"] != 0
        result = {}
        for variable in ["SST_cpl", "ice_cov"]:
            difference = new[variable] - old[variable]
            np.testing.assert_array_equal(new[variable][:, :, nonseam], old[variable][:, :, nonseam])
            assert np.isfinite(new[variable]).all()
            result[variable] = {
                "changed_month_grid_cells": int(np.count_nonzero(difference)),
                "max_abs_change": float(np.max(np.abs(difference))),
                "changes_outside_greenwich": int(np.count_nonzero(difference[:, :, nonseam])),
                "finite": True,
            }
        result["v1_greenwich_sst_equals_fixed_baseline_cells"] = int(np.count_nonzero(
            old["SST_cpl"][:, :, seam] == baseline["SST_cpl"][:, :, seam]
        ))
        result["greenwich_month_grid_cells"] = int(old["SST_cpl"][:, :, seam].size)
        result["fallback_and_clipping"] = {
            key: value for key, value in metadata["outputs"][period].items()
            if key.endswith("_filled") or key.endswith("_clipped")
        }
        assert result["fallback_and_clipping"]["sst_missing_filled"] > 0
        result["monthly_mean_preserving_adjustment_applied"] = not np.array_equal(
            new["SST_cpl"], new["SST_cpl_prediddle"]
        )
        result["v1_sha256"], result["v2_sha256"] = sha256(old_path), sha256(new_path)
        results[period] = result
    results["ocean_mask_and_temporal_interpolation_qc_complete"] = False
    results["complete_external_background_ready"] = False
    return results


def audit_source_contracts():
    archive = SOURCE / "cam_cesm2_1_rel_60.tar.gz"
    results = {"archive_sha256": sha256(archive), "files": {}}
    with tarfile.open(archive) as source:
        for relative in [
            "src/physics/waccm/qbo.F90",
            "bld/namelist_files/use_cases/waccm_tsmlt_hist_cam6.xml",
            "bld/namelist_files/use_cases/waccm_tsmlt_ssp245_cam6.xml",
        ]:
            content = source.extractfile(relative).read()
            text = content.decode("utf-8", errors="strict")
            item = {"sha256": hashlib.sha256(content).hexdigest()}
            if relative.endswith(".xml"):
                tree = ET.fromstring(text)
                tags = ["flbc_file", "solar_irrad_data_file", "tgcm_ubc_file", "ext_frc_type", "srf_emis_type"]
                item["settings"] = {tag: tree.findtext(tag) for tag in tags}
            else:
                tokens = ["allocate( u_inp(levsiz,timesiz)", "nf90_get_var( ncid, uqboid, u_inp )",
                          "tauz(kbot+1)", "tconst = 10._r8", "yr*365._r8"]
                item["contract_lines"] = {
                    token: [i for i, line in enumerate(text.splitlines(), 1) if token in line]
                    for token in tokens
                }
                assert all(item["contract_lines"].values())
            results["files"][relative] = item
    return results


def main():
    # These are archived short technical runs, not a production benchmark.
    runs = [("32406", 2, 4 * 3600 + 17 * 60 + 20), ("32415", 1, 2 * 3600 + 60 + 14)]
    costs = [{
        "archived_job_id": job, "model_days": days, "wall_seconds": seconds, "mpi_tasks": 16,
        "effective_sypd": days / 365 / (seconds / 86400),
        "linear_extrapolation_core_hours_per_model_year": 16 * seconds / 3600 / days * 365,
        "production_benchmark": False,
    } for job, days, seconds in runs]
    payload = {
        "status": "local_input_regressions_verified_not_scientific_pilot",
        "qbo": audit_qbo(),
        "backgrounds": audit_backgrounds(),
        "pinned_source_contracts": audit_source_contracts(),
        "technical_cost_screen": costs,
        "software": {"numpy": np.__version__, "netCDF4": netCDF4.__version__},
    }
    scripts = [
        SOURCE / name for name in CONFIG["original_script_sha256"]
    ] + [Path(__file__), ROOT / "tests/test_input_contracts.py"]
    payload["analysis_script_sha256"] = {str(path): sha256(path) for path in scripts}
    (OUT / "independent_audit.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

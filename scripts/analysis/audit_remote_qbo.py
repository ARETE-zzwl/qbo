"""Compare native WACCM input reads with the target NetCDF files."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import tarfile

import netCDF4
import numpy as np

from qbo.forcing import forcing_layers, check_native_values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(os.environ.get("QBO_WORKSPACE", Path(__file__).resolve().parents[2]))
    job = args.job_dir
    checkpoint = json.loads((job / "checkpoint.json").read_text())
    job_id = checkpoint["slurm_job_id"]
    if "NATIVE_QC_COMPLETED" not in (job.parent / f"qc-{job_id}.out").read_text():
        raise ValueError("Native job completion marker is absent.")
    native = {}
    for phase in ("W", "E"):
        target = root / f"outputs/20260913_input_audit/qbo_targets_v2/qbo_target_{phase}.nc"
        with netCDF4.Dataset(target) as data:
            pressure = data["level"][:]
            dates = data["date"][:]
            native[phase] = check_native_values(
                np.loadtxt(job / f"native_{phase}.txt"), pressure, dates,
                data["secs"][:], data["qbo"][:],
            )
        if not re.search(r"NETCDF_STATUS\s+-57\b", (job / f"old_{phase}.txt").read_text()):
            raise ValueError("Archived v1 negative control did not return NC_EEDGE.")
    source_root = Path(json.loads((root / "source_paths.json").read_text())["experiment"])
    source_checks = {}
    with tarfile.open(source_root / "cam_cesm2_1_rel_60.tar.gz") as archive:
        for remote, actual_hash in checkpoint["source_sha256"].items():
            member = "src/" + remote.split("/cam/src/", 1)[1]
            raw = archive.extractfile(member).read()
            normalized = raw.replace(b"\r\n", b"\n")
            source_checks[member] = {
                "pinned_raw_sha256": hashlib.sha256(raw).hexdigest(),
                "pinned_lf_sha256": hashlib.sha256(normalized).hexdigest(),
                "remote_sha256": actual_hash,
                "equivalent_after_lf_normalization": hashlib.sha256(normalized).hexdigest() == actual_hash,
            }
    if not all(row["equivalent_after_lf_normalization"] for row in source_checks.values()):
        raise ValueError("The deployed source differs from the pinned source beyond line endings.")
    grid = np.loadtxt(job / "vertical_grid.txt")
    layers = forcing_layers(grid[:, 1], pressure)
    for row in layers:
        level = grid[row["index"] - 1]
        row.update(top_interface_hpa=float(level[2]),
                   bottom_interface_hpa=float(level[3]), hybm=float(level[4]))
    result = {
        "scope": "native_input_QC_not_scientific_pilot",
        "job_id": job_id, "native_input_qc_passed": True,
        "native_values": native, "v1_negative_controls_NC_EEDGE": True,
        "source_equivalence": source_checks, "model_levels": len(grid),
        "forcing_layers": layers,
        "lowest_forced_midpoint_hpa": layers[-1]["pressure_hpa"],
        "lowest_forced_bottom_interface_hpa": layers[-1]["bottom_interface_hpa"],
        "all_forced_layers_are_pure_pressure": all(row["hybm"] == 0 for row in layers),
        "target_dates": [int(dates[0]), int(dates[-1])],
        "parent_case_start_date": checkpoint["case"]["RUN_STARTDATE"],
        "parent_case_requires_new_start_and_continue_false": True,
        "physical_tropopause_boundary_verified": False,
        "same_qbo_runtime_fidelity_verified": False,
        "personal_quota_verified": checkpoint["personal_quota"]["verified"],
    }
    output = job / "native_assessment.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: result[key] for key in (
        "job_id", "native_input_qc_passed", "native_values",
        "lowest_forced_midpoint_hpa", "lowest_forced_bottom_interface_hpa",
        "all_forced_layers_are_pure_pressure",
    )}, indent=2))


if __name__ == "__main__":
    main()

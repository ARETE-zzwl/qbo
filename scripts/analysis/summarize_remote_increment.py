"""Summarize observed runtime, diagnostics and apparent file sizes."""
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(os.environ.get("QBO_WORKSPACE", Path(__file__).resolve().parents[2]))
OUT = ROOT / "outputs/20260913_remote_resume"


def main():
    accounting = json.loads((OUT / "slurm_accounting.json").read_text())
    jobs = {row["JobID"]: row for row in accounting["records"]}
    cases = {}
    for phase, job_id, package in (
        ("W", "33240", "history_W_33242"), ("E", "33243", "history_E_33244")
    ):
        directory = OUT / package
        report = json.loads((directory / "runtime_assessment.json").read_text())
        sizes = {
            name: int(size) for name, size in (
                line.rsplit(" ", 1) for line in (directory / "run_file_sizes.txt").read_text().splitlines()
            )
        }
        integrated = [row for row in report["records"] if row["nstep"] > 0]
        assert jobs[job_id]["State"] == "COMPLETED" and jobs[job_id]["ExitCode"] == "0:0"
        cases[phase] = {
            "job_id": job_id, "wall_seconds": int(jobs[job_id]["ElapsedRaw"]),
            "allocated_cpus": int(jobs[job_id]["AllocCPUS"]),
            "mpi_step_maxrss_reported": jobs[job_id + ".0"]["MaxRSS"],
            "integrated_steps_interface_passed": report["integrated_steps_interface_passed"],
            "all_record_max_qbo_u0_error_mps": max(
                row["qbo_u0_max_abs_error_mps"] for row in report["records"]
            ),
            "primary_direct_nudging_area_fraction": [
                row["tropopause"]["combined"]["nonzero_tendency_below_fraction"]
                for row in integrated
            ],
            "cold_point_direct_nudging_area_fraction": [
                row["tropopause"]["cold_point"]["nonzero_tendency_below_fraction"]
                for row in integrated
            ],
            "run_regular_file_bytes": sum(sizes.values()),
            "cam_main_restart_bytes": sum(
                size for name, size in sizes.items() if ".cam.r." in name
            ),
            "compressed_h0_bytes": sum(item["bytes"] for item in report["files"]),
        }
    result = {
        "scope": "bounded_technical_increment_complete_scientific_gate_failed",
        "cases": cases,
        "total_model_job_core_hours": sum(
            case["wall_seconds"] * case["allocated_cpus"] / 3600 for case in cases.values()
        ),
        "total_run_regular_file_bytes": sum(case["run_regular_file_bytes"] for case in cases.values()),
        "storage_definition": "Apparent sizes of regular files in both run directories; excludes staged reference symlinks.",
        "rss_definition": "Slurm MPI-step MaxRSS, not an independently sampled sum of rank RSS.",
        "scientific_pilot_permitted": False, "seasonal_fidelity_verified": False,
        "scripts_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted((ROOT / "scripts").glob("*.py"))
        },
    }
    (OUT / "increment_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

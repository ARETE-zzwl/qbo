"""Audit archived v3 settings and screen zero-buffer geometry, without runs."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shlex
import sys
import xml.etree.ElementTree as ET

import netCDF4
import numpy as np

from qbo.preconditions import namelist_failures, hard_cutoff_strength, cutoff_overlap

from qbo.forcing import PRESSURE_TOLERANCE_HPA

REMOTE = "/public/home/zhushibang/qbo_factorial_20260910"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Use a new output path; archived results are immutable.")
    root = Path(os.environ.get("QBO_WORKSPACE", Path(__file__).resolve().parents[2]))
    evidence = args.evidence
    manifest = json.loads((evidence / "manifest.json").read_text(encoding="utf-8"))
    for item in manifest["files"]:
        if hashlib.sha256((evidence / item["local"]).read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError("Evidence hash changed: " + item["local"])

    source = Path(json.loads((root / "source_paths.json").read_text())["experiment"])
    cime = source / "cesm_source/cime/scripts/lib"
    sys.path.insert(0, str(cime))
    from CIME.namelist import parse

    cases = {}
    for phase in ("W", "E"):
        directory = evidence / phase
        nml = parse(in_file=str(directory / "atm_in"))
        names = (
            "qbo_use_forcing", "qbo_cyclic", "qbo_forcing_file",
            "nhtfrq", "mfilt", "avgflag_pertape", "empty_htapes", "fincl1",
        )
        values = {name: nml.get_value(name) for name in names}
        settings = {
            entry.attrib["id"]: entry.attrib.get("value")
            for filename in ("env_run.xml", "env_build.xml")
            for entry in ET.parse(directory / filename).iter("entry")
        }
        target = root / f"outputs/20260913_input_audit/qbo_targets_v2/qbo_target_{phase}.nc"
        with netCDF4.Dataset(target) as data:
            dates = [int(data["date"][0]), int(data["date"][-1])]
        start = int(settings["RUN_STARTDATE"].replace("-", ""))
        failures = namelist_failures(values, REMOTE + f"/inputs/qbo_targets_v2/qbo_target_{phase}.nc")
        cases[phase] = {
            "namelist_failures": failures, "namelist_gate_passed": not failures,
            "runtime_values": values,
            "settings": {key: settings[key] for key in (
                "EXEROOT", "RUN_STARTDATE", "STOP_OPTION", "STOP_N",
                "CONTINUE_RUN", "RESUBMIT", "CALENDAR",
            )},
            "candidate_target_dates": dates,
            "candidate_target_covers_start": dates[0] <= start < dates[1],
            "candidate_target_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
            "qbo_open_marker_found": any(
                "qbo_init: successfully opened" in gzip.open(path, "rt").read()
                for path in directory.glob("atm.log.*.gz")
            ),
            "coupler_success_marker_found": any(
                "SUCCESSFUL TERMINATION" in gzip.open(path, "rt").read()
                for path in directory.glob("cpl.log.*.gz")
            ),
            "boundary_safety": "unverified_forcing_disabled",
            "qbo_fidelity": "unverified_forcing_disabled",
            "forced_numerical_stability": "unverified_forcing_disabled",
        }

    with gzip.open(evidence / "build/atm.bldlog.260914-102555.gz", "rt") as handle:
        compile_commands = [line.strip() for line in handle
                            if line.startswith("mpiifort ") and "qbo.F90" in line]
    compiled_sources = sorted({token for line in compile_commands for token in shlex.split(line)
                               if token.endswith("/qbo.F90")})
    if not compiled_sources:
        raise ValueError("No QBO compilation command found.")

    records, inputs = [], {}
    packages = root / "outputs/20260913_remote_resume"
    for package in ("history_W_33242", "history_E_33244"):
        for path in sorted((packages / package).glob("*.nc")):
            inputs[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
            with netCDF4.Dataset(path) as data:
                p0 = float(data["P0"][()]) / 100
                mid = (data["hyam"][:] + data["hybm"][:]) * p0
                edge = (data["hyai"][:] + data["hybi"][:]) * p0
                active = hard_cutoff_strength(mid, [10, 100]) > 0
                if np.any(data["hybm"][:][active] != 0) or data["hybi"][np.flatnonzero(active)[-1] + 1] != 0:
                    raise ValueError("Reference pressures cannot screen terrain-dependent layers.")
                weights = data["gw"][:] * (np.abs(np.deg2rad(data["lat"][:])) <= .384)
                for index, step in enumerate(data["nsteph"][:]):
                    for definition in ("TROP", "TROPP", "TROPF"):
                        pressure = np.ma.asarray(data[definition + "_P"][index])
                        found = np.ma.asarray(data[definition + "_FD"][index])
                        valid = (~np.ma.getmaskarray(pressure) & ~np.ma.getmaskarray(found)
                                 & (found.filled(0) > .5))
                        records.append({
                            "package": package, "file": path.name, "nstep": int(step),
                            "tropopause_definition": definition,
                            **cutoff_overlap(mid, edge, pressure.filled(np.nan) / 100, valid, weights),
                        })
    summary = {}
    for definition in ("TROP", "TROPP", "TROPF"):
        rows = [row for row in records if row["tropopause_definition"] == definition]
        if len(rows) != 6:
            raise ValueError("Expected all six archived records, including both startup records.")
        summary[definition] = {
            "records": len(rows),
            "minimum_valid_area_fraction": min(row["valid_area_fraction"] for row in rows),
            **{key: max((row[key] for row in rows if row[key] is not None), default=None)
               for key in ("midpoint_below_fraction", "interface_below_fraction")},
        }
    result = {
        "scope": "v3_configuration_and_offline_hard_cutoff_screen_only",
        "cases": cases, "compiled_qbo_sources": compiled_sources,
        "v3_override_compiled": any("/SourceMods/" in path for path in compiled_sources),
        "hard_cutoff_screen": {
            "target_hpa": [10, 100], "pressure_tolerance_hpa": PRESSURE_TOLERANCE_HPA,
            "summary": summary, "records": records,
            "not_actual_v3_tendencies": True, "not_a_seasonal_safety_test": True,
        },
        "history_sha256": inputs,
        "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cime_namelist_parser_sha256": hashlib.sha256((cime / "CIME/namelist.py").read_bytes()).hexdigest(),
        "scientific_pilot_approved": False, "new_model_jobs_submitted": [],
    }
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({
        "namelist_failures": {phase: row["namelist_failures"] for phase, row in cases.items()},
        "compiled_qbo_sources": compiled_sources, "hard_cutoff_screen": summary,
    }, indent=2))


if __name__ == "__main__":
    main()

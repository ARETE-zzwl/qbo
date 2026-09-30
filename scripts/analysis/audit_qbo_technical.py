"""Reconstruct QBO diagnostics and quantify tropopause overlap in short runs."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path

import cftime
import netCDF4
import numpy as np

from qbo.forcing import (
    PRESSURE_TOLERANCE_HPA, expected_qbo_diagnostic, forcing_layers, tropopause_overlap,
)

TOLERANCE_MPS = 2e-5
DYNAMIC_FIELDS = ("U", "T", "PS", "QBOTEND", "QBO_U0")


def model_dates(dates, seconds):
    return [
        cftime.DatetimeNoLeap(
            int(date) // 10000, (int(date) // 100) % 100, int(date) % 100,
            int(sec) // 3600, (int(sec) // 60) % 60, int(sec) % 60,
        )
        for date, sec in zip(dates, seconds)
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--package", type=Path, required=True)
    parser.add_argument("--phase", choices=("W", "E"), required=True)
    args = parser.parse_args()
    root = Path(os.environ.get("QBO_WORKSPACE", Path(__file__).resolve().parents[2]))
    target = root / ("outputs/20260913_input_audit/qbo_targets_v2/qbo_target_" + args.phase + ".nc")
    with netCDF4.Dataset(target) as source:
        target_pressure = np.asarray(source["level"][:])
        target_wind = np.asarray(source["qbo"][:])
        dates = model_dates(source["date"][:], source["secs"][:])
    logpath = next(args.package.glob("cpl.log.*"))
    log = logpath.read_bytes()
    if logpath.suffix == ".gz":
        log = gzip.decompress(log)
    if b"SUCCESSFUL TERMINATION" not in log:
        raise ValueError("Actual coupler termination marker is absent.")
    expected_hashes = {}
    for line in (args.package / "history.sha256").read_text().splitlines():
        digest, filename = line.split(None, 1)
        if Path(filename).parent.name == args.package.name:
            expected_hashes[Path(filename).name] = digest
    records, files = [], []
    for path in sorted(args.package.glob("*.nc")):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != expected_hashes[path.name]:
            raise ValueError("Compressed history transfer checksum mismatch.")
        files.append({"name": path.name, "sha256": digest, "bytes": path.stat().st_size})
        with netCDF4.Dataset(path) as data:
            assert data["QBO_U0"].dimensions == ("time", "lev", "lat", "lon")
            assert data["time"].calendar in ("noleap", "365_day")
            lat = np.asarray(data["lat"][:])
            mid = (data["hyam"][:] + data["hybm"][:]) * float(data["P0"][()]) / 100
            edges = (data["hyai"][:] + data["hybi"][:]) * float(data["P0"][()]) / 100
            layers = forcing_layers(mid, target_pressure)
            indices = [row["index"] - 1 for row in layers]
            assert np.all(data["hybm"][indices] == 0)
            last = indices[-1]
            latitude_band = np.abs(np.deg2rad(lat)) <= .384
            weights = np.asarray(data["gw"][:]) * latitude_band
            support = np.zeros((len(mid), len(lat)), dtype=bool)
            support[indices, :] = latitude_band
            nodes = netCDF4.date2num(dates, data["time"].units, calendar="noleap")
            for index, point in enumerate(data["time"][:]):
                assert nodes[0] <= point <= nodes[-1]
                profile = [
                    np.interp(point, nodes, target_wind[:, level])
                    for level in range(len(target_pressure))
                ]
                expected = expected_qbo_diagnostic(mid, lat, target_pressure, profile)
                values, ranges = {}, {}
                for name in DYNAMIC_FIELDS:
                    array = np.ma.asarray(data[name][index]).filled(np.nan)
                    values[name] = array
                    good = np.isfinite(array)
                    ranges[name] = {
                        "finite": bool(good.all()),
                        "minimum": float(array[good].min()) if good.any() else None,
                        "maximum": float(array[good].max()) if good.any() else None,
                        "units": getattr(data[name], "units", ""),
                    }
                finite = all(value["finite"] for value in ranges.values())
                error = float(np.max(np.abs(values["QBO_U0"] - expected[:, :, None]))) if finite else None
                outside = values["QBOTEND"][~support, :]
                outside_max = float(np.max(np.abs(outside))) if finite else None
                record = {
                    "nstep": int(data["nsteph"][index]),
                    "date": int(data["date"][index]), "seconds": int(data["datesec"][index]),
                    "dynamic_fields": ranges,
                    "qbo_u0_max_abs_error_mps": error,
                    "qbo_u0_reconstruction_passed": error is not None and error <= TOLERANCE_MPS,
                    "outside_support_max_abs_tendency": outside_max,
                    "nonzero_tendency_present": bool(np.any(values["QBOTEND"][support, :] != 0)),
                    "tropopause": {},
                }
                for label, prefix in (("combined", "TROP"), ("primary", "TROPP"), ("cold_point", "TROPF")):
                    pressure = np.ma.asarray(data[prefix + "_P"][index])
                    found = np.ma.asarray(data[prefix + "_FD"][index])
                    valid = (~np.ma.getmaskarray(pressure) & ~np.ma.getmaskarray(found)
                             & (found.filled(0) > .5))
                    overlap = tropopause_overlap(
                        pressure.filled(np.nan) / 100, valid, weights,
                        float(mid[last]), float(edges[last + 1]),
                        np.isfinite(values["QBOTEND"][last]) & (values["QBOTEND"][last] != 0),
                    )
                    record["tropopause"][label] = overlap
                records.append(record)
    records.sort(key=lambda row: (row["date"], row["seconds"], row["nstep"]))
    integrated = [row for row in records if row["nstep"] > 0]
    passed = (
        [(row["nstep"], row["date"], row["seconds"]) for row in integrated]
        == [(1, 10901, 1800), (2, 10901, 3600)]
        and all(
            row["qbo_u0_reconstruction_passed"] and row["outside_support_max_abs_tendency"] == 0
            and row["nonzero_tendency_present"]
            and all(field["finite"] for field in row["dynamic_fields"].values())
            for row in integrated
        )
    )
    direct = any(
        (diagnosis["nonzero_tendency_below_fraction"] or 0) > 0
        for row in integrated for diagnosis in row["tropopause"].values()
    )
    result = {
        "scope": "two_step_technical_diagnostic_not_scientific_pilot", "phase": args.phase,
        "coupler_termination_verified": True, "integrated_steps_interface_passed": passed,
        "qbo_u0_absolute_tolerance_mps": TOLERANCE_MPS,
        "pressure_equality_tolerance_hpa": PRESSURE_TOLERANCE_HPA,
        "all_records_reconstruction_passed": all(row["qbo_u0_reconstruction_passed"] for row in records),
        "direct_tropospheric_nudging_observed": direct,
        "scientific_use_permitted": False, "seasonal_qbo_fidelity_verified": False,
        "files": files, "records": records,
        "target_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
        "analysis_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "helper_sha256": hashlib.sha256(
            Path(__file__).with_name("audit_remote_qbo.py").read_bytes()
        ).hexdigest(),
        "software": {"numpy": np.__version__, "netCDF4": netCDF4.__version__,
                     "cftime": cftime.__version__},
    }
    output = args.package / "runtime_assessment.json"
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: result[key] for key in (
        "phase", "integrated_steps_interface_passed", "all_records_reconstruction_passed",
        "direct_tropospheric_nudging_observed",
    )}, indent=2))
    for row in records:
        print(row["nstep"], row["seconds"], row["qbo_u0_max_abs_error_mps"], row["tropopause"]["combined"])
    if not passed:
        raise SystemExit("Integrated-step interface checks failed; see the saved assessment.")


if __name__ == "__main__":
    main()

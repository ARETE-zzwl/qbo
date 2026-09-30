"""Offline operator/control-cost assessment, not an integrated QBO experiment."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import netCDF4
import numpy as np

from qbo.metrics import retention_metrics

from qbo.forcing import forcing_layers, PRESSURE_TOLERANCE_HPA
from qbo.boundary import column_taper, safe_pressure, LOG_WIDTH, MARGIN_HPA


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    csv_path = args.output.with_suffix(".csv")
    if args.output.exists() or csv_path.exists():
        raise FileExistsError("Use a new JSON/CSV output pair.")
    root = Path(os.environ.get("QBO_WORKSPACE", Path(__file__).resolve().parents[2]))
    previous = json.loads((root / "outputs/20260914_v3_audit/assessment.json").read_text())
    records, profiles, inputs = [], [], {}
    for relative, expected_hash in previous["history_sha256"].items():
        path = root / relative
        actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual_hash != expected_hash:
            raise ValueError("Archived history changed: " + relative)
        inputs[relative] = actual_hash
        with netCDF4.Dataset(path) as data:
            p0 = float(data["P0"][()])
            hyam, hybm = np.asarray(data["hyam"][:]), np.asarray(data["hybm"][:])
            hyai, hybi = np.asarray(data["hyai"][:]), np.asarray(data["hybi"][:])
            reference_mid = (hyam + hybm) * p0 / 100
            layers = forcing_layers(reference_mid, [10, 100])
            base = np.zeros(len(reference_mid))
            for layer in layers:
                base[layer["index"] - 1] = layer["relative_strength"]
            lat = np.asarray(data["lat"][:])
            radians = np.deg2rad(lat)
            latitude_rate = np.exp(-radians ** 2 / (2 * .174532925 ** 2))
            latitude_rate[np.abs(radians) <= .035] = 1
            latitude_rate[np.abs(radians) > .384] = 0
            bands = {
                "forced_band": np.abs(radians) <= .384,
                "equatorial_5deg": np.abs(lat) <= 5,
                "central_2deg": np.abs(lat) <= 2,
            }
            gw = np.asarray(data["gw"][:])
            eligible = (base[:, None, None] > 0) & (latitude_rate[None, :, None] > 0)
            for index, step in enumerate(data["nsteph"][:]):
                ps = np.ma.asarray(data["PS"][index], dtype=float).filled(np.nan)
                if not np.isfinite(ps).all():
                    raise ValueError("Nonfinite surface pressure.")
                mid = (hyam[:, None, None] * p0 + hybm[:, None, None] * ps) / 100
                bottom = (hyai[1:, None, None] * p0 + hybi[1:, None, None] * ps) / 100
                trop = np.stack([
                    np.ma.asarray(data[name + "_P"][index], dtype=float).filled(np.nan) / 100
                    for name in ("TROP", "TROPP", "TROPF")
                ])
                flags = np.stack([
                    np.ma.asarray(data[name + "_FD"][index], dtype=float).filled(np.nan)
                    for name in ("TROP", "TROPP", "TROPF")
                ])
                safe = safe_pressure(trop, flags)
                taper = column_taper(bottom, safe)
                hard = np.isfinite(safe) & (bottom + MARGIN_HPA < safe)
                weight = np.where(eligible, taper, 0)
                if not np.isfinite(weight).all() or np.any((weight < 0) | (weight > 1)):
                    raise ValueError("Invalid operator weight.")
                forbidden = ~np.isfinite(safe) | (bottom + MARGIN_HPA >= safe)
                if np.any(weight[forbidden] != 0):
                    raise ValueError("Nonzero weight on a forbidden layer.")
                area = np.broadcast_to((gw * bands["forced_band"])[:, None], safe.shape)
                boundary = {}
                for d, name in enumerate(("TROP", "TROPP", "TROPF")):
                    valid = np.isfinite(trop[d]) & (trop[d] > 0) & (flags[d] > .5)
                    boundary[name] = {
                        "valid_area_fraction": float(np.sum(area * valid) / np.sum(area)),
                        "nonzero_midpoint_overlap_cells": int(np.count_nonzero(
                            (weight > 0) & valid & (mid > trop[d] + PRESSURE_TOLERANCE_HPA))),
                        "nonzero_interface_overlap_cells": int(np.count_nonzero(
                            (weight > 0) & valid & (bottom > trop[d] + PRESSURE_TOLERANCE_HPA))),
                    }
                    if boundary[name]["nonzero_midpoint_overlap_cells"] or boundary[name]["nonzero_interface_overlap_cells"]:
                        raise ValueError("Boundary overlap in candidate.")
                dt_rate = weight * base[:, None, None] * latitude_rate[None, :, None] * (1800 / 864000)
                if np.max(dt_rate) > 1800 / 864000:
                    raise ValueError("Scalar relaxation bound exceeded.")
                replay = None
                if int(step) > 0:
                    old = np.ma.asarray(data["QBOTEND"][index], dtype=float).filled(np.nan)
                    if not np.isfinite(old).all():
                        raise ValueError("Nonfinite old tendency cannot be replayed.")
                    new = old * weight
                    replay = {
                        "finite": bool(np.isfinite(new).all()),
                        "no_amplification": bool(np.all(np.abs(new) <= np.abs(old))),
                        "forbidden_max_abs_mps2": float(np.max(np.abs(new[forbidden]))),
                    }
                    if not replay["finite"] or not replay["no_amplification"] or replay["forbidden_max_abs_mps2"] != 0:
                        raise ValueError("Tendency replay check failed.")
                record_id = f"{path.parent.name}:{int(step)}"
                records.append({
                    "record": record_id, "file": relative, "nstep": int(step),
                    "all_definitions_valid_area_fraction": float(np.sum(area * np.isfinite(safe)) / np.sum(area)),
                    "boundary": boundary, "maximum_scalar_dt_rate": float(np.max(dt_rate)),
                    "offline_tendency_replay": replay,
                })
                for layer in layers:
                    k = layer["index"] - 1
                    for band, selection in bands.items():
                        profiles.append({
                            "record": record_id, "band": band, "level": k + 1,
                            "midpoint_hpa": float(reference_mid[k]),
                            "reference_bottom_hpa": float((hyai[k + 1] + hybi[k + 1]) * p0 / 100),
                            **retention_metrics(
                                taper[k, selection], hard[k, selection], latitude_rate[selection],
                                gw[selection], base[k]),
                        })
    if len(records) != 6 or sum(row["offline_tendency_replay"] is not None for row in records) != 4:
        raise ValueError("Expected six geometry records and four integrated tendency replays.")
    near70 = []
    for band in ("forced_band", "equatorial_5deg", "central_2deg"):
        for level in (46, 47):
            rows = [row for row in profiles if row["band"] == band and row["level"] == level]
            near70.append({
                "band": band, "level": level, "midpoint_hpa": rows[0]["midpoint_hpa"],
                **{name + "_range": [min(row[name] for row in rows), max(row[name] for row in rows)]
                   for name in ("retained_rate_fraction", "hard_mask_upper_bound",
                                "equivalent_relaxation_days", "longitude_modulation_rms_ratio")},
            })
    result = {
        "scope": "offline_operator_replay_on_v2_states_not_v4_model_integration",
        "formula": {"margin_hpa": MARGIN_HPA, "log_pressure_width": LOG_WIDTH},
        "geometry_and_scalar_checks_passed": True,
        "full_model_boundary_verified": False,
        "seasonal_qbo_fidelity_verified": False,
        "full_model_numerical_stability_verified": False,
        "scientific_pilot_approved": False,
        "near70_control": near70, "records": records,
        "profiles": profiles, "history_sha256": inputs,
        "code_sha256": {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                       for name in ("scripts/qbo_boundary_v4.py", "scripts/evaluate_boundary_v4.py",
                                    "BOUNDARY_TAPER_V4_PROTOCOL.md")},
    }
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    with csv_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(profiles[0]))
        writer.writeheader()
        writer.writerows(profiles)
    print(json.dumps({"records": len(records), "near70_control": near70}, indent=2))


if __name__ == "__main__":
    main()

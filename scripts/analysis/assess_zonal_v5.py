"""One reviewed uniform candidate and a necessary bound, not a climate run."""
import argparse
import csv
import json
import os
from pathlib import Path

import netCDF4
import numpy as np

from audit_qbo_operator_scope import digest, operator_metrics, RATE_TOLERANCE
from qbo.forcing import forcing_layers, PRESSURE_TOLERANCE_HPA
from qbo.metrics import retention_metrics
from qbo.boundary import column_taper, safe_pressure, MARGIN_HPA
from qbo.zonal import zonal_envelope

ROOT = Path(os.environ.get("QBO_WORKSPACE", Path(__file__).resolve().parents[2]))


def assess():
    previous_path = ROOT / "outputs/20260914_boundary_taper_v4/offline_assessment.json"
    previous = json.loads(previous_path.read_text(encoding="utf-8"))
    for name, expected in previous["code_sha256"].items():
        if digest(ROOT / name) != expected:
            raise ValueError("Frozen v4 material changed: " + name)
    old_profiles = {(p["record"], p["band"], p["level"]): p for p in previous["profiles"]}
    records, profiles, rings, native = [], [], [], []
    longitude_grid = None
    for relative, expected in previous["history_sha256"].items():
        path = ROOT / relative
        if digest(path) != expected:
            raise ValueError("Archived history changed: " + relative)
        with netCDF4.Dataset(path) as data:
            lon, lat = np.asarray(data["lon"][:]), np.asarray(data["lat"][:])
            if not np.allclose(np.diff(lon), 360 / len(lon), atol=1e-10, rtol=0):
                raise ValueError("Complete uniformly spaced longitude rings required.")
            if longitude_grid is not None and not np.array_equal(lon, longitude_grid):
                raise ValueError("History longitude grids differ.")
            longitude_grid = lon
            p0 = float(data["P0"][()])
            hyam, hybm = np.asarray(data["hyam"][:]), np.asarray(data["hybm"][:])
            hyai, hybi = np.asarray(data["hyai"][:]), np.asarray(data["hybi"][:])
            reference_mid = (hyam + hybm) * p0 / 100
            layers = forcing_layers(reference_mid, [10, 100])
            strength = np.zeros(len(reference_mid))
            for layer in layers:
                strength[layer["index"] - 1] = layer["relative_strength"]
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
            base_rate = strength[:, None, None] * latitude_rate[None, :, None] / 864000
            for index, step in enumerate(data["nsteph"][:]):
                record = f"{path.parent.name}:{int(step)}"
                ps = np.ma.asarray(data["PS"][index], dtype=float).filled(np.nan)
                if not np.isfinite(ps).all():
                    raise ValueError("Nonfinite surface pressure.")
                mid = (hyam[:, None, None] * p0 + hybm[:, None, None] * ps) / 100
                bottom = (hyai[1:, None, None] * p0 + hybi[1:, None, None] * ps) / 100
                trop = np.stack([
                    np.ma.asarray(data[name + "_P"][index], dtype=float).filled(np.nan) / 100
                    for name in ("TROP", "TROPP", "TROPF")
                ])
                found = np.stack([
                    np.ma.asarray(data[name + "_FD"][index], dtype=float).filled(np.nan)
                    for name in ("TROP", "TROPP", "TROPF")
                ])
                safe = safe_pressure(trop, found)
                column = column_taper(bottom, safe)
                uniform = np.broadcast_to(zonal_envelope(column), column.shape)
                hard = np.isfinite(safe) & (bottom > 0) & (bottom + MARGIN_HPA < safe)
                hard_uniform = np.broadcast_to(np.all(hard, axis=-1, keepdims=True), hard.shape)
                rate = base_rate * uniform
                if (np.any(uniform > column) or np.any(uniform[~hard] != 0)
                        or np.any(np.ptp(uniform, axis=-1) != 0)
                        or np.max(1800 * rate) > 1800 / 864000):
                    raise ValueError("Envelope, boundary or scalar check failed.")
                boundary = {}
                for d, name in enumerate(("TROP", "TROPP", "TROPF")):
                    valid = np.isfinite(trop[d]) & (trop[d] > 0) & (found[d] > .5)
                    boundary[name] = {
                        "nonzero_midpoint_overlap_cells": int(np.count_nonzero(
                            (rate > 0) & valid & (mid > trop[d] + PRESSURE_TOLERANCE_HPA))),
                        "nonzero_interface_overlap_cells": int(np.count_nonzero(
                            (rate > 0) & valid & (bottom > trop[d] + PRESSURE_TOLERANCE_HPA))),
                    }
                    if any(boundary[name].values()):
                        raise ValueError("Nonzero forbidden tendency support.")
                records.append({
                    "record": record, "file": relative, "nstep": int(step),
                    "boundary": boundary, "max_scalar_dt_rate": float(np.max(1800 * rate)),
                })
                for layer in layers:
                    k = layer["index"] - 1
                    selected = np.flatnonzero(bands["forced_band"])
                    native.extend(column[k, selected])
                    for j in selected:
                        clearance = safe[j] - bottom[k, j] - MARGIN_HPA
                        rings.append({
                            "record": record, "level": k + 1, "latitude_deg": float(lat[j]),
                            "midpoint_hpa": float(reference_mid[k]),
                            "v5_weight": float(uniform[k, j, 0]),
                            "hard_uniform_upper_bound": int(hard_uniform[k, j, 0]),
                            "valid_column_fraction": float(np.mean(np.isfinite(safe[j]))),
                            "safe_column_fraction": float(np.mean(hard[k, j])),
                            "minimum_clearance_hpa": float(np.min(clearance))
                            if np.isfinite(clearance).all() else None,
                        })
                    for band, selection in bands.items():
                        metrics = retention_metrics(
                            uniform[k, selection], hard_uniform[k, selection],
                            latitude_rate[selection], gw[selection], strength[k],
                        )
                        metrics["hard_uniform_rate_upper_bound"] = metrics.pop("hard_mask_upper_bound")
                        old = old_profiles[record, band, k + 1]
                        check_old = retention_metrics(
                            column[k, selection], hard[k, selection],
                            latitude_rate[selection], gw[selection], strength[k],
                        )
                        if abs(check_old["retained_rate_fraction"] - old["retained_rate_fraction"]) > 1e-14:
                            raise ValueError("Inherited column audit does not reproduce.")
                        decomposition = operator_metrics(rate[k, selection], gw[selection])
                        if decomposition["eddy_rms_per_s"] > RATE_TOLERANCE:
                            raise ValueError("Nonzonal rate in uniform envelope.")
                        profiles.append({
                            "record": record, "band": band, "level": k + 1,
                            "midpoint_hpa": float(reference_mid[k]), **metrics,
                            "v4_retained_rate_fraction": old["retained_rate_fraction"],
                            "v5_to_v4_rate_ratio": metrics["retained_rate_fraction"] / old["retained_rate_fraction"]
                            if old["retained_rate_fraction"] else None,
                            "v4_column_hard_upper_bound": old["hard_mask_upper_bound"],
                            "latitude_count": int(np.count_nonzero(selection)),
                            "positive_weight_latitudes": int(np.count_nonzero(uniform[k, selection, 0])),
                            "hard_safe_latitudes": int(np.count_nonzero(hard_uniform[k, selection, 0])),
                            "eddy_rate_rms_per_s": decomposition["eddy_rms_per_s"],
                        })
    if len(records) != 6 or len(profiles) != 252 or len(native) != len(rings):
        raise ValueError("Incomplete archived-state coverage.")
    summaries = []
    for band in bands:
        for level in sorted({row["level"] for row in profiles}):
            rows = [row for row in profiles if row["band"] == band and row["level"] == level]
            summary = {
                "band": band, "level": level, "midpoint_hpa": rows[0]["midpoint_hpa"],
                "records": len(rows),
                "zero_control_records": sum(row["retained_rate_fraction"] == 0 for row in rows),
            }
            for name in ("retained_rate_fraction", "hard_uniform_rate_upper_bound",
                         "equivalent_relaxation_days", "positive_weight_latitudes", "hard_safe_latitudes"):
                values = [row[name] for row in rows if row[name] is not None]
                summary[name + "_range"] = [min(values), max(values)] if values else None
            summaries.append(summary)
    source_names = [
        "scripts/qbo_zonal_v5.py", "scripts/assess_zonal_v5.py",
        "scripts/audit_qbo_operator_scope.py", "scripts/audit_remote_qbo.py",
        "ZONAL_ENVELOPE_V5_PROTOCOL.md",
    ]
    return {
        "scope": "offline_longitude_uniform_candidate_and_necessary_support_bound_not_climate_run",
        "candidate": "min_longitude_of_frozen_v4_weight",
        "boundary_uniformity_dominance_scalar_checks_passed": True,
        "full_model_boundary_verified": False,
        "seasonal_qbo_fidelity_verified": False,
        "coupled_stability_verified": False,
        "scientific_pilot_approved": False,
        "records": records, "profiles": profiles, "rings": rings, "summaries": summaries,
        "history_sha256": previous["history_sha256"],
        "source_sha256": {
            **previous["code_sha256"], **{name: digest(ROOT / name) for name in source_names},
            str(previous_path.relative_to(ROOT)): digest(previous_path),
        },
    }, np.asarray(native)


def write_output(output):
    if output.exists():
        raise FileExistsError("Use a fresh v5 output directory.")
    result, archived = assess()
    nlon = archived.shape[1]
    controls = np.array([np.zeros(nlon), np.ones(nlon), np.full(nlon, .5),
                         np.r_[0, np.ones(nlon - 1)], np.r_[np.ones(nlon - 1), 0],
                         np.linspace(0, 1, nlon)])
    synthetic = np.vstack([controls, np.random.default_rng(20260914).uniform(size=(128, nlon))])
    vectors = np.vstack([archived, synthetic])
    output.mkdir(parents=True)
    for name, values in (("native_weights.txt", vectors), ("native_expected.txt", vectors.min(axis=1))):
        with (output / name).open("x", encoding="ascii", newline="\n") as handle:
            if values.ndim == 2:
                handle.write(f"{len(values)} {nlon}\n")
            np.savetxt(handle, values, fmt="%.17e")
    result["native_probe"] = {
        "archived_rings": len(archived), "synthetic_rings": len(synthetic),
        "total_rings": len(vectors), "longitudes": nlon, "seed": 20260914,
        "files_sha256": {name: digest(output / name)
                         for name in ("native_weights.txt", "native_expected.txt")},
    }
    with (output / "assessment.json").open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    for name, rows in (("profiles.csv", result["profiles"]), ("rings.csv", result["rings"])):
        with (output / name).open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    result = write_output(parser.parse_args().output)
    print(json.dumps({
        "records": len(result["records"]), "native_probe": result["native_probe"],
        "equatorial_5deg": [row for row in result["summaries"] if row["band"] == "equatorial_5deg"],
    }, indent=2))

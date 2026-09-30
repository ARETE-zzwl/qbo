"""Zonal/eddy forcing decomposition on archived states, not a v4 simulation."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import netCDF4
import numpy as np

from qbo.metrics import operator_metrics

from qbo.forcing import forcing_layers
from qbo.boundary import column_taper, safe_pressure

ROOT = Path(os.environ.get("QBO_WORKSPACE", Path(__file__).resolve().parents[2]))
RATE_TOLERANCE = 1e-18
RELATIVE_TOLERANCE = 1e-12


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assess():
    previous_path = ROOT / "outputs/20260914_boundary_taper_v4/offline_assessment.json"
    previous = json.loads(previous_path.read_text(encoding="utf-8"))
    for name, expected in previous["code_sha256"].items():
        if digest(ROOT / name) != expected:
            raise ValueError("Frozen source/protocol changed: " + name)
    qbo = ROOT / "outputs/20260914_cam_adapter_v4/payload_r2/SourceMods/src.cam/qbo.F90"
    if digest(qbo) != "a1874a9c264cc7d2786286de9d686e7059700f14697d593cf5e84871000b66fd":
        raise ValueError("Frozen generated QBO source changed.")
    old_profiles = {(p["record"], p["band"], p["level"]): p for p in previous["profiles"]}
    profiles, records = [], []
    for relative, expected in previous["history_sha256"].items():
        path = ROOT / relative
        if digest(path) != expected:
            raise ValueError("Archived history changed: " + relative)
        with netCDF4.Dataset(path) as data:
            lon, lat = np.asarray(data["lon"][:]), np.asarray(data["lat"][:])
            if not np.allclose(np.diff(lon), 360 / len(lon), atol=1e-10, rtol=0):
                raise ValueError("Expected a complete uniformly spaced longitude grid.")
            p0 = float(data["P0"][()])
            midpoint = (np.asarray(data["hyam"][:]) + np.asarray(data["hybm"][:])) * p0 / 100
            hyai, hybi = np.asarray(data["hyai"][:]), np.asarray(data["hybi"][:])
            layers = forcing_layers(midpoint, [10, 100])
            radians = np.deg2rad(lat)
            lat_rate = np.exp(-radians ** 2 / (2 * .174532925 ** 2))
            lat_rate[np.abs(radians) <= .035] = 1
            lat_rate[np.abs(radians) > .384] = 0
            bands = {
                "forced_band": np.abs(radians) <= .384,
                "equatorial_5deg": np.abs(lat) <= 5,
                "central_2deg": np.abs(lat) <= 2,
            }
            area_weights = np.asarray(data["gw"][:])
            for index, step in enumerate(data["nsteph"][:]):
                record = f"{path.parent.name}:{int(step)}"
                records.append({"record": record, "file": relative, "nstep": int(step)})
                ps = np.ma.asarray(data["PS"][index], dtype=float).filled(np.nan)
                if not np.isfinite(ps).all():
                    raise ValueError("Nonfinite surface pressure.")
                bottom = (hyai[1:, None, None] * p0 + hybi[1:, None, None] * ps) / 100
                pressures = np.stack([
                    np.ma.asarray(data[name + "_P"][index], dtype=float).filled(np.nan) / 100
                    for name in ("TROP", "TROPP", "TROPF")
                ])
                found = np.stack([
                    np.ma.asarray(data[name + "_FD"][index], dtype=float).filled(np.nan)
                    for name in ("TROP", "TROPP", "TROPF")
                ])
                taper = column_taper(bottom, safe_pressure(pressures, found))
                for layer in layers:
                    k = layer["index"] - 1
                    original = np.broadcast_to(
                        layer["relative_strength"] * lat_rate[:, None] / 864000,
                        taper[k].shape,
                    )
                    rate = original * taper[k]
                    for band, selection in bands.items():
                        metrics = operator_metrics(rate[selection], area_weights[selection])
                        control = operator_metrics(original[selection], area_weights[selection])
                        if (metrics["maximum_reconstruction_error_per_s"] > RATE_TOLERANCE
                                or metrics["maximum_eddy_zonal_mean_per_s"] > RATE_TOLERANCE
                                or metrics["squared_norm_relative_error"] > RELATIVE_TOLERANCE
                                or control["eddy_rms_per_s"] > RATE_TOLERANCE):
                            raise ValueError("Projection or uniform-control check failed.")
                        old_ratio = old_profiles[record, band, k + 1]["longitude_modulation_rms_ratio"]
                        ratio = metrics["eddy_to_zonal_rms_ratio"]
                        if ((ratio is None) != (old_ratio is None)
                                or (ratio is not None and abs(ratio - old_ratio) > RELATIVE_TOLERANCE)):
                            raise ValueError("Independent decomposition disagrees with the archived audit.")
                        profiles.append({
                            "record": record, "nstep": int(step), "band": band,
                            "level": k + 1, "midpoint_hpa": float(midpoint[k]),
                            "relative_strength": layer["relative_strength"],
                            **metrics,
                            "original_eddy_rms_per_s": control["eddy_rms_per_s"],
                            "nonzonal_above_roundoff": metrics["eddy_rms_per_s"] > RATE_TOLERANCE,
                        })
    if len(records) != 6 or len(profiles) != 252:
        raise ValueError("Expected all six records, fourteen levels and three bands.")
    summaries = []
    for band in bands:
        for level in sorted({p["level"] for p in profiles}):
            rows = [p for p in profiles if p["band"] == band and p["level"] == level]
            summary = {
                "band": band, "level": level, "midpoint_hpa": rows[0]["midpoint_hpa"],
                "records": len(rows), "nonzonal_records": sum(p["nonzonal_above_roundoff"] for p in rows),
            }
            for name in ("eddy_rms_per_s", "eddy_to_zonal_rms_ratio", "eddy_squared_amplitude_fraction"):
                values = [p[name] for p in rows if p[name] is not None]
                summary[name + "_range"] = [min(values), max(values)] if values else None
            summaries.append(summary)
    review = any(p["nonzonal_above_roundoff"] for p in profiles)
    return {
        "scope": "unit_error_operator_decomposition_on_old_states_not_v4_model_output",
        "unit_target_minus_zonal_buffer_error_m_s": 1,
        "assumes_target_does_not_trigger_original_50_m_s_skip": True,
        "area_weights": "archived model gw; uniform longitude spacing",
        "numerical_tolerances": {"rate_per_s": RATE_TOLERANCE, "relative": RELATIVE_TOLERANCE},
        "decomposition_and_uniform_controls_passed": True,
        "agrees_with_previous_longitude_modulation_metric": True,
        "purely_zonal_imposed_momentum_tendency_preserved": not review,
        "requires_user_review_before_model_run_or_redesign": review,
        "actual_climate_effect_measured": False,
        "background_interaction_estimated": False,
        "records": records, "profiles": profiles, "summaries": summaries,
        "history_sha256": previous["history_sha256"],
        "source_sha256": {
            **previous["code_sha256"],
            str(qbo.relative_to(ROOT)): digest(qbo),
            str(previous_path.relative_to(ROOT)): digest(previous_path),
            "scripts/audit_qbo_operator_scope.py": digest(Path(__file__)),
            "DIRECTION_AND_OPERATOR_AUDIT_PROTOCOL.md": digest(ROOT / "DIRECTION_AND_OPERATOR_AUDIT_PROTOCOL.md"),
        },
    }


def write_assessment(output):
    csv_path = output.with_suffix(".csv")
    if output.exists() or csv_path.exists():
        raise FileExistsError("Use a new operator assessment JSON/CSV pair.")
    result = assess()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    with csv_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result["profiles"][0]))
        writer.writeheader()
        writer.writerows(result["profiles"])
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    result = write_assessment(parser.parse_args().output)
    print(json.dumps({
        "records": len(result["records"]), "profiles": len(result["profiles"]),
        "requires_user_review": result["requires_user_review_before_model_run_or_redesign"],
        "equatorial_5deg": [row for row in result["summaries"] if row["band"] == "equatorial_5deg"],
    }, indent=2))

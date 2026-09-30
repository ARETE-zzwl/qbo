"""Offline geometry screen only; do not generate or approve new QBO inputs."""
import argparse
import json
from pathlib import Path

import netCDF4
import numpy as np

from qbo.forcing import forcing_layers, tropopause_overlap, PRESSURE_TOLERANCE_HPA

CANDIDATE_BOTTOMS_HPA = (100, 70, 50)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("packages", type=Path, nargs="+")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    results = []
    for package in args.packages:
        for path in sorted(package.glob("*.nc")):
            with netCDF4.Dataset(path) as data:
                p0 = float(data["P0"][()]) / 100
                mid = (data["hyam"][:] + data["hybm"][:]) * p0
                edge = (data["hyai"][:] + data["hybi"][:]) * p0
                weights = data["gw"][:] * (np.abs(np.deg2rad(data["lat"][:])) <= .384)
                for index, step in enumerate(data["nsteph"][:]):
                    for cap in CANDIDATE_BOTTOMS_HPA:
                        layers = forcing_layers(mid, [10, cap])
                        full = [row for row in layers if row["relative_strength"] == 1]
                        last = layers[-1]["index"] - 1
                        for definition in ("TROP", "TROPP", "TROPF"):
                            pressure = np.ma.asarray(data[definition + "_P"][index])
                            found = np.ma.asarray(data[definition + "_FD"][index])
                            valid = (~np.ma.getmaskarray(pressure) & ~np.ma.getmaskarray(found)
                                     & (found.filled(0) > .5))
                            geometry = tropopause_overlap(
                                pressure.filled(np.nan) / 100, valid, weights,
                                float(mid[last]), float(edge[last + 1]), np.zeros(valid.shape, bool),
                            )
                            geometry.pop("nonzero_tendency_below_fraction")
                            results.append({
                                "package": package.name, "nstep": int(step),
                                "target_bottom_hpa": cap, "tropopause_definition": definition,
                                "full_strength_bottom_hpa": full[-1]["pressure_hpa"],
                                "buffer_midpoint_hpa": float(mid[last]),
                                "buffer_bottom_hpa": float(edge[last + 1]),
                                **geometry,
                            })
    summary = []
    for cap in CANDIDATE_BOTTOMS_HPA:
        for definition in ("TROP", "TROPP", "TROPF"):
            rows = [row for row in results if row["target_bottom_hpa"] == cap
                    and row["tropopause_definition"] == definition]
            if not rows:
                raise ValueError("No history records for boundary screening.")
            summary.append({
                "target_bottom_hpa": cap, "tropopause_definition": definition,
                "records": len(rows), "full_strength_bottom_hpa": rows[0]["full_strength_bottom_hpa"],
                "buffer_midpoint_hpa": rows[0]["buffer_midpoint_hpa"],
                "buffer_bottom_hpa": rows[0]["buffer_bottom_hpa"],
                "minimum_valid_area_fraction": min(row["valid_area_fraction"] for row in rows),
                "maximum_midpoint_below_fraction": max(
                    (row["midpoint_below_fraction"] for row in rows if row["midpoint_below_fraction"] is not None),
                    default=None),
                "maximum_interface_below_fraction": max(
                    (row["interface_below_fraction"] for row in rows if row["interface_below_fraction"] is not None),
                    default=None),
            })
    result = {
        "scope": "offline_geometry_only_no_new_inputs_no_new_integrations",
        "pressure_equality_tolerance_hpa": PRESSURE_TOLERANCE_HPA,
        "not_a_seasonal_safety_test": True, "no_candidate_approved": True,
        "summary": summary, "records": results,
    }
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

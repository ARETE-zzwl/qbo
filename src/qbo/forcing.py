"""WACCM layer support, reference diagnostics and tropopause overlap."""
import numpy as np

PRESSURE_TOLERANCE_HPA = 1e-3


def forcing_layers(mid_hpa, target_hpa):
    """Mirror pinned qbo.F90 layer selection and its two half-strength buffers."""
    pressure = np.asarray(mid_hpa, dtype=float)
    target = np.asarray(target_hpa, dtype=float)
    if (
        pressure.ndim != 1 or target.ndim != 1
        or pressure.size < 1 or target.size < 2
        or not np.isfinite(pressure).all() or not np.isfinite(target).all()
        or not (np.diff(pressure) > 0).all()
        or not (np.diff(target) > 0).all()
    ):
        raise ValueError("Pressure coordinates must be finite and increasing.")
    inside = np.flatnonzero((pressure >= target[0]) & (pressure <= target[-1]))
    if not inside.size:
        raise ValueError("No model layer lies inside the target pressure range.")
    first, last = int(inside[0]), int(inside[-1])
    return [
        {
            "index": k + 1,
            "pressure_hpa": float(pressure[k]),
            "relative_strength": 1.0 if first <= k <= last else 0.5,
        }
        for k in range(max(0, first - 1), min(pressure.size, last + 2))
    ]

def check_native_values(rows, pressure, dates, seconds, wind):
    wind = np.asarray(wind)
    nt, nz = wind.shape
    coordinates = np.column_stack([
        np.repeat(np.arange(1, nt + 1), nz),
        np.tile(np.arange(1, nz + 1), nt),
        np.repeat(dates, nz), np.repeat(seconds, nz), np.tile(pressure, nt),
    ])
    rows = np.asarray(rows)
    if rows.shape != (nt * nz, 6) or not np.array_equal(rows[:, :5], coordinates):
        raise ValueError("Native coordinate or array-shape mismatch.")
    if not np.isfinite(rows).all() or not np.array_equal(rows[:, 5], wind.ravel()):
        raise ValueError("Native wind values differ from the candidate file.")
    return {"values_checked": int(wind.size), "max_abs_error_mps": 0.0}

def expected_qbo_diagnostic(mid_hpa, lat_degrees, target_hpa, profile):
    """Pinned QBO_U0: copied edge wind, half-strength buffers, latitude taper."""
    mid = np.asarray(mid_hpa)
    layers = forcing_layers(mid, target_hpa)
    core = [row["index"] - 1 for row in layers if row["relative_strength"] == 1]
    wind = np.interp(np.clip(mid, mid[core[0]], mid[core[-1]]), target_hpa, profile)
    strength = np.zeros(mid.size)
    for row in layers:
        strength[row["index"] - 1] = row["relative_strength"]
    radians = np.deg2rad(lat_degrees)
    latitude = np.exp(-radians ** 2 / (2 * .174532925 ** 2))
    latitude[np.abs(radians) <= .035] = 1
    latitude[np.abs(radians) > .384] = 0
    strength[wind >= 50] = 0
    return (wind * strength)[:, None] * latitude[None, :]

def tropopause_overlap(trop_hpa, valid, lat_weights, midpoint, bottom, active):
    trop = np.asarray(trop_hpa)
    usable = np.asarray(valid) & np.isfinite(trop) & (trop > 0)
    weights = np.broadcast_to(np.asarray(lat_weights)[:, None], trop.shape)
    available = float(np.sum(weights[usable]))
    result = {
        "valid_area_fraction": available / float(np.sum(weights)),
        "midpoint_below_fraction": None,
        "interface_below_fraction": None,
        "nonzero_tendency_below_fraction": None,
    }
    if available:
        below = usable & (trop < midpoint - PRESSURE_TOLERANCE_HPA)
        result.update(
            midpoint_below_fraction=float(np.sum(weights[below]) / available),
            interface_below_fraction=float(np.sum(
                weights[usable & (trop < bottom - PRESSURE_TOLERANCE_HPA)]
            ) / available),
            nonzero_tendency_below_fraction=float(np.sum(weights[below & active]) / available),
        )
    return result

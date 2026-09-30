"""Frozen column-wise QBO mask; no target or relaxation-rate compensation."""
import numpy as np

MARGIN_HPA = 0.001
LOG_WIDTH = 0.2


def safe_pressure(pressure_hpa, found):
    pressure = np.ma.asarray(pressure_hpa, dtype=float).filled(np.nan)
    flags = np.ma.asarray(found, dtype=float).filled(np.nan)
    if pressure.ndim < 1 or pressure.shape[0] != 3 or flags.shape != pressure.shape:
        raise ValueError("Provide matching TROP, TROPP and TROPF arrays.")
    valid = np.all(
        np.isfinite(pressure) & (pressure > 0) & np.isfinite(flags) & (flags > .5),
        axis=0,
    )
    return np.where(valid, np.min(pressure, axis=0), np.nan)


def column_taper(bottom_hpa, safe_hpa):
    bottom, safe = np.broadcast_arrays(
        np.ma.asarray(bottom_hpa, dtype=float).filled(np.nan),
        np.ma.asarray(safe_hpa, dtype=float).filled(np.nan),
    )
    weight = np.zeros(bottom.shape)
    valid = (np.isfinite(bottom) & (bottom > 0) & np.isfinite(safe)
             & (safe > bottom + MARGIN_HPA))
    distance = np.log(safe[valid]) - np.log(bottom[valid] + MARGIN_HPA)
    x = np.clip(distance / LOG_WIDTH, 0, 1)
    weight[valid] = x * x * (3 - 2 * x)
    return weight

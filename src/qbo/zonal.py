"""Longitude-uniform envelope dominated by every frozen column weight."""
import numpy as np


def zonal_envelope(column_weights):
    weights = np.ma.asarray(column_weights, dtype=float).filled(np.nan)
    if (weights.ndim < 1 or weights.shape[-1] == 0
            or not np.isfinite(weights).all()
            or np.any((weights < 0) | (weights > 1))):
        raise ValueError("Provide finite weights in [0,1] with a longitude axis.")
    return np.min(weights, axis=-1, keepdims=True)

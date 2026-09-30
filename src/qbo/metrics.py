"""Area-weighted control retention and zonal/eddy decomposition."""
import numpy as np


def retention_metrics(taper, hard, latitude_rate, area_weights, original_strength):
    area = np.broadcast_to(np.asarray(area_weights)[:, None], taper.shape)
    rate = np.broadcast_to(np.asarray(latitude_rate)[:, None], taper.shape)
    denominator = float(np.sum(area * rate))
    mean_rate = float(np.sum(area * rate * taper) / np.sum(area))
    local_rate = rate * taper
    zonal_mean = local_rate.mean(axis=1, keepdims=True)
    modulation = float(np.sum(area * (local_rate - zonal_mean) ** 2))
    mean_square = float(np.sum(area * zonal_mean ** 2))
    return {
        "positive_area_fraction": float(np.sum(area * (taper > 0)) / np.sum(area)),
        "full_weight_area_fraction": float(np.sum(area * (taper == 1)) / np.sum(area)),
        "retained_rate_fraction": float(np.sum(area * local_rate) / denominator),
        "hard_mask_upper_bound": float(np.sum(area * rate * hard) / denominator),
        "equivalent_relaxation_days": 10 / (original_strength * mean_rate) if mean_rate else None,
        "longitude_modulation_rms_ratio": np.sqrt(modulation / mean_square) if mean_square else None,
    }


def operator_metrics(rate, latitude_weights):
    rate = np.asarray(rate, dtype=float)
    weights = np.asarray(latitude_weights, dtype=float)
    if (rate.ndim != 2 or weights.shape != (rate.shape[0],)
            or not np.isfinite(rate).all() or not np.isfinite(weights).all()
            or np.any(weights <= 0)):
        raise ValueError("Provide finite rates and positive matching latitude weights.")
    area = np.broadcast_to(weights[:, None], rate.shape)
    area_sum = area.sum()
    zonal = rate.mean(axis=1, keepdims=True)
    eddy = rate - zonal
    total_sq, zonal_sq, eddy_sq = [
        float(np.sum(area * field ** 2) / area_sum)
        for field in (rate, zonal, eddy)
    ]
    return {
        "area_mean_rate_per_s": float(np.sum(area * rate) / area_sum),
        "zonal_rms_per_s": float(np.sqrt(zonal_sq)),
        "eddy_rms_per_s": float(np.sqrt(eddy_sq)),
        "eddy_to_zonal_rms_ratio": float(np.sqrt(eddy_sq / zonal_sq)) if zonal_sq else None,
        "eddy_squared_amplitude_fraction": eddy_sq / total_sq if total_sq else None,
        "maximum_reconstruction_error_per_s": float(np.max(np.abs(rate - (zonal + eddy)))),
        "maximum_eddy_zonal_mean_per_s": float(np.max(np.abs(eddy.mean(axis=1)))),
        "squared_norm_relative_error": abs(total_sq - zonal_sq - eddy_sq) / total_sq if total_sq else 0.,
    }

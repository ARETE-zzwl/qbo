"""Diagnostic namelist checks and hard-cutoff geometry."""
import numpy as np

from .forcing import forcing_layers, tropopause_overlap

FIELDS = (
    "U", "T", "PS", "QBOTEND", "QBO_U0", "TROP_P", "TROP_FD",
    "TROPP_P", "TROPP_FD", "TROPF_P", "TROPF_FD",
)

def namelist_failures(values, target):
    def first(name):
        return next(iter(values.get(name, [])), "").strip("'\"")

    checks = {
        "forcing_disabled": first("qbo_use_forcing").lower() != ".true.",
        "noncyclic_not_explicit": first("qbo_cyclic").lower() != ".false.",
        "wrong_target": first("qbo_forcing_file") != target,
        "history_frequency": first("nhtfrq") != "1",
        "history_record_count": first("mfilt") != "2",
        "not_instantaneous": first("avgflag_pertape") != "I",
        "default_history_enabled": first("empty_htapes").lower() != ".true.",
        "missing_diagnostics": not set(FIELDS).issubset(
            {value.strip("'\"") for value in values.get("fincl1", [])}
        ),
    }
    return [name for name, failed in checks.items() if failed]

def hard_cutoff_strength(mid_hpa, target_hpa):
    strength = np.zeros(len(mid_hpa))
    for layer in forcing_layers(mid_hpa, target_hpa):
        if layer["relative_strength"] == 1:
            strength[layer["index"] - 1] = 1
    return strength

def cutoff_overlap(mid, edge, trop, valid, weights):
    last = np.flatnonzero(hard_cutoff_strength(mid, [10, 100]))[-1]
    result = tropopause_overlap(
        trop, valid, weights, float(mid[last]), float(edge[last + 1]),
        np.zeros(np.shape(valid), bool),
    )
    result.pop("nonzero_tendency_below_fraction")
    return {
        "lowest_active_level": int(last + 1),
        "lowest_active_midpoint_hpa": float(mid[last]),
        "lowest_active_bottom_hpa": float(edge[last + 1]),
        **result,
    }

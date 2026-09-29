"""Validate the Python HABR controller against the archived MATLAB output."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat

from habr_0d import (
    HABRParameters,
    baseline_currents,
    current_to_ml_min,
    simulate_habR_transition,
    update_after_portal_change,
)


ARCHIVE = Path(
    r"E:\GoogleDrive_backup\Student_Projects\Summer Projects\summer_2012\Liver_Circulation\HABR\archived_habr_run_output.mat"
)


def archived_values(path: Path) -> dict[str, float]:
    out = loadmat(path, squeeze_me=True, struct_as_record=False)["out"]
    return {field: float(getattr(out, field)) for field in out._fieldnames}


def compare(expected: float, observed: float, tolerance: float = 1e-9) -> dict[str, float | bool]:
    error = observed - expected
    scale = max(abs(expected), 1.0)
    return {"expected": float(expected), "observed": float(observed), "absolute_error": float(error), "pass": bool(abs(error) <= tolerance * scale)}


def main() -> int:
    expected = archived_values(ARCHIVE)
    params = HABRParameters()
    baseline = baseline_currents(params)
    update = update_after_portal_change(expected["newQpv_mL_min"], params, baseline)
    result = simulate_habR_transition(expected["newQpv_mL_min"], params=params)

    observed = {
        "oldQpv_mL_min": current_to_ml_min(baseline.i_pv),
        "oldQha_mL_min": current_to_ml_min(baseline.i_ha),
        "newQpv_mL_min": expected["newQpv_mL_min"],
        "newQha_mL_min": update["new_ha_ml_min"],
        "changePV_percent": update["percent_pv_change"],
        "changeHA_percent": update["percent_ha_change"],
        "Rpv_post": update["r_pv"],
        "Rha_post": update["r_ha"],
    }
    checks = {name: compare(expected[name], value) for name, value in observed.items()}
    checks["HV_conservation"] = {
        "maximum_absolute_error": float(np.max(np.abs(result["post"]["hv_ml_min"] - result["post"]["pv_ml_min"] - result["post"]["ha_ml_min"]))),
        "pass": bool(np.allclose(result["post"]["hv_ml_min"], result["post"]["pv_ml_min"] + result["post"]["ha_ml_min"])),
    }
    report = {
        "archive": str(ARCHIVE),
        "checks": checks,
        "trajectory_available": False,
        "note": "The archived MAT file contains scalar controller outputs only; full Simulink trajectories were not available for comparison.",
    }
    output = Path(__file__).with_name("habr_validation_report.json")
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if all(check["pass"] for check in checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())

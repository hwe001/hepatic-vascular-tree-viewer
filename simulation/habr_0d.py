"""Small, auditable Python translation of the archived hepatic HABR controller.

The source controller is ``HABR/runmodel.m`` plus the Simulink model
``HABRresistorsdebbaut.mdl``.  The two equations below are the dynamic core
implied by the algebraic flow calculation in ``runmodel.m``::

    L_pv dI_pv/dt = V_pv - R_pv I_pv - R_s (I_pv + I_ha)
    L_ha dI_ha/dt = V_ha + A sin(2 pi f t) - R_ha I_ha - R_s (I_pv + I_ha)

The conserved hepatic-vein current is ``I_hv = I_pv + I_ha``.  Currents use
the legacy Simulink convention (m^3/s); helper functions expose mL/min.

This is a faithful, solver-independent translation of the controller core,
not a claim that every internal Simscape RC state has been recreated.  The
archived model contains additional RC branches and initial-state machinery;
those are retained as parameter provenance below and should be added only
after trajectory comparison against a Simulink export.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Dict, Iterable, Optional

import numpy as np
from scipy.integrate import solve_ivp


ML_MIN_TO_M3_S = 1.0 / 60.0e6
M3_S_TO_ML_MIN = 60.0e6


@dataclass(frozen=True)
class HABRParameters:
    """Parameters copied from the archived MATLAB/Simulink model."""

    # Pressure-source values in the legacy model's pressure units.
    pv_source: float = 1.33322368e3
    ha_source: float = 13.3322368e3
    ha_pulse_amplitude: float = 2666.447368421
    ha_frequency_hz: float = 1.2

    # Dynamic PV/HA branches.
    r_pv: float = 37.807e6
    r_ha: float = 1.8931e9
    r_shared: float = 27.357e6
    l_pv: float = 50.0e6
    l_ha: float = 200.0e6

    # HV branch provenance from the Simulink model.  HV is conserved in the
    # reduced core; these values are available for a future explicit outlet.
    r_hv: float = 10.857e6
    l_hv: float = 30.0e6
    r_ivc: float = 16.5e6

    # RC branch provenance, copied without pretending they are used here.
    c_pv: float = 0.0094995e-6
    c_ha: float = 0.00012803e-6
    c_hv: float = 0.14224e-6
    r_cpv: float = 210.54e6
    r_cha: float = 1.562e10
    r_chv: float = 14.0607e6


@dataclass(frozen=True)
class HABRState:
    """A single reduced-model state in legacy current units."""

    i_pv: float
    i_ha: float


def ml_min_to_current(flow_ml_min: float) -> float:
    return np.asarray(flow_ml_min) * ML_MIN_TO_M3_S


def current_to_ml_min(current: float) -> float:
    converted = np.asarray(current) * M3_S_TO_ML_MIN
    return float(converted) if converted.ndim == 0 else converted


def baseline_currents(params: HABRParameters) -> HABRState:
    """Solve the same steady two-inlet closure used in ``runmodel.m``."""

    denominator = ((params.r_ha + params.r_pv) * params.r_shared
                   + params.r_ha * params.r_pv)
    i_ha = (params.ha_source * (params.r_shared + params.r_pv)
            - params.pv_source * params.r_shared) / denominator
    i_pv = (params.pv_source * (params.r_shared + params.r_ha)
            - params.ha_source * params.r_shared) / denominator
    return HABRState(i_pv=i_pv, i_ha=i_ha)


def habr_arterial_change(percent_pv_change: float) -> float:
    """Return the archived quadratic HA response in percent."""

    x = float(percent_pv_change)
    return 0.0007102 * x * x + 0.5492 * x


def update_after_portal_change(
    new_pv_ml_min: float,
    params: HABRParameters,
    old: Optional[HABRState] = None,
) -> Dict[str, float]:
    """Apply the exact HABR update logic from ``runmodel.m``.

    Returns updated flow and resistance values suitable for a post-resection
    solve.  The percentage convention follows the original script exactly:
    ``(old - new) / old * 100``.
    """

    old = old or baseline_currents(params)
    new_i_pv = ml_min_to_current(new_pv_ml_min)
    percent_pv_change = (old.i_pv - new_i_pv) / old.i_pv * 100.0
    percent_ha_change = habr_arterial_change(percent_pv_change)
    new_i_ha = old.i_ha * (1.0 + percent_ha_change / 100.0)

    r_pv = ((params.pv_source - new_i_pv * params.r_shared
             - new_i_ha * params.r_shared) / new_i_pv
            if new_i_pv else 9.0e19)
    r_ha = ((params.ha_source - new_i_pv * params.r_shared
             - new_i_ha * params.r_shared) / new_i_ha
            if new_i_ha else 9.0e19)
    return {
        "old_pv_ml_min": current_to_ml_min(old.i_pv),
        "old_ha_ml_min": current_to_ml_min(old.i_ha),
        "new_pv_ml_min": current_to_ml_min(new_i_pv),
        "new_ha_ml_min": current_to_ml_min(new_i_ha),
        "percent_pv_change": percent_pv_change,
        "percent_ha_change": percent_ha_change,
        "r_pv": r_pv,
        "r_ha": r_ha,
    }


def _rhs(t: float, y: np.ndarray, params: HABRParameters) -> np.ndarray:
    i_pv, i_ha = y
    shared = params.r_shared * (i_pv + i_ha)
    ha_pulse = params.ha_pulse_amplitude * np.sin(2.0 * np.pi * params.ha_frequency_hz * t)
    return np.array([
        (params.pv_source - params.r_pv * i_pv - shared) / params.l_pv,
        (params.ha_source + ha_pulse - params.r_ha * i_ha - shared) / params.l_ha,
    ])


def simulate(
    params: HABRParameters,
    duration_s: float = 60.0,
    sample_period_s: float = 0.02,
    initial: Optional[HABRState] = None,
) -> Dict[str, np.ndarray]:
    """Integrate the reduced PV/HA model and return all three flow traces."""

    initial = initial or baseline_currents(params)
    t_eval = np.arange(0.0, duration_s + sample_period_s * 0.5, sample_period_s)
    result = solve_ivp(
        lambda t, y: _rhs(t, y, params),
        (0.0, duration_s),
        [initial.i_pv, initial.i_ha],
        t_eval=t_eval,
        method="BDF",
        rtol=1.0e-8,
        atol=1.0e-14,
    )
    if not result.success:
        raise RuntimeError(result.message)
    i_pv, i_ha = result.y
    return {
        "time_s": result.t,
        "pv_ml_min": current_to_ml_min(i_pv),
        "ha_ml_min": current_to_ml_min(i_ha),
        "hv_ml_min": current_to_ml_min(i_pv + i_ha),
        "pv_pressure_legacy": params.pv_source - params.r_pv * i_pv,
        "ha_pressure_legacy": params.ha_source - params.r_ha * i_ha,
        "shared_pressure_legacy": params.r_shared * (i_pv + i_ha),
    }


def simulate_with_compliance(
    params: HABRParameters,
    duration_s: float = 60.0,
    sample_period_s: float = 0.02,
    initial: Optional[HABRState] = None,
) -> Dict[str, np.ndarray]:
    """Integrate an explicit RC-extended reconstruction of the circuit.

    The archived model contains paired RC branches for the PV, HA, and HV
    beds. Their exact Simscape connection graph is not represented in the
    scalar archive, so this function uses the standard series-inertance plus
    compliant-compartment interpretation:

    * PV and HA inertive branches feed compliant compartments.
    * Each compartment drains toward a shared sinusoidal-bed pressure through
      its archived RC resistance.
    * The HV branch drains the shared bed toward the IVC reference.

    This is an auditable reconstruction for transient exploration. It is not
    used by the exact scalar validator until its trajectories are compared
    with an exported Simulink trace.
    """

    initial = initial or baseline_currents(params)
    q_hv0 = initial.i_pv + initial.i_ha
    p_shared0 = (params.r_hv + params.r_ivc) * q_hv0
    p_pv0 = p_shared0 + params.r_cpv * initial.i_pv
    p_ha0 = p_shared0 + params.r_cha * initial.i_ha
    y0 = np.array([initial.i_pv, initial.i_ha, q_hv0, p_pv0, p_ha0, p_shared0])
    t_eval = np.arange(0.0, duration_s + sample_period_s * 0.5, sample_period_s)

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        i_pv, i_ha, i_hv, p_pv, p_ha, p_shared = y
        ha_pulse = params.ha_pulse_amplitude * np.sin(2.0 * np.pi * params.ha_frequency_hz * t)
        return np.array([
            (params.pv_source - params.r_pv * i_pv - p_pv) / params.l_pv,
            (params.ha_source + ha_pulse - params.r_ha * i_ha - p_ha) / params.l_ha,
            (p_shared - (params.r_hv + params.r_ivc) * i_hv) / params.l_hv,
            (i_pv - (p_pv - p_shared) / params.r_cpv) / params.c_pv,
            (i_ha - (p_ha - p_shared) / params.r_cha) / params.c_ha,
            ((p_pv - p_shared) / params.r_cpv
             + (p_ha - p_shared) / params.r_cha - i_hv) / params.c_hv,
        ])

    result = solve_ivp(
        rhs,
        (0.0, duration_s),
        y0,
        t_eval=t_eval,
        method="BDF",
        rtol=1.0e-7,
        atol=1.0e-14,
    )
    if not result.success:
        raise RuntimeError(result.message)
    i_pv, i_ha, i_hv, p_pv, p_ha, p_shared = result.y
    return {
        "time_s": result.t,
        "pv_ml_min": current_to_ml_min(i_pv),
        "ha_ml_min": current_to_ml_min(i_ha),
        "hv_ml_min": current_to_ml_min(i_hv),
        "pv_pressure_legacy": p_pv,
        "ha_pressure_legacy": p_ha,
        "shared_pressure_legacy": p_shared,
        "hv_balance_residual_ml_min": current_to_ml_min(i_pv + i_ha - i_hv),
    }


def simulate_habR_transition(
    new_pv_ml_min: float,
    duration_s: float = 150.0,
    sample_period_s: float = 0.02,
    params: Optional[HABRParameters] = None,
) -> Dict[str, np.ndarray]:
    """Run baseline and post-change traces using the archived HABR update."""

    params = params or HABRParameters()
    old = baseline_currents(params)
    update = update_after_portal_change(new_pv_ml_min, params, old)
    post_params = replace(params, r_pv=update["r_pv"], r_ha=update["r_ha"])
    baseline = simulate(params, duration_s=60.0, sample_period_s=sample_period_s, initial=old)
    post_initial = HABRState(
        i_pv=ml_min_to_current(update["new_pv_ml_min"]),
        i_ha=ml_min_to_current(update["new_ha_ml_min"]),
    )
    post = simulate(post_params, duration_s=duration_s, sample_period_s=sample_period_s, initial=post_initial)
    return {"baseline": baseline, "post": post, "habr_update": update}


def _print_summary(result: Dict[str, object]) -> None:
    update = result["habr_update"]
    baseline = result["baseline"]
    post = result["post"]
    print("Archived HABR controller translated to Python")
    # Report a short-window mean because the HA source includes a 1.2 Hz pulse.
    n = min(250, len(baseline["time_s"]))
    print(f"Baseline mean (last {n} samples): PV={np.mean(baseline['pv_ml_min'][-n:]):.3f}, HA={np.mean(baseline['ha_ml_min'][-n:]):.3f}, HV={np.mean(baseline['hv_ml_min'][-n:]):.3f} mL/min")
    print(f"HABR update: PV={update['new_pv_ml_min']:.3f}, HA={update['new_ha_ml_min']:.3f} mL/min")
    print(f"Updated resistances: Rpv={update['r_pv']:.6g}, Rha={update['r_ha']:.6g}")
    n = min(250, len(post["time_s"]))
    print(f"Post mean (last {n} samples): PV={np.mean(post['pv_ml_min'][-n:]):.3f}, HA={np.mean(post['ha_ml_min'][-n:]):.3f}, HV={np.mean(post['hv_ml_min'][-n:]):.3f} mL/min")


if __name__ == "__main__":
    _print_summary(simulate_habR_transition(new_pv_ml_min=1060.0))

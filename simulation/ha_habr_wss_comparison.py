"""Compare HA WSS before and after right hepatectomy with/without HABR.

This is a transparent reduced-order post-processing model. Mean flow is
allocated to each edge in proportion to the number of terminal sites below
that edge; WSS is then estimated with fully developed laminar flow,

    tau = 4 mu Q / (pi r^3).

The geometry is the Murray-regularised HA CCO asset. The calculation is not a
replacement for 3-D CFD: it reports a geometry-aware screening quantity for
the digital-twin scenarios and records all assumptions in the JSON output.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, Iterable

import numpy as np

from ha_remnant_postresection import retained_tree


ML_MIN_TO_M3_S = 1.0 / 60.0e6
PA_TO_DYN_CM2 = 10.0


def children_for(data: dict) -> Dict[int, list[dict]]:
    children: Dict[int, list[dict]] = {}
    for edge in data["edges"]:
        children.setdefault(int(edge["a"]), []).append(edge)
    return children


def terminal_counts(data: dict) -> Dict[int, int]:
    children = children_for(data)
    counts: Dict[int, int] = {}

    def visit(node_id: int) -> int:
        if node_id in counts:
            return counts[node_id]
        outgoing = children.get(node_id, [])
        counts[node_id] = 1 if not outgoing else sum(visit(int(edge["b"])) for edge in outgoing)
        return counts[node_id]

    root = int(data.get("rootId", 0))
    visit(root)
    return counts


def edge_rows(data: dict, root_flow_ml_min: float, label: str, mu: float, pulsatility: float, radius_floor_mm: float) -> list[dict]:
    children = children_for(data)
    counts = terminal_counts(data)
    root = int(data.get("rootId", 0))
    total_terminals = counts[root]
    root_radius_mm = max(float(children[root][0].get("radius", 0.0)) if children.get(root) else radius_floor_mm, radius_floor_mm)
    rows = []

    def walk(node_id: int, q_parent_ml_min: float) -> None:
        outgoing = children.get(node_id, [])
        if not outgoing:
            return
        weights = np.asarray([max(float(edge.get("radius", 0.0)), radius_floor_mm) ** 3 for edge in outgoing], dtype=float)
        weights /= weights.sum()
        for edge, branch_fraction in zip(outgoing, weights):
            child = int(edge["b"])
            # At a bifurcation, enforce the Murray-style Q proportional to r^3
            # split. Unary segments carry their parent flow unchanged.
            q_mean = q_parent_ml_min * float(branch_fraction) if len(outgoing) > 1 else q_parent_ml_min
            q_amp = q_mean * pulsatility
            raw_radius_mm = float(edge.get("radius", 0.0))
            # The CCO asset stops above the unresolved arteriolar bed. A
            # floor is used only for WSS post-processing, not to rewrite the
            # published geometry or its Murray-radius audit.
            flow_consistent_radius_mm = root_radius_mm * max(q_mean / root_flow_ml_min, 0.0) ** (1.0 / 3.0)
            effective_radius_mm = max(raw_radius_mm, radius_floor_mm, flow_consistent_radius_mm)
            radius_m = max(effective_radius_mm * 1.0e-3, 1.0e-9)
            coefficient = 4.0 * mu / (np.pi * radius_m**3)
            tau_mean_pa = coefficient * q_mean * ML_MIN_TO_M3_S
            tau_peak_pa = coefficient * (q_mean + q_amp) * ML_MIN_TO_M3_S
            rows.append({
                "scenario": label,
                "edge_index": len(rows),
                "a": int(edge["a"]),
                "b": child,
                "length_mm": float(edge.get("length", 0.0)),
                "radius_mm": raw_radius_mm,
                "effective_radius_mm_for_wss": effective_radius_mm,
                "terminal_sites_below": counts[child],
                "root_flow_ml_min": root_flow_ml_min,
                "mean_flow_ml_min": q_mean,
                "pulsatile_amplitude_ml_min": q_amp,
                "mean_wss_pa": tau_mean_pa,
                "peak_wss_pa": tau_peak_pa,
                "mean_wss_dyn_cm2": tau_mean_pa * PA_TO_DYN_CM2,
                "peak_wss_dyn_cm2": tau_peak_pa * PA_TO_DYN_CM2,
            })
            walk(child, q_mean)

    walk(root, root_flow_ml_min)
    return rows


def summary(rows: Iterable[dict]) -> dict:
    rows = list(rows)
    def stats(key: str) -> dict:
        values = np.asarray([row[key] for row in rows], dtype=float)
        return {
            "min": float(values.min()),
            "p05": float(np.percentile(values, 5)),
            "median": float(np.median(values)),
            "p95": float(np.percentile(values, 95)),
            "max": float(values.max()),
            "mean": float(values.mean()),
        }
    return {
        "edge_count": len(rows),
        "mean_flow_ml_min": stats("mean_flow_ml_min"),
        "mean_wss_pa": stats("mean_wss_pa"),
        "peak_wss_pa": stats("peak_wss_pa"),
        "mean_wss_dyn_cm2": stats("mean_wss_dyn_cm2"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", type=Path, default=Path("cco/data/cco_arterial15_repaired_murray.json"))
    parser.add_argument("--zones", type=Path, default=Path("data/portal_zones_volume.json"))
    parser.add_argument("--output", type=Path, default=Path("simulation/ha_habr_wss_comparison.json"))
    parser.add_argument("--csv", type=Path, default=Path("simulation/ha_habr_wss_edges.csv"))
    parser.add_argument("--mu", type=float, default=0.0035, help="dynamic blood viscosity in Pa s")
    parser.add_argument("--intact-flow", type=float, default=250.0)
    parser.add_argument("--post-habr-flow", type=float, default=150.0)
    parser.add_argument("--portal-flow", type=float, default=800.0)
    parser.add_argument("--pulsatility", type=float, default=0.25)
    parser.add_argument("--wss-radius-floor-mm", type=float, default=0.15, help="effective distal radius floor for unresolved arterioles")
    args = parser.parse_args()

    data = json.loads(args.geometry.read_text(encoding="utf-8"))
    zones = json.loads(args.zones.read_text(encoding="utf-8"))
    remnant, mask = retained_tree(data, zones, {2, 3, 4})
    remnant_fraction = mask["retained_terminal_count"] / mask["total_terminal_count"]
    remnant_pre_flow = args.intact_flow * remnant_fraction

    scenarios = {
        "intact_pre_hepatectomy": (data, args.intact_flow),
        "remnant_pre_hepatectomy_no_HABR": (remnant, remnant_pre_flow),
        "remnant_post_hepatectomy_no_HABR": (remnant, remnant_pre_flow),
        "remnant_post_hepatectomy_with_HABR": (remnant, args.post_habr_flow),
    }
    all_rows = []
    summaries = {}
    for label, (geometry, flow) in scenarios.items():
        rows = edge_rows(geometry, flow, label, args.mu, args.pulsatility, args.wss_radius_floor_mm)
        all_rows.extend(rows)
        summaries[label] = summary(rows)

    result = {
        "model": "reduced_order_HA_WSS_with_terminal_demand_allocation",
        "geometry": str(args.geometry),
        "radius_model": "Murray-regularised production asset",
        "assumptions": {
            "dynamic_viscosity_pa_s": args.mu,
            "wss_equation": "tau = 4*mu*Q/(pi*r^3)",
            "flow_allocation": "conservative branch split proportional to effective radius cubed; unary segments carry parent flow",
            "wss_radius_regularisation": "effective radius is also floored by the Murray-consistent radius implied by local flow fraction, preventing unresolved unary contractions from creating artificial WSS spikes",
            "ha_intact_mean_flow_ml_min": args.intact_flow,
            "portal_total_flow_ml_min": args.portal_flow,
            "post_habr_ha_mean_flow_ml_min": args.post_habr_flow,
            "pulsatility_fraction": args.pulsatility,
            "wss_radius_floor_mm": args.wss_radius_floor_mm,
            "wss_radius_floor_scope": "post-processing only; raw geometry remains unchanged",
            "habr_interpretation": "post-resection HA flow is an illustrative buffered response, not a calibrated patient measurement",
        },
        "remnant_mask": mask,
        "remnant_fraction_of_terminal_sites": remnant_fraction,
        "remnant_pre_flow_ml_min": remnant_pre_flow,
        "scenarios": summaries,
        "literature_anchor_note": "Published liver-flow studies commonly place HA near 20-25% of total afferent liver inflow and describe reciprocal HA recruitment after portal-flow reduction; the WSS values here remain model-derived because no subject-specific hepatic-artery WSS measurement is imposed.",
        "edges": all_rows,
    }
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_rows[0].keys())
        writer.writeheader()
        writer.writerows(all_rows)

    print(json.dumps({
        "remnant_terminal_sites": mask["retained_terminal_count"],
        "total_terminal_sites": mask["total_terminal_count"],
        "remnant_pre_flow_ml_min": remnant_pre_flow,
        "scenarios": summaries,
    }, indent=2))


if __name__ == "__main__":
    main()

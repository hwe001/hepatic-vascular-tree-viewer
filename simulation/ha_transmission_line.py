"""Linearized 1-D transmission-line solve for the hepatic-arterial CCO tree.

The geometry is read from the HA CCO JSON asset. A harmonic root-flow
perturbation is propagated through every branch using the linearized
continuity and momentum equations::

    dP/dx = -(R' + i*w*L') Q
    dQ/dx = -(G' + i*w*C') P

The prototype uses G'=0, a thin-wall compliance estimate, and terminal
Windkessel resistance/capacitance. These are explicit assumptions, not
subject-specific measurements. The root mean flow is kept separate from the
pulsatile perturbation so the waveform has a prescribed mean.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


ML_MIN_TO_M3_S = 1.0 / 60.0e6
M3_S_TO_ML_MIN = 60.0e6


@dataclass(frozen=True)
class Parameters:
    frequency_hz: float = 1.2
    root_mean_ml_min: float = 250.0
    root_pulsatility: float = 0.25
    blood_density_kg_m3: float = 1060.0
    dynamic_viscosity_pa_s: float = 0.0035
    wall_modulus_pa: float = 400_000.0
    wall_thickness_ratio: float = 0.10
    terminal_resistance_pa_s_m3: float = 3.0e8
    terminal_compliance_m3_pa: float = 2.0e-10


def _edge_geometry(edge: dict) -> Tuple[float, float]:
    return (
        max(float(edge.get("length", 0.0)) * 1.0e-3, 1.0e-6),
        max(float(edge.get("radius", 0.0)) * 1.0e-3, 1.0e-6),
    )


def _matrix(edge: dict, params: Parameters) -> np.ndarray:
    length_m, radius_m = _edge_geometry(edge)
    area_m2 = np.pi * radius_m**2
    omega = 2.0 * np.pi * params.frequency_hz
    resistance_per_m = 8.0 * params.dynamic_viscosity_pa_s / (np.pi * radius_m**4)
    inertance_per_m = params.blood_density_kg_m3 / area_m2
    compliance_per_m = 2.0 * np.pi * radius_m**3 / (
        params.wall_modulus_pa * params.wall_thickness_ratio * radius_m
    )
    z = resistance_per_m + 1j * omega * inertance_per_m
    y = 1j * omega * compliance_per_m
    gamma = np.sqrt(z * y)
    zc = np.sqrt(z / y)
    gl = gamma * length_m
    c, s = np.cosh(gl), np.sinh(gl)
    return np.array([[c, -zc * s], [-(s / zc), c]], dtype=complex)


def _terminal_impedance(params: Parameters) -> complex:
    omega = 2.0 * np.pi * params.frequency_hz
    return params.terminal_resistance_pa_s_m3 + 1.0 / (
        1j * omega * params.terminal_compliance_m3_pa
    )


def _root_candidates(data: dict) -> list[int]:
    incoming = {int(node["id"]): 0 for node in data["nodes"]}
    for edge in data["edges"]:
        incoming[int(edge["b"])] = incoming.get(int(edge["b"]), 0) + 1
    return [node_id for node_id, count in incoming.items() if count == 0]


def _root_id(data: dict, children: Dict[int, list]) -> tuple[int, list[int]]:
    if data.get("rootId") is not None:
        return int(data["rootId"]), [int(data["rootId"])]
    roots = _root_candidates(data)
    if not roots:
        raise ValueError("The HA asset has no root candidate")

    def component_size(root: int) -> int:
        seen = set()
        stack = [root]
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            stack.extend(child for _, child in children.get(node, []))
        return len(seen)

    # The current provisional asset contains several disconnected components.
    # Select the largest rooted component, and report the others in the result.
    return max(roots, key=component_size), roots


def solve(data: dict, params: Parameters) -> dict:
    edges = data["edges"]
    children: Dict[int, list] = {}
    for index, edge in enumerate(edges):
        children.setdefault(int(edge["a"]), []).append((index, int(edge["b"])))

    matrices = [_matrix(edge, params) for edge in edges]
    input_impedances = np.zeros(len(edges), dtype=complex)
    node_impedances: Dict[int, complex] = {}
    terminal_z = _terminal_impedance(params)

    def node_impedance(node_id: int) -> complex:
        if node_id in node_impedances:
            return node_impedances[node_id]
        outgoing = children.get(node_id, [])
        if not outgoing:
            node_impedances[node_id] = terminal_z
            return terminal_z
        child_z = []
        for edge_index, child in outgoing:
            load_z = node_impedance(child)
            a, b, c, d = matrices[edge_index].ravel()
            z_in = (a * load_z + b) / (c * load_z + d)
            input_impedances[edge_index] = z_in
            child_z.append(z_in)
        node_impedances[node_id] = 1.0 / sum(1.0 / z for z in child_z)
        return node_impedances[node_id]

    root, roots = _root_id(data, children)
    root_z = node_impedance(root)
    root_q = params.root_mean_ml_min * params.root_pulsatility * ML_MIN_TO_M3_S
    root_p = root_z * root_q
    segments = []

    def walk(node_id: int, p_node: complex) -> None:
        for edge_index, child in children.get(node_id, []):
            edge = edges[edge_index]
            q_up = p_node / input_impedances[edge_index]
            state_down = np.linalg.solve(
                matrices[edge_index], np.array([p_node, q_up], dtype=complex)
            )
            segments.append({
                "edge_index": edge_index,
                "a": int(edge["a"]),
                "b": int(edge["b"]),
                "length_mm": float(edge.get("length", 0.0)),
                "radius_mm": float(edge.get("radius", 0.0)),
                "flow_amplitude_ml_min": float(abs(q_up) * M3_S_TO_ML_MIN),
                "flow_phase_deg": float(np.angle(q_up, deg=True)),
                "pressure_amplitude_mmHg": float(abs(p_node) / 133.322),
                "pressure_phase_deg": float(np.angle(p_node, deg=True)),
                "flow_downstream_amplitude_ml_min": float(abs(state_down[1]) * M3_S_TO_ML_MIN),
            })
            walk(child, state_down[0])

    walk(root, root_p)
    flow_amplitudes = np.array([item["flow_amplitude_ml_min"] for item in segments])
    pressure_amplitudes = np.array([item["pressure_amplitude_mmHg"] for item in segments])
    component_nodes = set()
    stack = [root]
    while stack:
        node = stack.pop()
        if node in component_nodes:
            continue
        component_nodes.add(node)
        stack.extend(child for _, child in children.get(node, []))
    return {
        "model": "linearized_1d_frequency_domain_transmission_line",
        "root_id": root,
        "root_candidates": roots,
        "topology_note": "The provisional HA asset has disconnected rooted components; this result solves the largest component and does not silently connect the others.",
        "parameters": {key: getattr(params, key) for key in params.__dataclass_fields__},
        "root_mean_flow_ml_min": params.root_mean_ml_min,
        "root_pulsatile_amplitude_ml_min": params.root_mean_ml_min * params.root_pulsatility,
        "root_pressure_amplitude_mmHg": float(abs(root_p) / 133.322),
        "segment_count": len(segments),
        "component_node_count": len(component_nodes),
        "terminal_count": sum(1 for node in component_nodes if not children.get(node)),
        "flow_amplitude_range_ml_min": [float(flow_amplitudes.min()), float(flow_amplitudes.max())],
        "pressure_amplitude_range_mmHg": [float(pressure_amplitudes.min()), float(pressure_amplitudes.max())],
        "segments": segments,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", type=Path, default=Path("cco/data/cco_arterial15_provisional.json"))
    parser.add_argument("--output", type=Path, default=Path("simulation/ha_transmission_line_result.json"))
    args = parser.parse_args()
    data = json.loads(args.geometry.read_text(encoding="utf-8"))
    result = solve(data, Parameters())
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({key: result[key] for key in result if key != "segments"}, indent=2))


if __name__ == "__main__":
    main()

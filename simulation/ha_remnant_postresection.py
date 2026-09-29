"""Build and solve the Z2+Z3+Z4 post-resection HA sub-tree.

The zone mask is deliberately labelled exploratory: terminal sites are assigned
to the nearest portal-zone centroid, and every root-to-retained-terminal path
is retained. The resulting connected tree is then passed to the same
frequency-domain transmission-line solver used for the intact HA tree.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np

from ha_transmission_line import Parameters, solve


def centres(zone_data: dict) -> list[np.ndarray]:
    result = []
    for zone in zone_data["zones"]:
        values = np.asarray(zone["vertices"], dtype=float).reshape(-1, 3)
        result.append(values.mean(axis=0))
    return result


def retained_tree(data: dict, zone_data: dict, retained_zones: set[int]) -> tuple[dict, dict]:
    nodes = {int(node["id"]): node for node in data["nodes"]}
    children = {node_id: [] for node_id in nodes}
    parents = {node_id: [] for node_id in nodes}
    for edge in data["edges"]:
        children[int(edge["a"])].append(edge)
        parents[int(edge["b"])].append(int(edge["a"]))

    zone_centres = centres(zone_data)
    zone_of = {}
    for node_id, node in nodes.items():
        point = np.asarray(node["p"], dtype=float)
        zone_of[node_id] = 1 + min(
            range(len(zone_centres)),
            key=lambda index: float(np.sum((point - zone_centres[index]) ** 2)),
        )

    terminals = [node_id for node_id, outgoing in children.items() if not outgoing]
    retained_terminals = [node_id for node_id in terminals if zone_of[node_id] in retained_zones]
    needed = set(retained_terminals)
    stack = list(retained_terminals)
    while stack:
        node_id = stack.pop()
        for parent in parents[node_id]:
            if parent not in needed:
                needed.add(parent)
                stack.append(parent)

    filtered = copy.deepcopy(data)
    filtered["nodes"] = [node for node in data["nodes"] if int(node["id"]) in needed]
    filtered["edges"] = [
        edge for edge in data["edges"]
        if int(edge["a"]) in needed and int(edge["b"]) in needed
    ]
    filtered["rootId"] = int(data.get("rootId", 0))
    meta = {
        "retained_zones": sorted(retained_zones),
        "zone_assignment": "nearest portal-zone centroid",
        "total_terminal_count": len(terminals),
        "retained_terminal_count": len(retained_terminals),
        "retained_node_count": len(filtered["nodes"]),
        "retained_edge_count": len(filtered["edges"]),
        "retained_terminal_ids": retained_terminals,
    }
    return filtered, meta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geometry", type=Path, default=Path("cco/data/cco_arterial15_repaired.json"))
    parser.add_argument("--zones", type=Path, default=Path("data/portal_zones_volume.json"))
    parser.add_argument("--output", type=Path, default=Path("simulation/ha_remnant_postresection_result.json"))
    parser.add_argument("--post-ha-mean", type=float, default=150.0)
    args = parser.parse_args()

    data = json.loads(args.geometry.read_text(encoding="utf-8"))
    zone_data = json.loads(args.zones.read_text(encoding="utf-8"))
    filtered, mask = retained_tree(data, zone_data, {2, 3, 4})

    params = Parameters(root_mean_ml_min=args.post_ha_mean)
    result = solve(filtered, params)
    result["post_resection_mask"] = mask
    result["scenario"] = {
        "portal_total_flow_ml_min": 800.0,
        "intact_ha_mean_ml_min": 250.0,
        "post_resection_habr_ha_mean_ml_min": args.post_ha_mean,
        "interpretation": "Illustrative HABR scenario; values are not subject-specific calibration.",
    }
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({
        "root_id": result["root_id"],
        "retained_terminals": mask["retained_terminal_count"],
        "retained_edges": mask["retained_edge_count"],
        "post_ha_mean_ml_min": args.post_ha_mean,
        "flow_amplitude_range_ml_min": result["flow_amplitude_range_ml_min"],
        "pressure_amplitude_range_mmHg": result["pressure_amplitude_range_mmHg"],
    }, indent=2))


if __name__ == "__main__":
    main()

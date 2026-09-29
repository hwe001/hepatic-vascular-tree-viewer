"""Apply a conservative Murray-style radius regularisation to the repaired HA tree."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path


def smooth(data: dict) -> tuple[dict, dict]:
    out = copy.deepcopy(data)
    edges = out["edges"]
    children: dict[int, list[dict]] = {}
    incoming: dict[int, list[dict]] = {}
    for edge in edges:
        children.setdefault(int(edge["a"]), []).append(edge)
        incoming.setdefault(int(edge["b"]), []).append(edge)

    roots = [int(node["id"]) for node in out["nodes"] if not incoming.get(int(node["id"]))]
    queue = roots[:]
    seen = set()
    branch_count = 0
    unary_caps = 0
    while queue:
        node = queue.pop(0)
        if node in seen:
            continue
        seen.add(node)
        parent_edges = incoming.get(node, [])
        parent_radius = max((float(edge.get("radius", 0.0)) for edge in parent_edges), default=None)
        child_edges = children.get(node, [])
        if parent_radius and parent_radius > 0 and len(child_edges) == 1:
            child_radius = float(child_edges[0].get("radius", 0.0))
            if child_radius > parent_radius:
                child_edges[0]["radius"] = parent_radius * 0.98
                unary_caps += 1
        elif parent_radius and parent_radius > 0 and len(child_edges) >= 2:
            cubes = [max(float(edge.get("radius", 0.0)), 1e-6) ** 3 for edge in child_edges]
            scale = (parent_radius**3 / sum(cubes)) ** (1.0 / 3.0)
            for edge, cube in zip(child_edges, cubes):
                edge["radius"] = (cube ** (1.0 / 3.0)) * scale
            branch_count += 1
        queue.extend(int(edge["b"]) for edge in child_edges)

    return out, {"roots": roots, "branch_nodes_regularised": branch_count, "unary_radius_caps": unary_caps}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path("cco/data/cco_arterial15_repaired.json"))
    parser.add_argument("--output", type=Path, default=Path("cco/data/cco_arterial15_repaired_murray.json"))
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    smoothed, metadata = smooth(data)
    smoothed.setdefault("metadata", {})["radius_regularisation"] = {
        "method": "Murray cubic bifurcation scaling with unary-segment cap",
        **metadata,
    }
    args.output.write_text(json.dumps(smoothed, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()

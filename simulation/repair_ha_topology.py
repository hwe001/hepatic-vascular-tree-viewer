"""Attach disconnected HA components to the largest proximal component.

This does not modify the original provisional asset. Each added segment is
marked ``syntheticConnector`` with its source component and distance so the
repair remains auditable.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path


SOURCE = Path("cco/data/cco_arterial15_provisional.json")
OUTPUT = Path("cco/data/cco_arterial15_unregularised.json")


def components(data):
    outgoing = defaultdict(list)
    incoming = defaultdict(set)
    for edge in data["edges"]:
        outgoing[int(edge["a"])].append(edge)
        incoming[int(edge["b"])].add(int(edge["a"]))
    roots = [int(node["id"]) for node in data["nodes"] if not incoming[int(node["id"])] ]
    result = {}
    for root in roots:
        seen = set()
        stack = [root]
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            stack.extend(int(edge["b"]) for edge in outgoing.get(node, []))
        result[root] = seen
    return roots, result, outgoing


def main():
    data = json.loads(SOURCE.read_text(encoding="utf-8"))
    by_id = {int(node["id"]): node for node in data["nodes"]}
    roots, comps, outgoing = components(data)
    # Node 0 is the anatomical inlet: it has the largest proximal radius and
    # the short proximal trunk, whereas node 6 is a downstream component.
    main_root = 0
    connected_nodes = set(comps[main_root])
    pending_roots = set(roots) - {main_root}
    connectors = []

    while pending_roots:
        root = min(
            pending_roots,
            key=lambda candidate: min(
                math.dist(by_id[a]["p"], by_id[candidate]["p"])
                for a in connected_nodes
            ),
        )
        distance, main_node, component_node = min(
            (
                math.dist(by_id[a]["p"], by_id[root]["p"]),
                a,
                root,
            )
            for a in connected_nodes
        )
        first_edges = outgoing.get(root, [])
        if not first_edges:
            raise ValueError(f"Root {root} has no outgoing edge")
        radius = float(first_edges[0].get("radius", 0.3))
        by_id[main_node]["terminal"] = False
        connectors.append({
            "a": main_node,
            "b": root,
            "length": distance,
            "radius": radius,
            "resistance": 0.0,
            "flow": 0.0,
            "curve": [by_id[main_node]["p"], by_id[root]["p"]],
            "scaffold": False,
            "syntheticConnector": True,
            "sourceComponentRoot": root,
            "connectionDistanceMm": distance,
        })
        connected_nodes.update(comps[root])
        pending_roots.remove(root)

    repaired = dict(data)
    repaired["id"] = "cco_arterial15_repaired"
    repaired["label"] = "Hepatic-artery CCO tree with audited proximal topology repair"
    repaired["rootId"] = main_root
    repaired["edges"] = list(data["edges"]) + connectors
    repaired["assumptions"] = dict(data.get("assumptions", {}))
    repaired["assumptions"].update({
        "topologyRepair": "largest rooted component retained; five disconnected components attached by nearest-node synthetic connectors",
        "originalRootCandidates": roots,
        "primaryRoot": main_root,
        "syntheticConnectorCount": len(connectors),
        "syntheticConnectorDistancesMm": [c["connectionDistanceMm"] for c in connectors],
        "originalAssetPreserved": True,
    })
    repaired["summary"] = dict(data.get("summary", {}))
    repaired["summary"].update({
        "topology": "single connected rooted graph after explicit proximal repair",
        "syntheticConnectorCount": len(connectors),
    })
    OUTPUT.write_text(json.dumps(repaired, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "primaryRoot": main_root,
        "originalRoots": roots,
        "connectors": connectors,
    }, indent=2))


if __name__ == "__main__":
    main()

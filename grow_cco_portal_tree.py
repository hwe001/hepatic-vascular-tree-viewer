import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"


def flat_vertices(path):
    values = json.loads(path.read_text(encoding="utf-8"))["vertices"]
    return np.asarray(values, dtype=float).reshape(-1, 3)


def farthest(points, count):
    selected = [points[len(points) // 2]]
    while len(selected) < count:
        distances = np.min(((points[:, None] - np.asarray(selected)[None, :]) ** 2).sum(axis=2), axis=1)
        selected.append(points[np.argmax(distances)])
    return np.asarray(selected)


def nearest_parent(point, nodes):
    distances = ((nodes - point) ** 2).sum(axis=1)
    return int(np.argmin(distances))


portal = flat_vertices(DATA / "libzinc-json" / "portal15_1.json")
volume = json.loads((DATA / "portal_zones_volume.json").read_text(encoding="utf-8"))

# The existing mesh is used as a geometric scaffold. Its mesh vertices are not
# treated as a validated centreline; the added network is therefore synthetic.
root = portal.mean(axis=0)
scaffold = farthest(portal, 48)
nodes = [root]
parents = [-1]
segments = []

for point in sorted(scaffold, key=lambda p: float(np.linalg.norm(p - root))):
    parent = nearest_parent(point, np.asarray(nodes))
    nodes.append(point)
    parents.append(parent)
    segments.append({"a": nodes[parent].tolist(), "b": point.tolist(), "kind": "scaffold"})

# Add terminal perfusion targets from each current portal zone. Connecting the
# targets one at a time is the constructive part of this prototype; the
# objective is short, non-crossing-ish connections with balanced terminal load.
for zone_index, zone in enumerate(volume["zones"]):
    vertices = np.asarray(zone["vertices"], dtype=float).reshape(-1, 3)
    targets = farthest(vertices, 18)
    for target in sorted(targets, key=lambda p: float(np.linalg.norm(p - root))):
        parent = nearest_parent(target, np.asarray(nodes))
        nodes.append(target)
        parents.append(parent)
        segments.append({"a": nodes[parent].tolist(), "b": target.tolist(), "kind": "grown", "zone": zone_index + 1})

# Murray-style relative radius assignment from terminal load. Each new target
# contributes one unit of terminal demand; parent radii scale with downstream
# load. These are relative visualization radii, not calibrated portal radii.
children = [[] for _ in nodes]
for i in range(1, len(nodes)):
    children[parents[i]].append(i)
load = [1.0 if not children[i] else 0.0 for i in range(len(nodes))]
for i in range(len(nodes) - 1, -1, -1):
    if children[i]:
        load[i] = sum(load[j] for j in children[i])
edge_cursor = 0
for i in range(1, len(nodes)):
    radius = 0.22 * max(load[i], 1.0) ** 0.35
    segments[edge_cursor]["radius"] = round(float(radius), 4)
    edge_cursor += 1

(DATA / "cco_portal_growth.json").write_text(json.dumps({
    "method": "CCO-style synthetic extension from the portal15 geometric scaffold",
    "coordinate_context": "15 asset set",
    "scaffold_nodes": len(scaffold),
    "grown_terminals": 8 * 18,
    "segments": segments,
    "limitations": "The source portal JSON is a triangulated surface mesh without an annotated centreline. Added branches are synthetic visualization geometry, not a validated subject-specific reconstruction.",
}, separators=(",", ":")), encoding="utf-8")
print(f"Wrote {len(segments)} scaffold/grown segments to data/cco_portal_growth.json")

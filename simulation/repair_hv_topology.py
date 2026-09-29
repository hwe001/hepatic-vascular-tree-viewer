import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "cco" / "data" / "cco_hepatic_vein15_provisional.json"
OUTPUT = ROOT / "simulation" / "hv_topology_repaired.json"


def main():
    data = json.loads(GEOMETRY.read_text())
    nodes = {n["id"]: n for n in data["nodes"]}
    adjacency = defaultdict(set)
    for edge in data["edges"]:
        adjacency[edge["a"]].add(edge["b"])
        adjacency[edge["b"]].add(edge["a"])

    seen = set()
    components = []
    for root in nodes:
        if root in seen:
            continue
        stack = [root]
        seen.add(root)
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            for other in adjacency[current]:
                if other not in seen:
                    seen.add(other)
                    stack.append(other)
        components.append(component)

    origin = nodes[0]["p"]
    junctions = {
        "left": {"name": "Left hepatic vein", "node": 1, "color": "#7b1fa2"},
        "middle": {"name": "Middle hepatic vein", "node": 10, "color": "#00897b"},
        "right": {"name": "Right hepatic vein", "node": 255, "color": "#ef6c00"},
    }
    branch_vectors = {}
    for key, item in junctions.items():
        v = [nodes[item["node"]]["p"][i] - origin[i] for i in range(3)]
        norm = math.sqrt(sum(x * x for x in v))
        branch_vectors[key] = [x / norm for x in v]

    assignments = []
    virtual_links = []
    for component in sorted(components, key=lambda c: min(c)):
        anchor = min(component)
        if 0 in component:
            assignments.append({"component_anchor": anchor, "size": len(component), "branch": "ivc_scaffold", "virtual_link": None})
            continue
        centroid = [sum(nodes[node]["p"][i] for node in component) / len(component) for i in range(3)]
        direction = [centroid[i] - origin[i] for i in range(3)]
        norm = math.sqrt(sum(x * x for x in direction))
        direction = [x / norm for x in direction]
        branch = max(branch_vectors, key=lambda key: sum(direction[i] * branch_vectors[key][i] for i in range(3)))
        target = junctions[branch]["node"]
        length = math.dist(nodes[anchor]["p"], nodes[target]["p"])
        link = {"from_component_anchor": anchor, "to_junction_node": target, "branch": junctions[branch]["name"], "length": length, "inferred": True}
        assignments.append({"component_anchor": anchor, "size": len(component), "branch": branch, "virtual_link": link})
        virtual_links.append(link)

    result = {
        "model": "repaired three-hepatic-vein-to-IVC topology",
        "outlet": "IVC",
        "ivc_path_edges": [[58, 0], [0, 59], [59, 1331]],
        "junctions": [
            {"name": item["name"], "edge": [0, item["node"]], "color": item["color"]}
            for item in junctions.values()
        ],
        "component_assignments": assignments,
        "virtual_links": virtual_links,
        "flow_direction": "distal component anchors drain through inferred virtual links into one of three hepatic-vein junctions and then into the IVC",
        "status": "Derived topology repair; original geometry and measured scaffold are unchanged, and virtual links require anatomical validation",
    }
    OUTPUT.write_text(json.dumps(result, indent=2))
    print(f"wrote {OUTPUT} with {len(virtual_links)} inferred links")


if __name__ == "__main__":
    main()

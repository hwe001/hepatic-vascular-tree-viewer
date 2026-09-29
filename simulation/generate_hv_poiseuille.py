import json
import math
from collections import defaultdict, deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT / "cco" / "data" / "cco_hepatic_vein15_provisional.json"
ZONES = ROOT / "data" / "portal_zones_volume.json"
OUTPUT = ROOT / "simulation" / "hv_poiseuille_resection.json"

MU = 0.0035
INTACT_FLOW = 1050.0
REMNANT_FLOW = 950.0
RADIUS_SCALE = 3.1


def zone_centers():
    data = json.loads(ZONES.read_text())
    centers = []
    for zone in data["zones"]:
        v = zone["vertices"]
        n = len(v) // 3
        centers.append([sum(v[i::3]) / n for i in range(3)])
    return centers


def make_plane(centers):
    # Same left-versus-right partition used by the existing hepatectomy viewer.
    keep = centers[:4]
    remove = centers[4:]
    a = [sum(p[i] for p in keep) / len(keep) for i in range(3)]
    b = [sum(p[i] for p in remove) / len(remove) for i in range(3)]
    point = [(a[i] + b[i]) / 2 for i in range(3)]
    normal = [b[i] - a[i] for i in range(3)]
    norm = math.sqrt(sum(x * x for x in normal))
    return point, [x / norm for x in normal]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def components(nodes, edges):
    adjacency = defaultdict(list)
    for index, edge in enumerate(edges):
        adjacency[edge["a"]].append((edge["b"], index))
        adjacency[edge["b"]].append((edge["a"], index))
    seen = set()
    result = []
    for node in nodes:
        if node["id"] in seen:
            continue
        seen.add(node["id"])
        q = [node["id"]]
        ids = []
        while q:
            current = q.pop()
            ids.append(current)
            for other, _ in adjacency[current]:
                if other not in seen:
                    seen.add(other)
                    q.append(other)
        result.append((ids, adjacency))
    return result, adjacency


def build_scenario(nodes, edges, adjacency, plane, total_flow, remnant):
    node_by_id = {n["id"]: n for n in nodes}
    retained = {}
    for node in nodes:
        if not remnant:
            retained[node["id"]] = True
        else:
            p, normal = plane
            # The HV asset uses the opposite viewing orientation from the PV
            # and HA assets: the right-lobe drainage occupies the negative
            # side of this partition. Keep the positive side for the remnant.
            retained[node["id"]] = dot([x - y for x, y in zip(node["p"], p)], normal) >= 0

    components_data, _ = components(nodes, edges)
    edge_rows = []
    total_sites = 0
    component_models = []

    for ids, _ in components_data:
        id_set = set(ids)
        scaffold = [x for x in ids if node_by_id[x].get("scaffold")]
        anchor = 0 if 0 in id_set else (min(scaffold) if scaffold else min(ids))
        parent = {anchor: None}
        parent_edge = {}
        order = [anchor]
        for current in order:
            for other, edge_index in adjacency[current]:
                if other in id_set and other not in parent:
                    parent[other] = current
                    parent_edge[other] = edge_index
                    order.append(other)

        children = defaultdict(list)
        for child, par in parent.items():
            if par is not None:
                children[par].append(child)

        site_count = {}
        for current in reversed(order):
            kids = children[current]
            is_site = bool(node_by_id[current].get("terminal")) or (not kids and current != anchor)
            site_count[current] = (1 if is_site and retained[current] else 0) + sum(site_count[k] for k in kids)
        component_models.append((anchor, parent, parent_edge, site_count))
        total_sites += site_count[anchor]

    for anchor, parent, parent_edge, site_count in component_models:
        for child, edge_index in parent_edge.items():
            edge = edges[edge_index]
            a = node_by_id[edge["a"]]
            b = node_by_id[edge["b"]]
            visible = retained[a["id"]] and retained[b["id"]] and site_count[child] > 0
            q = total_flow * site_count[child] / total_sites if total_sites else 0.0
            # The legacy CCO radius field is rendered with this scale in the
            # existing viewer; use the same conversion for the reduced-order
            # calculation while preserving the source geometry unchanged.
            radius_mm = max(float(edge.get("radius") or 0.18) * RADIUS_SCALE, 0.05)
            length_mm = max(float(edge.get("length") or 0.0), 1.0e-6)
            q_m3_s = q * 1.0e-6 / 60.0
            radius_m = radius_mm * 1.0e-3
            resistance = 8.0 * MU * (length_mm * 1.0e-3) / (math.pi * radius_m**4)
            wss = 4.0 * MU * q_m3_s / (math.pi * radius_m**3)
            edge_rows.append({
                "a": edge["a"], "b": edge["b"], "q_ml_min": q,
                "wss_pa": wss, "pressure_drop_mmhg": resistance * q_m3_s / 133.322,
                "visible": visible, "radius_mm": radius_mm,
                "length_mm": length_mm, "component_anchor": anchor,
            })

    visible_rows = [r for r in edge_rows if r["visible"]]
    drop_by_key = {(r["a"], r["b"]): r["pressure_drop_mmhg"] for r in edge_rows}
    max_path_drop = 0.0
    for anchor, parent, parent_edge, _ in component_models:
        cumulative = {anchor: 0.0}
        for child, par in parent.items():
            if par is None:
                continue
            e = edges[parent_edge[child]]
            cumulative[child] = cumulative[par] + drop_by_key[(e["a"], e["b"])]
            if retained[child]:
                max_path_drop = max(max_path_drop, cumulative[child])
    return {
        "total_flow_ml_min": total_flow,
        "visible_edge_count": len(visible_rows),
        "all_edge_count": len(edge_rows),
        "terminal_sites_retained": total_sites,
        "edge_values": edge_rows,
        "wss_pa": {"min": min((r["wss_pa"] for r in visible_rows), default=0), "max": max((r["wss_pa"] for r in visible_rows), default=0)},
        "pressure_drop_mmhg": {"max_path": max_path_drop, "sum_segments": sum(r["pressure_drop_mmhg"] for r in visible_rows)},
    }


def main():
    geometry = json.loads(GEOMETRY.read_text())
    nodes, edges = geometry["nodes"], geometry["edges"]
    plane = make_plane(zone_centers())
    _, adjacency = components(nodes, edges)
    report = {
        "model": "steady hepatic-venous Poiseuille diagnostic",
        "assumptions": {
            "dynamic_viscosity_pa_s": MU,
            "radius_conversion_factor": RADIUS_SCALE,
            "intact_root_outflow_ml_min": INTACT_FLOW,
            "remnant_root_outflow_ml_min": REMNANT_FLOW,
            "flow_allocation": "prescribed total outflow distributed by retained terminal-site fraction in each connected component",
            "resection": "standard right hepatectomy plane; the HV-specific negative half-space is hidden because this asset is oppositely oriented to the PV and HA assets",
            "caveat": "HV geometry is provisional and its anatomical outlet annotation remains pending; values are scenario outputs, not subject-specific measurements",
        },
        "scenarios": {
            "intact_pre_hepatectomy": build_scenario(nodes, edges, adjacency, plane, INTACT_FLOW, False),
            "remnant_post_right_hepatectomy": build_scenario(nodes, edges, adjacency, plane, REMNANT_FLOW, True),
        },
    }
    OUTPUT.write_text(json.dumps(report, indent=2))
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()

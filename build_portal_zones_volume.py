import json
from pathlib import Path

import numpy as np
import trimesh
from skimage.measure import marching_cubes

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "libzinc-json"
OUT = ROOT / "data" / "portal_zones_volume.json"


def load_legacy_geometry(path):
    data = json.loads(path.read_text(encoding="utf-8"))
    vertices = [tuple(data["vertices"][i:i + 3]) for i in range(0, len(data["vertices"]), 3)]
    raw = data["faces"]
    faces = []
    i = 0
    while i < len(raw):
        flags = raw[i]
        i += 1
        count = 4 if flags & 1 else 3
        face = tuple(raw[i:i + count])
        i += count
        if flags & 2: i += 1
        if flags & 4: i += 1
        if flags & 8: i += count
        if flags & 16: i += 1
        if flags & 32: i += count
        if flags & 64: i += 1
        if flags & 128: i += count
        if len(face) == 3:
            faces.append(face)
        else:
            faces.extend(((face[0], face[1], face[3]), (face[1], face[2], face[3])))
    return np.asarray(vertices, dtype=float), np.asarray(faces, dtype=int)


def farthest_seeds(points, count=8):
    seeds = [points[len(points) // 2]]
    while len(seeds) < count:
        distances = np.min(((points[:, None, :] - np.asarray(seeds)[None, :, :]) ** 2).sum(axis=2), axis=1)
        seeds.append(points[np.argmax(distances)])
    return np.asarray(seeds)


surface_vertices, surface_faces = load_legacy_geometry(DATA / "surface15_1.json")
portal_vertices, _ = load_legacy_geometry(DATA / "portal15_1.json")
surface_mesh = trimesh.Trimesh(surface_vertices, surface_faces, process=False)
pitch = 4.0
occupied = surface_mesh.voxelized(pitch=pitch).fill()
points = np.asarray(occupied.points)
seeds = farthest_seeds(portal_vertices)
labels = ((points[:, None, :] - seeds[None, :, :]) ** 2).sum(axis=2).argmin(axis=1)

zones = []
for zone in range(8):
    mask = np.zeros(occupied.matrix.shape, dtype=np.float32)
    mask[occupied.sparse_indices[labels == zone, 0], occupied.sparse_indices[labels == zone, 1], occupied.sparse_indices[labels == zone, 2]] = 1.0
    vertices, faces, _, _ = marching_cubes(mask, level=.5, spacing=(pitch, pitch, pitch), step_size=1)
    vertices += occupied.transform[:3, 3]
    zones.append({
        "label": zone + 1,
        "vertices": vertices.reshape(-1).round(4).tolist(),
        "faces": faces.reshape(-1).tolist(),
    })

OUT.write_text(json.dumps({
    "method": "liver-surface voxel occupancy with nearest portal-seed assignment",
    "voxel_pitch": pitch,
    "coordinate_context": "15 asset set",
    "zones": zones,
}, separators=(",", ":")), encoding="utf-8")
print(f"Wrote {OUT} with {len(points)} liver-interior voxels")

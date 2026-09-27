# Matched hepatic vascular-tree viewer

Separate Three.js viewer for inspecting three hepatic blood-vessel trees, the biliary tree, and the liver surface:

- hepatic artery: `data/libzinc-json/arterial15_1.json`
- portal vein: `data/libzinc-json/portal15_1.json`
- hepatic vein: `data/libzinc-json/hepatic15_1.json`
- biliary tree: `data/libzinc-json/bile15_1.json`
- liver surface: `data/libzinc-json/surface15_1.json`

The browser viewer is in `viewer/index.html`. It provides independent visibility and opacity controls for the liver surface and each vascular tree, together with rotation, right-drag panning, zooming, wireframe inspection, and a fit-view control.

A dedicated portal-segmentation viewer is available at `segments/index.html` and online at https://hwe001.github.io/hepatic-vascular-tree-viewer/segments/. It shows only the portal vein, a faint liver-surface reference, and eight independently generated colored 3D portal-territory surfaces. The surfaces are reconstructed from a voxelized occupancy mask derived from the actual liver surface rather than painted onto the liver-surface triangles or forced into an ellipsoid.

A CCO-output viewer is available at `cco/index.html`. It loads the JSON output of the separate `hwe001/liver-portal-cco` pipeline and displays the generated curved portal-tree centreline with the matched liver-surface reference. The public repository contains the viewer adapter and input/output documentation; the patient-specific generator and its required segment inputs remain in the separate private repository.


It also includes an exploratory eight-zone portal-territory overlay. The overlay derives eight representative regions from the portal-vein mesh and assigns each liver-surface vertex to its nearest portal region. The displayed I-VIII labels are provisional Couinaud labels: the source JSON does not provide explicit portal branch annotations or radiological landmarks, so this overlay is intended for hypothesis generation and visualization, not clinical segmentation.

## Provenance and scope

The assets were identified in the public repository [`hwe001/hyu754.github.io`](https://github.com/hwe001/hyu754.github.io), under `liverModelPatient1/`. The `15` asset set was selected because the arterial, portal, hepatic-vein, biliary, and surface files share a common coordinate context in the source repository. The files are retained in their original LibZinc/Three.js-compatible JSON format.

This is a geometry and visualization resource. It does not contain pressure, flow, velocity, or calibrated physiological predictions, and no vascular-flow simulation is implied by the overlay. The displayed biliary tree and vascular trees are matched geometry assets from the same `15` coordinate context; users should still treat the combined display as an anatomical registration resource rather than a validated subject-specific reconstruction.

To run locally, serve the repository directory over HTTP and open `viewer/index.html`. A live viewer is available at https://hwe001.github.io/hepatic-vascular-tree-viewer/viewer/ after GitHub Pages is enabled.

## Suggested citation

Ho H, Bartlett A, Jalu M. Matched hepatic vascular-tree geometry and Three.js viewer. Public data repository. Release 1.0.

## License

See `LICENSE`. Users should also review the provenance and licensing terms of the original source repository before redistributing modified derivatives.

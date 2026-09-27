# Matched hepatic vascular-tree viewer

Separate Three.js viewer for inspecting three hepatic blood-vessel trees, the biliary tree, and the liver surface:

- hepatic artery: `data/libzinc-json/arterial15_1.json`
- portal vein: `data/libzinc-json/portal15_1.json`
- hepatic vein: `data/libzinc-json/hepatic15_1.json`
- biliary tree: `data/libzinc-json/S01_bile_1.json`
- liver surface: `data/libzinc-json/surface15_1.json`

The browser viewer is in `viewer/index.html`. It provides independent visibility and opacity controls for the liver surface and each vascular tree, together with rotation, right-drag panning, zooming, wireframe inspection, and a fit-view control.

## Provenance and scope

The vascular and liver-surface assets were identified in the public repository [`hwe001/hyu754.github.io`](https://github.com/hwe001/hyu754.github.io), under `liverModelPatient1/`. The `15` asset set was selected because the arterial, portal, hepatic-vein, and surface files share a common coordinate context in the source repository. The biliary tree is included from the companion biliary-geometry resource and is shown as an independently toggleable geometry layer. The files are retained in their original LibZinc/Three.js-compatible JSON format.

This is a geometry and visualization resource. It does not contain pressure, flow, velocity, or calibrated physiological predictions, and no vascular-flow simulation is implied by the overlay. Because the biliary asset comes from a separate source geometry package, users should treat the combined display as an exploratory overlay unless they independently verify registration between the biliary and vascular coordinate systems.

To run locally, serve the repository directory over HTTP and open `viewer/index.html`. A live viewer is available at https://hwe001.github.io/hepatic-vascular-tree-viewer/viewer/ after GitHub Pages is enabled.

## Suggested citation

Ho H, Bartlett A, Jalu M. Matched hepatic vascular-tree geometry and Three.js viewer. Public data repository. Release 1.0.

## License

See `LICENSE`. Users should also review the provenance and licensing terms of the original source repository before redistributing modified derivatives.

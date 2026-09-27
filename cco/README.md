# CCO portal-tree viewer

This page is the public visualization adapter for the private `hwe001/liver-portal-cco` generator. It loads the generator's `cco_post_s*.json` output and displays its curved centreline tree against the matched liver-segment reference.

The bundled result is grown from the public `portal15_1.json` asset: 8,254 nodes, 8,253 edges, and 4,129 terminals across 32 generations. A controlled subset of 809 leaves that were already within 8 mm of the surface is locally anchored within 0.75 mm after validated interior-segment checks; distant interior leaves are left untouched to avoid radial spokes. The current liver surface asset is converted to STL as the containment surface. It is loaded automatically when the page opens.

The generator is not copied into this public repository because it is maintained separately and expects patient-specific inputs that are not public here. The private pipeline uses territory-partitioned growth, coverage-driven adaptive growth, Murray-law radius scaling, clearance checks, concave-boundary repair, long-tail densification, and a forward pressure/flow solve.

## Volume-filling comparison viewer

`cco-volume/index.html` loads `cco/data/cco_post_s2_volume_filling.json`, a separate CCO-inspired comparison result. It uses one connected global tree grown from the observed portal root, uniform interior perfusion targets, coverage-driven gap filling, and no forced surface-terminal shell. This is the current morphology experiment intended to reduce radial spokes and improve volume filling. The original surface-anchored result remains available at `cco/index.html` for side-by-side comparison.

The deep-learning/LFDO-CCO repository at [`pqpqpqpqpq/Generation`](https://github.com/pqpqpqpqpq/Generation) is a useful future extension for learning local bifurcation statistics from data. It is not used to claim a trained model for this subject-specific tree.

## Use

1. Run the private CCO pipeline with its `portal_post.json` and segment STL inputs.
2. Open `cco/index.html` over HTTP.
3. Choose the generated `cco_post_s*.json` file.

The viewer also supports mobile interaction: drag to rotate, two-finger/right-drag to pan where supported, and pinch/wheel to zoom.

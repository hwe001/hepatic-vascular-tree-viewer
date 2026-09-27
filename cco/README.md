# CCO portal-tree viewer

This page is the public visualization adapter for the private `hwe001/liver-portal-cco` generator. It loads the generator's `cco_post_s*.json` output and displays its curved centreline tree against the matched liver-segment reference.

The bundled result is grown from the public `portal15_1.json` asset: 8,254 nodes, 8,253 edges, and 4,129 terminals across 32 generations. A controlled subset of 1,032 distal endpoints is anchored within 0.75 mm of the liver surface after validated interior-segment checks. The current liver surface asset is converted to STL as the containment surface. It is loaded automatically when the page opens.

The generator is not copied into this public repository because it is maintained separately and expects patient-specific inputs that are not public here. The private pipeline uses territory-partitioned growth, coverage-driven adaptive growth, Murray-law radius scaling, clearance checks, concave-boundary repair, long-tail densification, and a forward pressure/flow solve.

## Use

1. Run the private CCO pipeline with its `portal_post.json` and segment STL inputs.
2. Open `cco/index.html` over HTTP.
3. Choose the generated `cco_post_s*.json` file.

The viewer also supports mobile interaction: drag to rotate, two-finger/right-drag to pan where supported, and pinch/wheel to zoom.

# CCO portal-tree viewer

This page is the public visualization adapter for the private `hwe001/liver-portal-cco` generator. It loads the generator's `cco_post_s*.json` output and displays its curved centreline tree against the matched liver-segment reference.

The bundled demonstration is the `demo-donor1-anatomy` CCO result: 9,766 nodes, 9,765 edges, and 4,885 terminals in segment S2. It is loaded automatically when the page opens.

The generator is not copied into this public repository because it is maintained separately and expects patient-specific inputs that are not public here. The private pipeline uses territory-partitioned growth, coverage-driven adaptive growth, Murray-law radius scaling, clearance checks, concave-boundary repair, long-tail densification, and a forward pressure/flow solve.

## Use

1. Run the private CCO pipeline with its `portal_post.json` and segment STL inputs.
2. Open `cco/index.html` over HTTP.
3. Choose the generated `cco_post_s*.json` file.

The viewer also supports mobile interaction: drag to rotate, two-finger/right-drag to pan where supported, and pinch/wheel to zoom.

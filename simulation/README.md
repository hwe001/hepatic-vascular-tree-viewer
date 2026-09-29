# Python hepatic circulation controller

`habr_0d.py` is a standalone Python translation of the dynamic controller
used in the archived hepatic circulation model:

- `HABR/runmodel.m`
- `HABR/HABRresistorsdebbaut.mdl`

The implementation keeps the original pressure-source values, resistances,
PV/HA inertances, hepatic-artery pulse, unit conversion, and quadratic HABR
response. The three reported flows are:

```text
PV = I_pv
HA = I_ha
HV = I_pv + I_ha
```

The last relation is the reduced-model conservation closure used for the
integrated viewer. The archived Simulink model also contains RC branches and
initial-state loading. Their parameter values are preserved in the Python
dataclass as provenance, but they are not silently represented as an exact
Simscape export. A direct trajectory comparison should be completed before
describing this code as a bit-for-bit replacement.

Run the archived-controller translation from the repository root:

```powershell
python simulation/habr_0d.py
```

The script prints the baseline mean flow, the HABR-updated arterial flow and
resistances, and the post-change mean flow. It requires NumPy and SciPy.

Validate the scalar outputs against the archived MATLAB `.mat` result:

```powershell
python simulation/validate_habr.py
```

This writes `habr_validation_report.json`. The current archive contains
controller summary values rather than full Simulink trajectories, so the
validator reports exact scalar agreement and explicitly marks trajectory
comparison as unavailable.

For transient exploration with the archived RC values, use
`simulate_with_compliance(...)` from `habr_0d.py`. This is a documented
topology reconstruction and currently a separate exploratory path; its
transient traces must be compared with exported Simulink signals before it
is used as the paper's quantitative result.

## Hepatic-artery transmission line

`ha_transmission_line.py` provides a linearized frequency-domain 1-D
transmission-line solve on the Murray-regularised
`cco/data/cco_arterial15_repaired_murray.json` asset. It
propagates a zero-mean harmonic perturbation around a root mean flow of
250 mL/min and reports segment flow amplitude, phase, pressure amplitude, and
downstream flow amplitude. Segment length and radius come from the geometry;
wall modulus, wall-thickness ratio, blood properties, and terminal loads are
explicit prototype assumptions in `Parameters` and are written into the JSON
result. The result is a reproducible mechanics layer, not yet a
patient-calibrated arterial solution.

The original provisional HA asset contained six disconnected root components.
`simulation/repair_ha_topology.py` creates a separate repaired asset with one
root (node 6), preserving the original segments and adding five explicitly
marked proximal connectors. Connector lengths and source component roots are
recorded in the repaired JSON metadata; the original asset is not modified.

Run it from the repository root:

```powershell
python simulation/ha_transmission_line.py
```

## HA HABR and WSS comparison

`ha_habr_wss_comparison.py` compares the intact HA tree with the retained
Z2+Z3+Z4 remnant after standard right hepatectomy. It uses conservative
Murray-style radius-cubed branch splits, a documented blood viscosity, and a
flow-consistent effective radius only for WSS post-processing. The generated
JSON and edge CSV record the assumptions and scenario summaries:

```powershell
python simulation/ha_habr_wss_comparison.py
```

The interactive result is available locally at `flow-ha-wss/index.html` and,
after GitHub Pages deployment, at
`https://hwe001.github.io/hepatic-vascular-tree-viewer/flow-ha-wss/`. The HABR
case is an illustrative sensitivity scenario, not a calibrated postoperative
prediction or a substitute for 3-D CFD.

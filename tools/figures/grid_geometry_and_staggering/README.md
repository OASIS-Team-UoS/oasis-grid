# Horizontal geometry and variable staggering

Publication figure for *Numerical Models of Planetary Atmosphere Dynamics*.
All deliverables, source, configuration and mesh provenance are in this directory.
The composite is 180 × 164 mm, with 8.5–10 pt DejaVu Sans text. The independent
geometry and staggering exports are 180 × 80 mm and 180 × 78 mm, respectively.

## Generate from the repository root

Use Python 3.10+ and the repository's existing dependencies:

```bash
python3 -m pip install -r requirements.txt
python3 tools/figures/grid_geometry_and_staggering/generate.py --reuse-native
```

This command uses the included, checksum-verified native snapshot and needs
neither an OASIS checkout nor a compiler. To refresh the native mesh from
source, use a C++17 compiler (`c++` by default) and an OASIS checkout,
discovered by default at `../oasis`:

```bash
python3 tools/figures/grid_geometry_and_staggering/generate.py
```

Other checkout locations work with:

```bash
python3 tools/figures/grid_geometry_and_staggering/generate.py --oasis-checkout /path/to/oasis
```

No OASIS build, MPI launch, climate integration or HPC job is needed. The utility
compiles three geometry methods in a temporary directory. Set `TMPDIR` to a
scratch directory if the system temporary filesystem is short of space.

Edit `config.json` for fonts, line weights, sphere/schematic sizes, projection
camera, angular sampling, densities and export dimensions. Page layout uses
millimetres in `generate.py`. Use `--config path.json` or `--output-dir path`
for variants. The default command regenerates only its own registered outputs;
it refuses to replace an unregistered file or an output modified since its last
generation. Exports are staged before replacement. Existing repository edits
and model defaults are untouched.

## Files

- `grid_geometry_and_staggering.pdf`, `.svg`, `.png`: main composite.
- `horizontal_grid_geometries.pdf`, `.svg`, `.png`: geometry group.
- `horizontal_variable_staggering.pdf`, `.svg`, `.png`: staggering group and shared legend.
- `caption.tex`: manuscript caption and label, with no invented citation keys.
- `generate.py`, `native_mesh.py`, `config.json`: reproducible workflow.
- `metadata.json`, `generation_config.json`: revisions, working-tree status,
  source hashes, dependencies, camera, mesh parameters and construction checks.
- `native_mesh.npz`, `native_provenance.json`: exact native mesh snapshot and provenance.
- `native_mesh_extract.cpp`: unchanged native method bodies plus the mesh-only wrapper.
  This generated source retains OASIS's GPL-3.0-or-later notice; the license
  text is in `LICENSE.native`.
- `generated_files.json`: SHA-256 ledger protecting generated outputs.
- `audit_exports.py`, `requirements-audit.txt`, `export_audit.json`: optional
  export rendering/inspection workflow and the report for the delivered files.

PDF and SVG mesh edges are vector paths; PDF fonts are embedded and SVG text
remains editable. All PNGs are 600 dpi. No raster mesh image is embedded in
the vector exports. SVG editors should have DejaVu Sans installed.

The publication exports and native snapshot are intentionally versioned here.
Repository metadata records revision, working-tree status counts and changes
to the figure inputs, without listing unrelated files elsewhere in a checkout.

## Mesh provenance

| Panel | Generator and parameters | Cells | Vertices | Edges |
| --- | --- | ---: | ---: | ---: |
| (a) | Figure reference generator: 24 longitude × 12 latitude intervals | 288 | 266 | 552 |
| (b) | Existing `oasis_grid.build_grid_model`, equiangular, level 2; 4 × 4 per face | 96 | 98 | 192 |
| (c) | Native OASIS methods, level 2, spring disabled | 162 | 320 | 480 |

Counts cover the full sphere. Density was selected for readable examples, with
no assertion of equal effective resolution. Latitude parallels are genuine
small circles, not great-circle interpolants. The cubed-sphere mapping is the
explorer's equiangular gnomonic construction; it is not claimed as an active
OASIS dynamical backend. Inspection of `../oasis/docs/experimental_grid_abstraction.md`
confirmed a read-only icosahedral adapter and a cubed-sphere placeholder;
`../oasis/experimental/native_grid` was absent in this checkout.

The explorer's icosahedral rhombi describe the construction mesh, not the
finite-volume polygons. `native_mesh.py` extracts and executes these actual
native methods, without changing their bodies:

- `src/grid/grid_build.cpp`: `Icogrid::sphere_ico`.
- `src/grid/grid_pointsq.cpp`: `Icogrid::find_qpoints_thor` and
  `Icogrid::relocate_centres_thor`.

Ten-panel coordinates/connectivity are globally deduplicated. Each rhombus
uses the native (+1,+1) triangular diagonal; native and Python centre
coordinates are matched bijectively. Clockwise neighbor rings are assembled
from that connectivity, and checked against the native pentagon indices.
The native q-point routine computes spherical triangle area-centroid vertices.
These vertices, shared by construction-triangle identity, form the polygonal
dual. Native centre relocation follows q-point construction; the q-points are
not recomputed afterwards. No generic icosphere or Voronoi library is used.

The chosen mesh has no spring relaxation. Native spring algorithms exist,
but this workflow does not invoke them or claim to show a production relaxed
grid. The full OASIS constructor requires level ≥ 4; the verified coarse
level-2 example calls its mesh-only methods directly. It is not a supported
full-model run configuration. The delivered native snapshot preserves these
choices, including relocated centres, construction triangles and neighbor rings.

The camera is longitude 22°, latitude 24°, with an equal sphere diameter of
48.78 mm. Circle/limb intersections are solved analytically before vector
sampling at 0.5°. Only visible intervals are drawn. Each globally shared edge
is drawn once, including cube-face seams. Pentagon boundaries use a heavier
neutral stroke, including visible parts of limb-crossing cells.

## Validation and inspection

Generation checks finite unit coordinates; valid indices; unique vertices;
reciprocal neighbors; consistent outward face winding; opposite shared-edge
orientations; two cells per edge; connectedness; Euler characteristic 2; and
4π surface coverage (absolute tolerance 1e-10 steradians). Latitude-longitude
areas use the exact small-circle formula. Other polygons use oriented
spherical-triangle solid angles. The native dual has twelve pentagons and
150 hexagons; six pentagons have visible boundary portions for this camera.

Native/Python coordinates agree within 1.3e-16, shared native q-points within
4.8e-16, and an independent spherical area-moment formula checks native
corners and relocated centres within 2.3e-15. Relocated centres lie inside
their native polygons. Actual residuals and tolerances are in `metadata.json`.
These checks validate figure construction, not solver accuracy or conservation.

Marker checks require identical A-grid anchors, correct C-grid face/component
positions, and a shared Z-grid scalar/divergence/vorticity location with no
independently stored velocity arrows. The drawing uses these checked arrays.
Additional clipping cases cover hidden edges, straddling edges, a visible
middle arc with hidden endpoints, and small-circle parallels. Text overlap,
minimum font size, crop boundaries and absence of SVG raster images are checked.

The existing repository tests and compile checks were run:

```bash
python3 -m pytest
python3 -m py_compile oasis_grid.py oasis-grid tools/figures/grid_geometry_and_staggering/generate.py tools/figures/grid_geometry_and_staggering/native_mesh.py
```

Seven tests passed. To independently rasterize every delivered format at
144 dpi (180 mm print width), check PDF images/embedded fonts, SVG editable
text, PNG dimensions/dpi, grayscale rendering and page margins:

```bash
python3 -m pip install -r tools/figures/grid_geometry_and_staggering/requirements-audit.txt
python3 tools/figures/grid_geometry_and_staggering/audit_exports.py > tools/figures/grid_geometry_and_staggering/export_audit.json
```

The report binds the nine inspected files by SHA-256 and records the temporary
inspection-image directory. Re-run the audit after regenerating the figure.
All nine exports were rendered and visually inspected at the intended print
width. Layout and hidden-surface checks passed. The groups use neutral styling;
the schematics are independent of the spheres above them.

## Literature context

No manuscript bibliography was supplied or found in the inspected repositories,
so the caption adds no citation keys. Relevant primary sources are
[Mendonça et al. (2016), THOR](https://arxiv.org/abs/1607.05535),
doi:10.3847/0004-637X/829/2/115, and
[Eldred & Randall (2017)](https://gmd.copernicus.org/articles/10/791/2017/),
doi:10.5194/gmd-10-791-2017. Code provenance, rather than an inferred
literature construction, determines the OASIS geometry plotted here.

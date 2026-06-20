# OASIS grid explorer

Python-only tools to build, inspect, and postprocess spherical grids for OASIS
3D-model development.

The icosahedral builder follows the ten-panel rhombic layout used by OASIS.
The cube-sphere builder provides an equivalent six-panel grid so both layouts
can be explored through the same domain and MPI postprocessing functions.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

For development and tests:

```bash
pip install -e ".[dev]"
python3 -m pytest
```

## Run the explorer

```bash
./oasis-grid
```

Open the cube-sphere directly:

```bash
./oasis-grid --grid cube-sphere --glevel 4 --divide-face 2
```

Save a non-interactive snapshot:

```bash
./oasis-grid --grid icosahedral --save grid.png --no-show
```

Use a more transparent grid in a saved view:

```bash
./oasis-grid --surface-alpha 0.25 --save transparent_grid.png --no-show
```

Start zoomed in:

```bash
./oasis-grid --zoom 1.6 --save zoomed_grid.png --no-show
```

## Postprocess MPI domains

Print a grid summary and MPI-domain table:

```bash
./oasis-grid --postprocess --grid icosahedral --glevel 4 --divide-face 2
```

Inspect one MPI rank/domain:

```bash
./oasis-grid --domain-info 12 --grid icosahedral --glevel 4 --divide-face 2
```

Export geometry and domain arrays:

```bash
./oasis-grid --export-npz oasis_ico_grid.npz --grid icosahedral --glevel 4 --divide-face 2
```

## Python API

```python
from oasis_grid import (
    build_grid_model,
    generate_oasis_icosahedral_points,
    postprocess_grid,
    postprocess_domain,
    postprocess_mpi_domains,
    extract_domain_mesh,
)

ico = build_grid_model("icosahedral", glevel=4, divide_face=2)
cube = build_grid_model("cube-sphere", glevel=4, divide_face=2)

summary = postprocess_grid(ico)
rank_12 = postprocess_domain(ico, 12)
mesh_12 = extract_domain_mesh(ico, 12)
mpi = postprocess_mpi_domains(ico)
xyz, rhomb, panel_grids = generate_oasis_icosahedral_points(4, return_coarse=True)
```

The default MPI rank order is panel-major, then `tile_y`, then `tile_x`.
Each `GridModel` exposes `points`, `panel_grids`, `panel_names`, and `domains`.

## Interactive controls

- Grid selector: switch between icosahedral and cube-sphere grids.
- Color selector: color by domain, base panel, or height.
- `level`: panel edge resolution, where each panel has `2**level` cells per edge.
- `split`: domain decomposition, where each panel has `4**split` domains.
- `domain`: focus one MPI domain/rank.
- `alpha`: grid/surface opacity, from transparent to opaque. With `surface`
  enabled and `alpha` at `1.0`, back-side grid lines, points, and labels are
  hidden behind the opaque sphere.
- `zoom`: zoom out/in without changing the grid data.
- Toggles: surface, grid lines, points, and labels.
- Keyboard: `n`/right arrow for next domain, `p`/left arrow for previous domain,
  `+`/`-` for zoom, `a`/Escape for all domains, and `g` to switch grid type.

## License

This project is released under the BSD-3-Clause license. See `LICENSE`.

## Contributing and citation

See `CONTRIBUTING.md` for development guidance and `CITATION.cff` for citation
metadata. Before publishing or mirroring the repository, review
`RELEASE_CHECKLIST.md`.

# Contributing

Contributions are welcome. Keep changes focused and include a short explanation
of the grid behavior or postprocessing workflow being changed.

## Development setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

## Checks

Run the lightweight test suite before opening a pull request:

```bash
python3 -m pytest
python3 -m py_compile oasis_grid.py oasis-grid
```

For visualization changes, also save at least one non-interactive render:

```bash
./oasis-grid --grid icosahedral --glevel 3 --divide-face 1 --save smoke.png --no-show
```

Generated `.png` and `.npz` files are ignored by git.

## Pull requests

- Keep numerical grid-generation changes separate from UI-only changes when possible.
- Mention which grid family is affected: icosahedral, cube-sphere, or both.
- Include test coverage for public API changes.
- For changes to OASIS 3D-model conventions, document the MPI-domain/rank ordering.

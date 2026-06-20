# Release checklist

Use this checklist before pushing a public release or mirroring to a group
account.

## Repository metadata

- Confirm the canonical GitHub repository URL.
- Add the repository URL to `pyproject.toml` under `[project.urls]`.
- Add the repository URL to `CITATION.cff` as `repository-code`.
- Confirm the BSD-3-Clause license and copyright holder are acceptable.
- Confirm the README wording is acceptable for public OASIS-related work.

## Local checks

```bash
python3 -m pytest
python3 -m py_compile oasis_grid.py oasis-grid
./oasis-grid --summary --grid icosahedral --glevel 3 --divide-face 1
./oasis-grid --domain-info 7 --grid cube-sphere --glevel 3 --divide-face 1
```

## Git setup

```bash
git init
git add .
git commit -m "Initial OASIS grid explorer"
git branch -M main
```

Add remotes:

```bash
git remote add origin git@github.com:<your-user>/oasis-grid.git
git remote add group git@github.com:<group-account>/oasis-grid.git
```

Push both:

```bash
git push -u origin main
git push -u group main
```

If the group repository should be the canonical one, swap the names so `origin`
points to the group repository and use a name such as `personal` for your fork.

# Release checklist

Use this checklist before pushing a public release or mirroring to a group
account.

## Repository metadata

- Canonical repository: `https://github.com/jmmendonca/oasis-grid`
- Group mirror: `https://github.com/OASIS-Team-UoS/oasis-grid`
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
git remote add origin https://github.com/jmmendonca/oasis-grid.git
git remote add mirror https://github.com/OASIS-Team-UoS/oasis-grid.git
```

Push both:

```bash
git push -u origin main
git push -u mirror main
```

Use `origin` for the personal canonical repository and `mirror` for the group
repository unless the project governance changes.

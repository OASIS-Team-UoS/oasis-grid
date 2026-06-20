"""OASIS spherical grid builder, explorer, and postprocessing tools.

This file is intentionally self-contained so it can be used both as a small
interactive application and as a reference module for OASIS 3D-model grid work.
It contains three layers:

1. Grid generation
   The icosahedral path follows the OASIS ten-panel rhombic construction. The
   cube-sphere path creates a six-panel equivalent projected to the unit sphere.

2. Postprocessing helpers
   Functions such as ``postprocess_domain`` and ``extract_domain_mesh`` expose
   MPI-domain/rank metadata, domain boundaries, neighbors, cell corner indices,
   and normalized cell centers.

3. Interactive exploration
   ``GridExplorer`` renders both grid families with controls for resolution,
   MPI-domain split, focus domain, opacity, zoom, labels, points, and linework.

Usage:
    ./oasis-grid
    ./oasis-grid --grid cube-sphere --glevel 4 --divide-face 2
    ./oasis-grid --summary --no-show
    ./oasis-grid --domain-info 12 --no-show
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401, needed to register 3D projection
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
from matplotlib.widgets import Button, CheckButtons, RadioButtons, Slider


# ---------------------------------------------------------------------------
# Public convenience entry points
# ---------------------------------------------------------------------------

def oasis_grid(
    glevel: int = 4,
    divide_face: int = 2,
    grid: str = "icosahedral",
    show: bool = True,
    surface_alpha: float = 0.50,
    zoom: float = 1.0,
) -> "GridModel":
    """Build a grid model and optionally open the interactive explorer."""
    model = build_grid_model(grid, glevel, divide_face)
    if show:
        GridExplorer(
            grid=grid,
            glevel=glevel,
            divide_face=divide_face,
            surface_alpha=surface_alpha,
            zoom=zoom,
        ).show()
    return model


def ico_2d_decompose(
    glevel: int = 4,
    divide_face: int = 2,
    grid: str = "icosahedral",
    show: bool = True,
    surface_alpha: float = 0.50,
    zoom: float = 1.0,
) -> "GridModel":
    """Compatibility wrapper for the older script name."""
    return oasis_grid(
        glevel=glevel,
        divide_face=divide_face,
        grid=grid,
        show=show,
        surface_alpha=surface_alpha,
        zoom=zoom,
    )


def generate_oasis_icosahedral_points(
    glevel: int,
    nl_edge: int | None = None,
    return_coarse: bool = False,
):
    """Generate the OASIS ten-panel icosahedral grid points and panel indices."""
    if nl_edge is None:
        nl_edge = 2 ** glevel
    return generate_points(glevel, nl_edge, return_coarse=return_coarse)


# ---------------------------------------------------------------------------
# OASIS icosahedral point generation
# ---------------------------------------------------------------------------

def generate_points(glevel: int, nl_edge: int, return_coarse: bool = False):
    """
    Build all points on the unit sphere and the “coarse” rhombus indexing.

    This function uses a two-pass approach:
      - Pass 1: Run through the same loops, increment a local counter to find out
                exactly how many points are needed.
      - Allocate arrays of exactly that many rows.
      - Pass 2: Run through the loops again, this time filling in xyz (shape=(N,3))
                and rhomb (shape=(10, nl_edge**2)).

    Any time two antipodal points average to [0,0,0], we copy one endpoint
    instead of normalizing [0,0,0].

    return_coarse : bool
        When True, also return the full per-panel index grids with boundaries.
        Existing callers can keep the default False behavior.

    Returns
    -------
    xyz : np.ndarray, shape=(N, 3)
        All generated unit‐sphere points.
    rhomb : np.ndarray, shape=(10, nl_edge**2)
        For each of the 10 faces, the “flattened” grid‐point indices into xyz.
    rhombi_coarse : np.ndarray, optional, shape=(10, nl_edge + 1, nl_edge + 1)
        Full face grids including the final row/column of boundary indices.
    """
    # --- PASS 1: Count how many points the OASIS construction creates ---
    count_estimate = 12  # 2 poles + 10 equatorial “starter” vertices

    for _face in range(10):
        # (A) Five‐point blocks for grd=1..glevel
        for grd in range(1, glevel + 1):
            blocks = 2 ** (grd - 1)
            count_estimate += 5 * blocks

        # (B) Triangular fills for grd=2..glevel
        for grd in range(2, glevel + 1):
            rows = 2 ** (grd - 1) - 1
            tri_points_each = 3 * (rows * (rows + 1) // 2)
            count_estimate += 2 * tri_points_each  # upper + lower

    total_points = count_estimate  # exact, no buffer needed

    # Allocate exactly the right size:
    xyz = np.zeros((total_points, 3), dtype=float)

    # We will build rhombi_coarse[face, i, j] first (shape=(10, nl_edge+1, nl_edge+1)),
    # then flatten into rhomb[face, nl_edge**2]. rhomb gets its final shape at the end.
    rhombi_coarse = np.zeros((10, nl_edge + 1, nl_edge + 1), dtype=int)

    # Now we know exactly how many xyz‐rows we need, so we can fill them in pass 2:
    idx = 0

    # 1) Two poles:
    xyz[idx, :] = [0.0, 0.0, 1.0]   # North
    idx += 1
    xyz[idx, :] = [0.0, 0.0, -1.0]  # South
    idx += 1

    # 2) Ten “equatorial” vertices:
    w = 2.0 * np.arccos(1.0 / (2.0 * np.sin(np.pi / 5.0)))

    # (a) Five “top” vertices (just below north pole):
    phi_vals_top = [
        -np.pi / 5.0,
        np.pi / 5.0,
        3.0 * np.pi / 5.0,
        np.pi,
        -3.0 * np.pi / 5.0
    ]
    for phi in phi_vals_top:
        x = np.cos(phi) * np.cos(np.pi / 2.0 - w)
        y = np.sin(phi) * np.cos(np.pi / 2.0 - w)
        z = np.sin(np.pi / 2.0 - w)
        xyz[idx, :] = [x, y, z]
        face_idx = (idx - 2) % 10
        rhombi_coarse[face_idx, 1, 1] = idx
        idx += 1

    # (b) Five “bottom” vertices (just above south pole):
    phi_vals_bot = [
        0.0,
        2.0 * np.pi / 5.0,
        4.0 * np.pi / 5.0,
        -4.0 * np.pi / 5.0,
        -2.0 * np.pi / 5.0
    ]
    for phi in phi_vals_bot:
        x = np.cos(phi) * np.cos(w - np.pi / 2.0)
        y = np.sin(phi) * np.cos(w - np.pi / 2.0)
        z = np.sin(w - np.pi / 2.0)
        xyz[idx, :] = [x, y, z]
        face_idx = (idx - 2) % 10
        rhombi_coarse[face_idx, 1, 1] = idx
        idx += 1

    # 3) Hard‐wire each face’s four corner indices (0‐based):
    # Face 1 (1-based [3,1,8,4] -> 0-based [2,0,7,3]):
    rhombi_coarse[0, 0, 0]       = 2
    rhombi_coarse[0, 0, nl_edge] = 0
    rhombi_coarse[0, nl_edge, 0] = 7
    rhombi_coarse[0, nl_edge, nl_edge] = 3

    # Face 2 ([4,1,9,5] → [3,0,8,4]):
    rhombi_coarse[1, 0, 0]       = 3
    rhombi_coarse[1, 0, nl_edge] = 0
    rhombi_coarse[1, nl_edge, 0] = 8
    rhombi_coarse[1, nl_edge, nl_edge] = 4

    # Face 3 ([5,1,10,6] → [4,0,9,5]):
    rhombi_coarse[2, 0, 0]       = 4
    rhombi_coarse[2, 0, nl_edge] = 0
    rhombi_coarse[2, nl_edge, 0] = 9
    rhombi_coarse[2, nl_edge, nl_edge] = 5

    # Face 4 ([6,1,11,7] → [5,0,10,6]):
    rhombi_coarse[3, 0, 0]       = 5
    rhombi_coarse[3, 0, nl_edge] = 0
    rhombi_coarse[3, nl_edge, 0] = 10
    rhombi_coarse[3, nl_edge, nl_edge] = 6

    # Face 5 ([7,1,12,3] → [6,0,11,2]):
    rhombi_coarse[4, 0, 0]       = 6
    rhombi_coarse[4, 0, nl_edge] = 0
    rhombi_coarse[4, nl_edge, 0] = 11
    rhombi_coarse[4, nl_edge, nl_edge] = 2

    # Face 6 ([8,4,2,9] → [7,3,1,8]):
    rhombi_coarse[5, 0, 0]       = 7
    rhombi_coarse[5, 0, nl_edge] = 3
    rhombi_coarse[5, nl_edge, 0] = 1
    rhombi_coarse[5, nl_edge, nl_edge] = 8

    # Face 7 ([9,5,2,10] → [8,4,1,9]):
    rhombi_coarse[6, 0, 0]       = 8
    rhombi_coarse[6, 0, nl_edge] = 4
    rhombi_coarse[6, nl_edge, 0] = 1
    rhombi_coarse[6, nl_edge, nl_edge] = 9

    # Face 8 ([10,6,2,11] → [9,5,1,10]):
    rhombi_coarse[7, 0, 0]       = 9
    rhombi_coarse[7, 0, nl_edge] = 5
    rhombi_coarse[7, nl_edge, 0] = 1
    rhombi_coarse[7, nl_edge, nl_edge] = 10

    # Face 9 ([11,7,2,12] → [10,6,1,11]):
    rhombi_coarse[8, 0, 0]       = 10
    rhombi_coarse[8, 0, nl_edge] = 6
    rhombi_coarse[8, nl_edge, 0] = 1
    rhombi_coarse[8, nl_edge, nl_edge] = 11

    # Face 10 ([12,3,2,8] → [11,2,1,7]):
    rhombi_coarse[9, 0, 0]       = 11
    rhombi_coarse[9, 0, nl_edge] = 2
    rhombi_coarse[9, nl_edge, 0] = 1
    rhombi_coarse[9, nl_edge, nl_edge] = 7

    # 4) Subdivide each of the 10 faces in two stages:
    #    (A) “Five‐point block” loop (grd=1..glevel)
    #    (B) “Upper & lower triangles” loop (grd=2..glevel)
    for face in range(10):
        # (A) Five‐point block
        for grd in range(1, glevel + 1):
            ind_jump = 2 ** (glevel - grd)
            ind = ind_jump
            iterations = 2 ** (grd - 1)
            for _ in range(iterations):
                # 1) midpoint on top edge: (0, ind)
                p1 = rhombi_coarse[face, 0, ind + ind_jump]
                p2 = rhombi_coarse[face, 0, ind - ind_jump]
                avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                norm = np.linalg.norm(avg)
                if norm < 1e-12:
                    xyz[idx, :] = xyz[p1, :]
                else:
                    xyz[idx, :] = avg / norm
                rhombi_coarse[face, 0, ind] = idx
                idx += 1

                # 2) midpoint on left edge: (ind, 0)
                p1 = rhombi_coarse[face, ind + ind_jump, 0]
                p2 = rhombi_coarse[face, ind - ind_jump, 0]
                avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                norm = np.linalg.norm(avg)
                if norm < 1e-12:
                    xyz[idx, :] = xyz[p1, :]
                else:
                    xyz[idx, :] = avg / norm
                rhombi_coarse[face, ind, 0] = idx
                idx += 1

                # 3) midpoint on center diagonal: (ind, ind)
                p1 = rhombi_coarse[face, ind + ind_jump, ind + ind_jump]
                p2 = rhombi_coarse[face, ind - ind_jump, ind - ind_jump]
                avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                norm = np.linalg.norm(avg)
                if norm < 1e-12:
                    xyz[idx, :] = xyz[p1, :]
                else:
                    xyz[idx, :] = avg / norm
                rhombi_coarse[face, ind, ind] = idx
                idx += 1

                # 4) midpoint on bottom edge: (nl_edge, ind)
                p1 = rhombi_coarse[face, nl_edge, ind + ind_jump]
                p2 = rhombi_coarse[face, nl_edge, ind - ind_jump]
                avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                norm = np.linalg.norm(avg)
                if norm < 1e-12:
                    xyz[idx, :] = xyz[p1, :]
                else:
                    xyz[idx, :] = avg / norm
                rhombi_coarse[face, nl_edge, ind] = idx
                idx += 1

                # 5) midpoint on right edge: (ind, nl_edge)
                p1 = rhombi_coarse[face, ind + ind_jump, nl_edge]
                p2 = rhombi_coarse[face, ind - ind_jump, nl_edge]
                avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                norm = np.linalg.norm(avg)
                if norm < 1e-12:
                    xyz[idx, :] = xyz[p1, :]
                else:
                    xyz[idx, :] = avg / norm
                rhombi_coarse[face, ind, nl_edge] = idx
                idx += 1

                ind += 2 * ind_jump

        # (B) Upper & lower triangle subdivisions, only for grd = 2..glevel
        for grd in range(2, glevel + 1):
            ind_jump = 2 ** (glevel - grd)

            # Upper triangle
            indj = ind_jump
            max_j = 2 ** (grd - 1) - 1
            for j in range(1, max_j + 1):
                indx = ind_jump
                indy = 0
                for _ in range(j):
                    rindx = indx
                    rindy = nl_edge - indj + indy

                    # (a) midpoint of (rindx±ind_jump, rindy±ind_jump)
                    p1 = rhombi_coarse[face, rindx - ind_jump, rindy - ind_jump]
                    p2 = rhombi_coarse[face, rindx + ind_jump, rindy + ind_jump]
                    avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                    norm = np.linalg.norm(avg)
                    if norm < 1e-12:
                        xyz[idx, :] = xyz[p1, :]
                    else:
                        xyz[idx, :] = avg / norm
                    rhombi_coarse[face, rindx, rindy] = idx
                    idx += 1

                    # (b) midpoint of (rindx±ind_jump, rindy-ind_jump)
                    p1 = rhombi_coarse[face, rindx - ind_jump, rindy - ind_jump]
                    p2 = rhombi_coarse[face, rindx + ind_jump, rindy - ind_jump]
                    avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                    norm = np.linalg.norm(avg)
                    if norm < 1e-12:
                        xyz[idx, :] = xyz[p1, :]
                    else:
                        xyz[idx, :] = avg / norm
                    rhombi_coarse[face, rindx, rindy - ind_jump] = idx
                    idx += 1

                    # (c) midpoint of (rindx + ind_jump, rindy±ind_jump)
                    p1 = rhombi_coarse[face, rindx + ind_jump, rindy + ind_jump]
                    p2 = rhombi_coarse[face, rindx + ind_jump, rindy - ind_jump]
                    avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                    norm = np.linalg.norm(avg)
                    if norm < 1e-12:
                        xyz[idx, :] = xyz[p1, :]
                    else:
                        xyz[idx, :] = avg / norm
                    rhombi_coarse[face, rindx + ind_jump, rindy] = idx
                    idx += 1

                    indx += 2 * ind_jump
                    indy += 2 * ind_jump
                indj += 2 * ind_jump

            # Lower triangle
            indj = ind_jump
            for j in range(1, max_j + 1):
                indx = 0
                indy = ind_jump
                for _ in range(j):
                    rindx = nl_edge - indj + indx
                    rindy = indy

                    # (a) midpoint of (rindx±ind_jump, rindy±ind_jump)
                    p1 = rhombi_coarse[face, rindx - ind_jump, rindy - ind_jump]
                    p2 = rhombi_coarse[face, rindx + ind_jump, rindy + ind_jump]
                    avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                    norm = np.linalg.norm(avg)
                    if norm < 1e-12:
                        xyz[idx, :] = xyz[p1, :]
                    else:
                        xyz[idx, :] = avg / norm
                    rhombi_coarse[face, rindx, rindy] = idx
                    idx += 1

                    # (b) midpoint of (rindx - ind_jump, rindy±ind_jump)
                    p1 = rhombi_coarse[face, rindx - ind_jump, rindy - ind_jump]
                    p2 = rhombi_coarse[face, rindx - ind_jump, rindy + ind_jump]
                    avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                    norm = np.linalg.norm(avg)
                    if norm < 1e-12:
                        xyz[idx, :] = xyz[p1, :]
                    else:
                        xyz[idx, :] = avg / norm
                    rhombi_coarse[face, rindx - ind_jump, rindy] = idx
                    idx += 1

                    # (c) midpoint of (rindx±ind_jump, rindy + ind_jump)
                    p1 = rhombi_coarse[face, rindx + ind_jump, rindy + ind_jump]
                    p2 = rhombi_coarse[face, rindx - ind_jump, rindy + ind_jump]
                    avg = 0.5 * (xyz[p1, :] + xyz[p2, :])
                    norm = np.linalg.norm(avg)
                    if norm < 1e-12:
                        xyz[idx, :] = xyz[p1, :]
                    else:
                        xyz[idx, :] = avg / norm
                    rhombi_coarse[face, rindx, rindy + ind_jump] = idx
                    idx += 1

                    indx += 2 * ind_jump
                    indy += 2 * ind_jump
                indj += 2 * ind_jump

    # At this point, idx == total_points exactly.

    # 5) Flatten rhombi_coarse → rhomb (shape = (10, nl_edge**2))
    rhomb = np.zeros((10, nl_edge * nl_edge), dtype=int)
    for face in range(10):
        for i in range(nl_edge):
            for j in range(nl_edge):
                rhomb[face, j * nl_edge + i] = rhombi_coarse[face, i, j]

    if return_coarse:
        return xyz, rhomb, rhombi_coarse

    return xyz, rhomb


def generate_rhombi(rhomb: np.ndarray, divide_face: int, n_region: int):
    """
    Subdivide each coarse rhombus into (divide_face × divide_face) smaller blocks.

    Parameters
    ----------
    rhomb : np.ndarray, shape=(10, nl_coarse**2)
        “Coarse” face→global index mapping (nl_coarse = sqrt(n_region)).
    divide_face : int
        How many times to subdivide each face (0 < divide_face < glevel).
    n_region : int
        Number of coarse subregions per face = (2**glevel)**2.

    Returns
    -------
    faces_rhombi : np.ndarray, shape=(10*4**divide_face, nl_sub**2)
        For each of the 10 faces and each small block, a flattened list
        of global indices (length = nl_sub**2).
    nl_sub : int
        Number of points along each small sub‐block edge:
        nl_sub = sqrt(n_region) // 2**divide_face.
    """
    faces_coarse = 10
    nl_coarse = int(np.sqrt(n_region))
    nl_sub = nl_coarse // (2 ** divide_face)
    infaces = 4 ** divide_face
    total_small = faces_coarse * infaces

    faces_rhombi = np.zeros((total_small, nl_sub * nl_sub), dtype=int)
    sqrt_in = int(np.sqrt(infaces))

    for f in range(faces_coarse):
        for kx in range(sqrt_in):
            for ky in range(sqrt_in):
                small_idx = f * infaces + ky * sqrt_in + kx
                for i in range(nl_sub):
                    for j in range(nl_sub):
                        coarse_row = ky * nl_sub + j
                        coarse_col = kx * nl_sub + i
                        big_index = coarse_row * nl_coarse + coarse_col
                        faces_rhombi[small_idx, j * nl_sub + i] = rhomb[f, big_index]

    return faces_rhombi, nl_sub


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GridDomain:
    """A contiguous MPI-style subdomain cut from one base panel.

    The default rank order is panel-major, then ``tile_y``, then ``tile_x``.
    ``point_grid`` stores global point indices with shape
    ``(cells_per_edge + 1, cells_per_edge + 1)``.
    """

    index: int
    panel_index: int
    panel_name: str
    tile_x: int
    tile_y: int
    point_grid: np.ndarray

    @property
    def cells_per_edge(self) -> int:
        return int(self.point_grid.shape[0] - 1)

    @property
    def rank(self) -> int:
        """MPI-style rank for this domain in the default decomposition."""
        return self.index

    @property
    def cell_count(self) -> int:
        return self.cells_per_edge * self.cells_per_edge

    @property
    def point_count(self) -> int:
        return int(self.point_grid.size)

    @property
    def shape(self) -> tuple[int, int]:
        return tuple(int(value) for value in self.point_grid.shape)


@dataclass(frozen=True)
class GridModel:
    """Complete spherical grid geometry and decomposition metadata.

    ``points`` stores normalized xyz coordinates. ``panel_grids`` maps each
    panel-local logical grid coordinate to a global point index. ``domains`` is
    derived from those panel grids and represents the MPI-style decomposition.
    """

    kind: str
    title: str
    points: np.ndarray
    panel_grids: list[np.ndarray]
    panel_names: list[str]
    domains: list[GridDomain]
    glevel: int
    divide_face: int
    projection: str = ""

    @property
    def edge_divisions(self) -> int:
        return 2 ** self.glevel

    @property
    def domains_per_panel_edge(self) -> int:
        return 2 ** self.divide_face

    @property
    def total_cells(self) -> int:
        return len(self.panel_grids) * self.edge_divisions * self.edge_divisions

    @property
    def panel_count(self) -> int:
        return len(self.panel_grids)

    @property
    def domain_count(self) -> int:
        return len(self.domains)


DOMAIN_EDGE_NAMES = ("south", "north", "west", "east")


# ---------------------------------------------------------------------------
# Postprocessing helpers for MPI-domain analysis
# ---------------------------------------------------------------------------

def postprocess_grid(model: GridModel) -> dict:
    """Return a compact summary of the grid and MPI-style decomposition."""
    domain_edge = 2 ** (model.glevel - model.divide_face)
    radii = np.linalg.norm(model.points, axis=1)
    summary = {
        "grid": model.kind,
        "title": model.title,
        "glevel": int(model.glevel),
        "divide_face": int(model.divide_face),
        "edge_divisions": int(model.edge_divisions),
        "panel_count": int(model.panel_count),
        "domains_per_panel_edge": int(model.domains_per_panel_edge),
        "domain_count": int(model.domain_count),
        "domain_cells_per_edge": int(domain_edge),
        "domain_cell_count": int(domain_edge * domain_edge),
        "total_cells": int(model.total_cells),
        "point_count": int(len(model.points)),
        "unit_sphere_max_error": float(np.max(np.abs(radii - 1.0))),
        "mpi_rank_order": "panel-major, tile_y, tile_x",
    }
    if model.projection:
        summary["projection"] = model.projection
    return summary


def build_mpi_domain_map(model: GridModel) -> np.ndarray:
    """Return panel/tile_y/tile_x -> MPI rank domain map."""
    tiles = model.domains_per_panel_edge
    rank_map = np.full((model.panel_count, tiles, tiles), -1, dtype=int)
    for domain in model.domains:
        rank_map[domain.panel_index, domain.tile_y, domain.tile_x] = domain.rank
    return rank_map


def postprocess_mpi_domains(model: GridModel) -> dict:
    """Return summary, rank map, edge adjacency, and per-domain metadata."""
    edge_adjacency = build_domain_edge_adjacency(model)
    adjacency = _domain_adjacency_from_edges(model, edge_adjacency)
    domains = [
        _with_domain_neighbors(
            postprocess_domain(model, domain.rank, include_neighbors=False),
            edge_adjacency,
            adjacency,
        )
        for domain in model.domains
    ]
    return {
        "summary": postprocess_grid(model),
        "rank_map": build_mpi_domain_map(model),
        "edge_adjacency": edge_adjacency,
        "adjacency": adjacency,
        "domains": domains,
    }


def postprocess_domains(
    model: GridModel,
    domain_indices: Sequence[int] | None = None,
    include_neighbors: bool = True,
) -> list[dict]:
    """Return metadata rows for all or selected MPI domains."""
    if domain_indices is None:
        domain_indices = range(model.domain_count)
    rows = [
        postprocess_domain(model, int(index), include_neighbors=False)
        for index in domain_indices
    ]
    if not include_neighbors:
        return rows

    edge_adjacency = build_domain_edge_adjacency(model)
    adjacency = _domain_adjacency_from_edges(model, edge_adjacency)
    return [_with_domain_neighbors(row, edge_adjacency, adjacency) for row in rows]


def postprocess_domain(
    model: GridModel,
    domain_index: int,
    include_neighbors: bool = True,
) -> dict:
    """Return a serializable summary for one MPI domain/rank."""
    domain = _get_domain(model, domain_index)
    center = _normalised_center(model.points[domain.point_grid.ravel()])
    lon, lat = _xyz_to_lon_lat(center)
    row = {
        "rank": int(domain.rank),
        "domain_index": int(domain.index),
        "panel_index": int(domain.panel_index),
        "panel_name": domain.panel_name,
        "tile_x": int(domain.tile_x),
        "tile_y": int(domain.tile_y),
        "point_grid_shape": [int(value) for value in domain.point_grid.shape],
        "cells_per_edge": int(domain.cells_per_edge),
        "cell_count": int(domain.cell_count),
        "point_count": int(domain.point_count),
        "unique_point_count": int(len(domain_point_indices(domain))),
        "boundary_point_count": int(len(domain_boundary_point_indices(domain))),
        "interior_point_count": int(len(domain_interior_point_indices(domain))),
        "center_xyz": [float(value) for value in center],
        "center_lon_lat_deg": [float(lon), float(lat)],
        "point_index_min": int(np.min(domain.point_grid)),
        "point_index_max": int(np.max(domain.point_grid)),
    }
    if include_neighbors:
        row["edge_neighbors"] = build_domain_edge_adjacency(model)[domain.rank]
        row["neighbors"] = build_domain_adjacency(model)[domain.rank]
    return row


def extract_domain_mesh(model: GridModel, domain_index: int) -> dict:
    """Return arrays needed to inspect or postprocess one domain mesh.

    The returned dictionary contains both point-index arrays and xyz arrays.
    This is the most useful function when comparing the explorer output against
    OASIS 3D-model fields stored per MPI rank.
    """
    domain = _get_domain(model, domain_index)
    cell_points = domain_cell_point_indices(domain)
    point_ids = domain_point_indices(domain)
    boundary_ids = domain_boundary_point_indices(domain)
    return {
        "rank": int(domain.rank),
        "panel_index": int(domain.panel_index),
        "panel_name": domain.panel_name,
        "tile_x": int(domain.tile_x),
        "tile_y": int(domain.tile_y),
        "point_grid": domain.point_grid.copy(),
        "point_indices": point_ids,
        "boundary_point_indices": boundary_ids,
        "interior_point_indices": domain_interior_point_indices(domain),
        "point_xyz": model.points[point_ids],
        "boundary_xyz": model.points[boundary_ids],
        "cell_point_indices": cell_points,
        "cell_centers_xyz": domain_cell_centers(model, domain),
    }


def domain_point_indices(domain: GridDomain) -> np.ndarray:
    """Return sorted unique point indices used by a domain."""
    return np.unique(domain.point_grid.ravel())


def domain_boundary_point_indices(domain: GridDomain) -> np.ndarray:
    """Return sorted unique point indices on a domain boundary."""
    grid = domain.point_grid
    boundary = np.concatenate([grid[0, :], grid[-1, :], grid[:, 0], grid[:, -1]])
    return np.unique(boundary)


def domain_interior_point_indices(domain: GridDomain) -> np.ndarray:
    """Return sorted unique point indices strictly inside a domain."""
    if min(domain.point_grid.shape) <= 2:
        return np.asarray([], dtype=int)
    return np.unique(domain.point_grid[1:-1, 1:-1].ravel())


def domain_edge_point_indices(domain: GridDomain) -> dict[str, np.ndarray]:
    """Return point-index arrays for each local domain edge."""
    grid = domain.point_grid
    return {
        "south": grid[0, :].copy(),
        "north": grid[-1, :].copy(),
        "west": grid[:, 0].copy(),
        "east": grid[:, -1].copy(),
    }


def domain_cell_point_indices(domain: GridDomain) -> np.ndarray:
    """Return cell corner point ids with shape (ny, nx, 4)."""
    grid = domain.point_grid
    south_west = grid[:-1, :-1]
    south_east = grid[:-1, 1:]
    north_east = grid[1:, 1:]
    north_west = grid[1:, :-1]
    return np.stack([south_west, south_east, north_east, north_west], axis=-1)


def domain_cell_centers(model: GridModel, domain: GridDomain | int) -> np.ndarray:
    """Return normalized xyz centers for each domain cell."""
    if isinstance(domain, int):
        domain = _get_domain(model, domain)
    cells = domain_cell_point_indices(domain)
    centers = np.mean(model.points[cells], axis=-2)
    norms = np.linalg.norm(centers, axis=-1, keepdims=True)
    return centers / np.maximum(norms, 1e-12)


def build_domain_adjacency(model: GridModel) -> dict[int, list[int]]:
    """Return MPI rank -> neighboring ranks based on shared domain edges."""
    edge_adjacency = build_domain_edge_adjacency(model)
    return _domain_adjacency_from_edges(model, edge_adjacency)


def _domain_adjacency_from_edges(
    model: GridModel,
    edge_adjacency: dict[int, dict[str, int | None]],
) -> dict[int, list[int]]:
    adjacency: dict[int, set[int]] = {domain.rank: set() for domain in model.domains}
    for rank, edges in edge_adjacency.items():
        for neighbor in edges.values():
            if neighbor is not None:
                adjacency[rank].add(int(neighbor))
    return {rank: sorted(neighbors) for rank, neighbors in adjacency.items()}


def _with_domain_neighbors(
    row: dict,
    edge_adjacency: dict[int, dict[str, int | None]],
    adjacency: dict[int, list[int]],
) -> dict:
    rank = int(row["rank"])
    row["edge_neighbors"] = edge_adjacency[rank]
    row["neighbors"] = adjacency[rank]
    return row


def build_domain_edge_adjacency(model: GridModel) -> dict[int, dict[str, int | None]]:
    """Return MPI rank -> local edge -> neighboring rank, if present.

    The icosahedral construction may hold separate global point ids on
    opposite sides of a panel edge. For that reason adjacency is keyed by
    rounded xyz coordinates rather than by raw point indices.
    """
    edge_lookup: dict[tuple[tuple[float, float, float], ...], list[tuple[int, str]]] = {}
    for domain in model.domains:
        for edge_name, point_ids in domain_edge_point_indices(domain).items():
            key = _canonical_edge_key(model.points[point_ids])
            edge_lookup.setdefault(key, []).append((domain.rank, edge_name))

    adjacency: dict[int, dict[str, int | None]] = {
        domain.rank: {edge_name: None for edge_name in DOMAIN_EDGE_NAMES}
        for domain in model.domains
    }
    for owners in edge_lookup.values():
        if len(owners) < 2:
            continue
        for rank, edge_name in owners:
            neighbor_ranks = sorted({other_rank for other_rank, _ in owners if other_rank != rank})
            if neighbor_ranks:
                adjacency[rank][edge_name] = int(neighbor_ranks[0])
    return adjacency


def export_postprocessing_npz(model: GridModel, path: str | Path) -> Path:
    """Save grid geometry and MPI-domain postprocessing arrays to an NPZ file."""
    output_path = Path(path)
    edge_adjacency = build_domain_edge_adjacency(model)
    edge_array = np.full((model.domain_count, len(DOMAIN_EDGE_NAMES)), -1, dtype=int)
    for domain in model.domains:
        for edge_index, edge_name in enumerate(DOMAIN_EDGE_NAMES):
            neighbor = edge_adjacency[domain.rank][edge_name]
            if neighbor is not None:
                edge_array[domain.rank, edge_index] = neighbor

    domain_table = np.asarray(
        [
            [
                domain.rank,
                domain.panel_index,
                domain.tile_x,
                domain.tile_y,
                domain.cells_per_edge,
                domain.cell_count,
                domain.point_count,
            ]
            for domain in model.domains
        ],
        dtype=int,
    )
    np.savez_compressed(
        output_path,
        grid_kind=np.asarray(model.kind),
        projection=np.asarray(model.projection),
        glevel=np.asarray(model.glevel),
        divide_face=np.asarray(model.divide_face),
        points=model.points,
        panel_grids=np.stack(model.panel_grids),
        domain_point_grids=np.stack([domain.point_grid for domain in model.domains]),
        domain_rank_map=build_mpi_domain_map(model),
        domain_table=domain_table,
        edge_adjacency=edge_array,
        edge_names=np.asarray(DOMAIN_EDGE_NAMES),
        domain_table_columns=np.asarray(
            [
                "rank",
                "panel_index",
                "tile_x",
                "tile_y",
                "cells_per_edge",
                "cell_count",
                "point_count",
            ]
        ),
    )
    return output_path


def _get_domain(model: GridModel, domain_index: int) -> GridDomain:
    if domain_index < 0 or domain_index >= model.domain_count:
        raise IndexError(
            f"domain index {domain_index} is outside 0..{model.domain_count - 1}"
        )
    return model.domains[domain_index]


def _point_key(point: np.ndarray, decimals: int = 10) -> tuple[float, float, float]:
    return tuple(float(value) for value in np.round(point, decimals))


def _canonical_edge_key(points: np.ndarray) -> tuple[tuple[float, float, float], ...]:
    keys = [_point_key(point) for point in points]
    return tuple(sorted(keys))


def _xyz_to_lon_lat(point: np.ndarray) -> tuple[float, float]:
    x, y, z = point
    lon = np.degrees(np.arctan2(y, x))
    lat = np.degrees(np.arcsin(np.clip(z, -1.0, 1.0)))
    return float(lon), float(lat)


# ---------------------------------------------------------------------------
# Grid model builders
# ---------------------------------------------------------------------------

def _validate_grid_parameters(glevel: int, divide_face: int) -> None:
    if glevel < 1:
        raise ValueError("glevel must be >= 1")
    if divide_face < 0 or divide_face > glevel:
        raise ValueError("divide_face must satisfy 0 <= divide_face <= glevel")


def _normalise_grid_name(grid: str) -> str:
    key = grid.strip().lower().replace("_", "-")
    if key in {"ico", "icosahedral", "icosphere", "icosahedral-rhombic"}:
        return "icosahedral"
    if key in {"cube", "cube-sphere", "cubed-sphere", "cubesphere"}:
        return "cube-sphere"
    raise ValueError(f"unknown grid type: {grid!r}")


def build_grid_model(
    grid: str = "icosahedral",
    glevel: int = 4,
    divide_face: int = 2,
    projection: str = "equiangular",
) -> GridModel:
    """Build either the OASIS icosahedral or cube-sphere grid model."""
    kind = _normalise_grid_name(grid)
    if kind == "icosahedral":
        return build_icosahedral_grid(glevel, divide_face)
    return build_cube_sphere_grid(glevel, divide_face, projection=projection)


def build_icosahedral_grid(glevel: int = 4, divide_face: int = 2) -> GridModel:
    """Create the OASIS ten-panel icosahedral rhombic grid.

    ``glevel`` controls panel resolution: each panel edge has ``2**glevel``
    cells. ``divide_face`` controls the MPI-style split: each panel contains
    ``4**divide_face`` domains.
    """
    _validate_grid_parameters(glevel, divide_face)
    nl_edge = 2 ** glevel
    xyz, _rhomb, rhombi_coarse = generate_oasis_icosahedral_points(
        glevel,
        nl_edge,
        return_coarse=True,
    )
    panel_grids = [rhombi_coarse[i, :, :].copy() for i in range(10)]
    panel_names = [f"Ico {i + 1:02d}" for i in range(10)]
    domains = _split_panel_domains(panel_grids, panel_names, divide_face)
    return GridModel(
        kind="icosahedral",
        title="Icosahedral rhombic grid",
        points=xyz,
        panel_grids=panel_grids,
        panel_names=panel_names,
        domains=domains,
        glevel=glevel,
        divide_face=divide_face,
    )


def build_cube_sphere_grid(
    glevel: int = 4,
    divide_face: int = 2,
    projection: str = "equiangular",
) -> GridModel:
    """Create a six-panel cube-sphere grid projected to the unit sphere.

    The panel/domain bookkeeping intentionally mirrors the icosahedral path so
    postprocessing functions can compare both grid families with the same API.
    ``projection='equiangular'`` spaces coordinates uniformly in cube-face
    angle; ``projection='normalized'`` spaces cube coordinates linearly before
    normalizing to the sphere.
    """
    _validate_grid_parameters(glevel, divide_face)
    projection = projection.strip().lower()
    if projection not in {"equiangular", "normalized"}:
        raise ValueError("projection must be either 'equiangular' or 'normalized'")

    nl_edge = 2 ** glevel
    coordinates = _cube_face_coordinates(nl_edge, projection)
    face_frames = [
        ("+X", np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0])),
        ("-X", np.array([-1.0, 0.0, 0.0]), np.array([0.0, -1.0, 0.0]), np.array([0.0, 0.0, 1.0])),
        ("+Y", np.array([0.0, 1.0, 0.0]), np.array([-1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])),
        ("-Y", np.array([0.0, -1.0, 0.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])),
        ("+Z", np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])),
        ("-Z", np.array([0.0, 0.0, -1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, -1.0, 0.0])),
    ]

    points: list[np.ndarray] = []
    point_lookup: dict[tuple[float, float, float], int] = {}
    panel_grids: list[np.ndarray] = []
    panel_names: list[str] = []

    for name, center, u_axis, v_axis in face_frames:
        grid = np.empty((nl_edge + 1, nl_edge + 1), dtype=int)
        for row, v_value in enumerate(coordinates):
            for col, u_value in enumerate(coordinates):
                vector = center + u_value * u_axis + v_value * v_axis
                point = vector / np.linalg.norm(vector)
                grid[row, col] = _deduplicated_point_index(
                    point,
                    points,
                    point_lookup,
                )
        panel_grids.append(grid)
        panel_names.append(f"Cube {name}")

    domains = _split_panel_domains(panel_grids, panel_names, divide_face)
    return GridModel(
        kind="cube-sphere",
        title="Cube-sphere grid",
        points=np.asarray(points, dtype=float),
        panel_grids=panel_grids,
        panel_names=panel_names,
        domains=domains,
        glevel=glevel,
        divide_face=divide_face,
        projection=projection,
    )


def _cube_face_coordinates(nl_edge: int, projection: str) -> np.ndarray:
    if projection == "normalized":
        return np.linspace(-1.0, 1.0, nl_edge + 1)
    angles = np.linspace(-np.pi / 4.0, np.pi / 4.0, nl_edge + 1)
    return np.tan(angles)


def _deduplicated_point_index(
    point: np.ndarray,
    points: list[np.ndarray],
    point_lookup: dict[tuple[float, float, float], int],
) -> int:
    key = tuple(float(value) for value in np.round(point, 12))
    if key in point_lookup:
        return point_lookup[key]
    point_lookup[key] = len(points)
    points.append(point)
    return point_lookup[key]


def _split_panel_domains(
    panel_grids: list[np.ndarray],
    panel_names: list[str],
    divide_face: int,
) -> list[GridDomain]:
    domains: list[GridDomain] = []
    tiles_per_edge = 2 ** divide_face

    for panel_index, panel_grid in enumerate(panel_grids):
        edge = int(panel_grid.shape[0] - 1)
        if edge % tiles_per_edge != 0:
            raise ValueError("panel edge divisions are not divisible by domains")
        tile_edge = edge // tiles_per_edge

        for tile_y in range(tiles_per_edge):
            for tile_x in range(tiles_per_edge):
                row0 = tile_y * tile_edge
                col0 = tile_x * tile_edge
                point_grid = panel_grid[
                    row0: row0 + tile_edge + 1,
                    col0: col0 + tile_edge + 1,
                ].copy()
                domains.append(
                    GridDomain(
                        index=len(domains),
                        panel_index=panel_index,
                        panel_name=panel_names[panel_index],
                        tile_x=tile_x,
                        tile_y=tile_y,
                        point_grid=point_grid,
                    )
                )

    return domains


# ---------------------------------------------------------------------------
# Interactive explorer
# ---------------------------------------------------------------------------

class GridExplorer:
    """Interactive Matplotlib explorer for both supported spherical grids.

    The explorer keeps the numerical grid model separate from rendering state.
    Changing the grid family, refinement level, or domain split rebuilds the
    model; changing opacity, zoom, labels, points, or linework only redraws the
    current model.
    """

    MAX_LEVEL = 6
    MIN_ZOOM = 0.55
    MAX_ZOOM = 2.75
    BASE_VIEW_LIMIT = 1.08
    GRID_LABELS = ("Icosahedral", "Cube-sphere")
    COLOR_LABELS = ("Domains", "Panels", "Height")

    def __init__(
        self,
        grid: str = "icosahedral",
        glevel: int = 4,
        divide_face: int = 2,
        projection: str = "equiangular",
        focus_domain: int | None = None,
        surface_alpha: float = 0.50,
        zoom: float = 1.0,
    ) -> None:
        self.grid_kind = _normalise_grid_name(grid)
        self.glevel = max(1, min(self.MAX_LEVEL, int(glevel)))
        self.divide_face = max(0, min(int(divide_face), self.glevel))
        self.projection = projection
        self.focus_domain = focus_domain

        self.show_surface = True
        self.show_lines = True
        self.show_points = False
        self.show_labels = True
        self.surface_alpha = float(np.clip(surface_alpha, 0.0, 1.0))
        self.zoom = float(np.clip(zoom, self.MIN_ZOOM, self.MAX_ZOOM))
        self.color_mode = "Domains"
        self._callbacks_suspended = False
        self._view_initialised = False

        self.model = build_grid_model(
            self.grid_kind,
            self.glevel,
            self.divide_face,
            projection=self.projection,
        )
        self._clamp_focus_domain()
        self._setup_figure()
        self._draw()

    def _setup_figure(self) -> None:
        self.fig = plt.figure(figsize=(14.0, 8.4), facecolor="#f5f7fb")
        self.ax = self.fig.add_axes([0.30, 0.06, 0.67, 0.88], projection="3d")

        self.fig.text(
            0.04,
            0.955,
            "OASIS grid explorer",
            fontsize=15,
            fontweight="bold",
            color="#111827",
        )
        self.fig.text(
            0.04,
            0.925,
            "Inspect panels, domains, grid lines, and focused subdomains.",
            fontsize=9,
            color="#4b5563",
        )
        self.info_text = self.fig.text(
            0.04,
            0.095,
            "",
            fontsize=7.6,
            color="#1f2937",
            linespacing=1.05,
            family="monospace",
            va="top",
        )

        grid_radio_ax = self.fig.add_axes([0.04, 0.77, 0.20, 0.11], facecolor="#eef2f7")
        active_grid = 0 if self.grid_kind == "icosahedral" else 1
        self.grid_radio = RadioButtons(grid_radio_ax, self.GRID_LABELS, active=active_grid)
        self.grid_radio.on_clicked(self._on_grid_change)
        self._style_radio(self.grid_radio)

        color_radio_ax = self.fig.add_axes([0.04, 0.61, 0.20, 0.12], facecolor="#eef2f7")
        self.color_radio = RadioButtons(color_radio_ax, self.COLOR_LABELS, active=0)
        self.color_radio.on_clicked(self._on_color_change)
        self._style_radio(self.color_radio)

        level_ax = self.fig.add_axes([0.065, 0.545, 0.17, 0.025], facecolor="#e5e7eb")
        self.level_slider = Slider(
            level_ax,
            "level",
            1,
            self.MAX_LEVEL,
            valinit=self.glevel,
            valstep=1,
            color="#2563eb",
        )
        self.level_slider.on_changed(self._on_level_change)

        divide_ax = self.fig.add_axes([0.065, 0.495, 0.17, 0.025], facecolor="#e5e7eb")
        self.divide_slider = Slider(
            divide_ax,
            "split",
            0,
            self.MAX_LEVEL,
            valinit=self.divide_face,
            valstep=1,
            color="#0f766e",
        )
        self.divide_slider.on_changed(self._on_divide_change)

        domain_ax = self.fig.add_axes([0.065, 0.445, 0.17, 0.025], facecolor="#e5e7eb")
        self.domain_slider = Slider(
            domain_ax,
            "domain",
            0,
            max(1, len(self.model.domains) - 1),
            valinit=0,
            valstep=1,
            color="#d97706",
        )
        self.domain_slider.on_changed(self._on_domain_change)

        alpha_ax = self.fig.add_axes([0.065, 0.405, 0.17, 0.025], facecolor="#e5e7eb")
        self.surface_alpha_slider = Slider(
            alpha_ax,
            "alpha",
            0.0,
            1.0,
            valinit=self.surface_alpha,
            valstep=0.05,
            color="#7c3aed",
        )
        self.surface_alpha_slider.on_changed(self._on_surface_alpha_change)

        zoom_ax = self.fig.add_axes([0.065, 0.365, 0.17, 0.025], facecolor="#e5e7eb")
        self.zoom_slider = Slider(
            zoom_ax,
            "zoom",
            self.MIN_ZOOM,
            self.MAX_ZOOM,
            valinit=self.zoom,
            valstep=0.05,
            color="#475569",
        )
        self.zoom_slider.on_changed(self._on_zoom_change)

        check_ax = self.fig.add_axes([0.04, 0.225, 0.20, 0.115], facecolor="#eef2f7")
        self.checks = CheckButtons(
            check_ax,
            ("surface", "grid lines", "points", "labels"),
            (self.show_surface, self.show_lines, self.show_points, self.show_labels),
        )
        self.checks.on_clicked(self._on_check_change)
        self._style_checks(self.checks)

        self.all_button = self._add_button([0.04, 0.175, 0.055, 0.035], "All", self._show_all)
        self.prev_button = self._add_button([0.105, 0.175, 0.055, 0.035], "Prev", self._previous_domain)
        self.next_button = self._add_button([0.17, 0.175, 0.055, 0.035], "Next", self._next_domain)
        self.save_button = self._add_button([0.04, 0.125, 0.185, 0.035], "Save view", self._save_from_button)

        self.fig.canvas.mpl_connect("key_press_event", self._on_key_press)
        self.fig.canvas.mpl_connect("button_release_event", self._on_view_release)
        self._sync_domain_slider_range()

    @staticmethod
    def _style_radio(radio: RadioButtons) -> None:
        for label in radio.labels:
            label.set_fontsize(9)
            label.set_color("#111827")
        radio.ax.set_xticks([])
        radio.ax.set_yticks([])
        for spine in radio.ax.spines.values():
            spine.set_edgecolor("#d1d5db")

    @staticmethod
    def _style_checks(checks: CheckButtons) -> None:
        for label in checks.labels:
            label.set_fontsize(9)
            label.set_color("#111827")
        checks.ax.set_xticks([])
        checks.ax.set_yticks([])
        for spine in checks.ax.spines.values():
            spine.set_edgecolor("#d1d5db")

    def _add_button(self, bounds: list[float], label: str, callback) -> Button:
        button_ax = self.fig.add_axes(bounds, facecolor="#e5e7eb")
        button = Button(button_ax, label, color="#e5e7eb", hovercolor="#dbeafe")
        button.label.set_fontsize(9)
        button.label.set_color("#111827")
        button.on_clicked(callback)
        return button

    def _on_grid_change(self, label: str) -> None:
        if self._callbacks_suspended:
            return
        self.grid_kind = _normalise_grid_name(label)
        self.focus_domain = None
        self._rebuild_model()

    def _on_color_change(self, label: str) -> None:
        if self._callbacks_suspended:
            return
        self.color_mode = label
        self._draw()

    def _on_level_change(self, _value: float) -> None:
        if self._callbacks_suspended:
            return
        self.glevel = int(self.level_slider.val)
        if self.divide_face > self.glevel:
            self.divide_face = self.glevel
            self._callbacks_suspended = True
            self.divide_slider.set_val(self.divide_face)
            self._callbacks_suspended = False
        self._rebuild_model()

    def _on_divide_change(self, _value: float) -> None:
        if self._callbacks_suspended:
            return
        self.divide_face = min(int(self.divide_slider.val), self.glevel)
        if int(self.divide_slider.val) != self.divide_face:
            self._callbacks_suspended = True
            self.divide_slider.set_val(self.divide_face)
            self._callbacks_suspended = False
        self.focus_domain = None
        self._rebuild_model()

    def _on_domain_change(self, value: float) -> None:
        if self._callbacks_suspended:
            return
        self.focus_domain = int(round(value))
        self._clamp_focus_domain()
        self._draw()

    def _on_surface_alpha_change(self, value: float) -> None:
        if self._callbacks_suspended:
            return
        self.surface_alpha = float(np.clip(value, 0.0, 1.0))
        self._draw()

    def _on_zoom_change(self, value: float) -> None:
        if self._callbacks_suspended:
            return
        self.zoom = float(np.clip(value, self.MIN_ZOOM, self.MAX_ZOOM))
        self._draw()

    def _on_check_change(self, _label: str) -> None:
        statuses = self.checks.get_status()
        self.show_surface = bool(statuses[0])
        self.show_lines = bool(statuses[1])
        self.show_points = bool(statuses[2])
        self.show_labels = bool(statuses[3])
        self._draw()

    def _on_key_press(self, event) -> None:
        if event.key in {"right", "n"}:
            self._next_domain(None)
        elif event.key in {"left", "p"}:
            self._previous_domain(None)
        elif event.key in {"a", "escape"}:
            self._show_all(None)
        elif event.key == "g":
            self.grid_kind = "cube-sphere" if self.grid_kind == "icosahedral" else "icosahedral"
            self.focus_domain = None
            self._callbacks_suspended = True
            self.grid_radio.set_active(0 if self.grid_kind == "icosahedral" else 1)
            self._callbacks_suspended = False
            self._rebuild_model()
        elif event.key in {"+", "="}:
            self._set_zoom(self.zoom + 0.10)
        elif event.key in {"-", "_"}:
            self._set_zoom(self.zoom - 0.10)

    def _on_view_release(self, event) -> None:
        """Refresh hidden-line culling after the user rotates the 3D view."""
        if event.inaxes is self.ax:
            self._draw()

    def _show_all(self, _event) -> None:
        self.focus_domain = None
        self._draw()

    def _previous_domain(self, _event) -> None:
        self._step_domain(-1)

    def _next_domain(self, _event) -> None:
        self._step_domain(1)

    def _step_domain(self, delta: int) -> None:
        domain_count = len(self.model.domains)
        current = 0 if self.focus_domain is None else self.focus_domain
        self.focus_domain = (current + delta) % domain_count
        self._callbacks_suspended = True
        self.domain_slider.set_val(self.focus_domain)
        self._callbacks_suspended = False
        self._draw()

    def _save_from_button(self, _event) -> None:
        self.save()

    def _set_zoom(self, zoom: float) -> None:
        self.zoom = float(np.clip(zoom, self.MIN_ZOOM, self.MAX_ZOOM))
        self._callbacks_suspended = True
        self.zoom_slider.set_val(self.zoom)
        self._callbacks_suspended = False
        self._draw()

    def _rebuild_model(self) -> None:
        self.model = build_grid_model(
            self.grid_kind,
            self.glevel,
            self.divide_face,
            projection=self.projection,
        )
        self._clamp_focus_domain()
        self._sync_domain_slider_range()
        self._draw()

    def _clamp_focus_domain(self) -> None:
        if self.focus_domain is None:
            return
        self.focus_domain = max(0, min(int(self.focus_domain), len(self.model.domains) - 1))

    def _sync_domain_slider_range(self) -> None:
        max_domain = max(1, len(self.model.domains) - 1)
        self.domain_slider.valmax = max_domain
        self.domain_slider.ax.set_xlim(self.domain_slider.valmin, max_domain)
        slider_value = 0 if self.focus_domain is None else self.focus_domain
        slider_value = min(slider_value, max_domain)
        self._callbacks_suspended = True
        self.domain_slider.set_val(slider_value)
        self._callbacks_suspended = False

    def _draw(self) -> None:
        if self._view_initialised:
            elev = self.ax.elev
            azim = self.ax.azim
        else:
            elev = 23
            azim = -42
            self._view_initialised = True

        self.ax.cla()
        self.ax.set_facecolor("#0b1020")
        self.ax.set_axis_off()
        view_limit = self.BASE_VIEW_LIMIT / self.zoom
        self.ax.set_xlim(-view_limit, view_limit)
        self.ax.set_ylim(-view_limit, view_limit)
        self.ax.set_zlim(-view_limit, view_limit)
        self.ax.set_box_aspect([1, 1, 1])
        self.ax.view_init(elev=elev, azim=azim)

        self._draw_reference_sphere()
        if self.show_surface:
            self._draw_domain_surfaces()
        if self.show_lines:
            self._draw_grid_lines()
        if self.show_points:
            self._draw_points()
        if self.show_labels:
            self._draw_labels()

        self._update_info_text()
        self.fig.canvas.draw_idle()

    def _hide_backside_grid(self) -> bool:
        """Return True when opaque surfaces should occlude back-side overlays."""
        return self.show_surface and self.surface_alpha >= 0.999

    def _camera_direction(self) -> np.ndarray:
        """Unit vector from the sphere center toward the current camera."""
        elev = np.deg2rad(self.ax.elev)
        azim = np.deg2rad(self.ax.azim)
        return np.array(
            [
                np.cos(elev) * np.cos(azim),
                np.cos(elev) * np.sin(azim),
                np.sin(elev),
            ],
            dtype=float,
        )

    def _draw_reference_sphere(self) -> None:
        # The reference sphere is only a subtle backdrop; the domain surface is
        # the authoritative visible grid. Keeping this alpha capped avoids
        # washing out fully opaque domain colors.
        u = np.linspace(0.0, 2.0 * np.pi, 60)
        v = np.linspace(0.0, np.pi, 30)
        x = np.outer(np.cos(u), np.sin(v))
        y = np.outer(np.sin(u), np.sin(v))
        z = np.outer(np.ones_like(u), np.cos(v))
        self.ax.plot_surface(
            x,
            y,
            z,
            color="#94a3b8",
            alpha=min(0.18, self.surface_alpha * 0.11),
            linewidth=0,
            shade=False,
        )

    def _draw_domain_surfaces(self) -> None:
        domains = self.model.domains
        background_alpha = self.surface_alpha
        step = _surface_step(domains, max_quads=14000)
        quads: list[np.ndarray] = []
        colors: list[tuple[float, float, float, float]] = []

        for domain in domains:
            domain_quads = list(_domain_quads(domain, self.model.points, step=step))
            quads.extend(domain_quads)
            colors.extend(
                [_domain_color(domain, self.model, self.color_mode, background_alpha)]
                * len(domain_quads)
            )

        if quads:
            self.ax.add_collection3d(
                Poly3DCollection(
                    quads,
                    facecolors=colors,
                    edgecolors="none",
                    linewidths=0.0,
                    zsort="average",
                )
            )

        if self.focus_domain is not None:
            selected = self.model.domains[self.focus_domain]
            selected_step = 1 if selected.cells_per_edge <= 48 else 2
            selected_quads = list(
                _domain_quads(selected, self.model.points, step=selected_step)
            )
            if selected_quads:
                selected_alpha = self.surface_alpha
                self.ax.add_collection3d(
                    Poly3DCollection(
                        selected_quads,
                        facecolors=[(1.0, 0.78, 0.18, selected_alpha)] * len(selected_quads),
                        edgecolors=_rgba("#ffffff", self.surface_alpha),
                        linewidths=0.25,
                        zsort="average",
                    )
                )

    def _draw_grid_lines(self) -> None:
        grid_alpha = self.surface_alpha
        panel_lines = list(_panel_grid_lines(self.model.panel_grids, self.model.points))
        if self._hide_backside_grid():
            panel_lines = _visible_line_segments(panel_lines, self._camera_direction())
        if panel_lines:
            self.ax.add_collection3d(
                Line3DCollection(
                    panel_lines,
                    colors=[_rgba("#edf5ff", grid_alpha)] * len(panel_lines),
                    linewidths=0.45,
                )
            )

        boundary_lines = [
            line
            for domain in self.model.domains
            for line in _domain_boundary_lines(domain, self.model.points)
        ]
        if self._hide_backside_grid():
            boundary_lines = _visible_line_segments(boundary_lines, self._camera_direction())
        if boundary_lines:
            self.ax.add_collection3d(
                Line3DCollection(
                    boundary_lines,
                    colors=[_rgba("#ffffff", grid_alpha)] * len(boundary_lines),
                    linewidths=0.75,
                )
            )

        if self.focus_domain is not None:
            selected = self.model.domains[self.focus_domain]
            selected_lines = list(_domain_grid_lines(selected, self.model.points))
            selected_lines.extend(_domain_boundary_lines(selected, self.model.points))
            if self._hide_backside_grid():
                selected_lines = _visible_line_segments(
                    selected_lines,
                    self._camera_direction(),
                )
            if selected_lines:
                self.ax.add_collection3d(
                    Line3DCollection(
                        selected_lines,
                        colors=[_rgba("#ffee73", grid_alpha)] * len(selected_lines),
                        linewidths=1.25,
                    )
                )

    def _draw_points(self) -> None:
        points = self.model.points
        if self._hide_backside_grid():
            points = points[_front_points_mask(points, self._camera_direction())]
        self.ax.scatter(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            s=7,
            color=_rgba("#f8fafc", self.surface_alpha),
            edgecolors="none",
            depthshade=False,
        )
        if self.focus_domain is not None:
            domain = self.model.domains[self.focus_domain]
            selected_ids = np.unique(domain.point_grid.ravel())
            selected_points = self.model.points[selected_ids]
            if self._hide_backside_grid():
                selected_points = selected_points[
                    _front_points_mask(selected_points, self._camera_direction())
                ]
            self.ax.scatter(
                selected_points[:, 0],
                selected_points[:, 1],
                selected_points[:, 2],
                s=16,
                color=_rgba("#facc15", self.surface_alpha),
                edgecolors="none",
                depthshade=False,
            )

    def _draw_labels(self) -> None:
        camera_direction = self._camera_direction()
        for panel_grid, panel_name in zip(self.model.panel_grids, self.model.panel_names):
            center = _normalised_center(self.model.points[panel_grid.ravel()])
            if self._hide_backside_grid() and np.dot(center, camera_direction) < 0.0:
                continue
            self.ax.text(
                *(1.12 * center),
                panel_name,
                color=_rgba("#f8fafc", self.surface_alpha),
                fontsize=8,
                ha="center",
                va="center",
            )

        if self.focus_domain is not None:
            domain = self.model.domains[self.focus_domain]
            center = _normalised_center(self.model.points[domain.point_grid.ravel()])
            if self._hide_backside_grid() and np.dot(center, camera_direction) < 0.0:
                return
            self.ax.text(
                *(1.22 * center),
                f"Domain {domain.index}\n{domain.panel_name}",
                color=_rgba("#fde68a", self.surface_alpha),
                fontsize=9,
                fontweight="bold",
                ha="center",
                va="center",
            )

    def _update_info_text(self) -> None:
        domain_edge = 2 ** (self.glevel - self.divide_face)
        focus = "all"
        if self.focus_domain is not None:
            domain = self.model.domains[self.focus_domain]
            focus = (
                f"{domain.index} "
                f"({domain.panel_name}, x={domain.tile_x}, y={domain.tile_y})"
            )

        projection = ""
        if self.model.kind == "cube-sphere":
            projection = f" | projection {self.model.projection}"

        self.info_text.set_text(
            f"{self.model.title}\n"
            f"{self.model.kind} | level {self.glevel} | split {self.divide_face}"
            f"{projection}\n"
            f"panels {len(self.model.panel_grids)} | domains {len(self.model.domains)}"
            f" | domain edge {domain_edge}\n"
            f"cells {self.model.total_cells} | points {len(self.model.points)}"
            f" | alpha {self.surface_alpha:.2f} | zoom {self.zoom:.2f}\n"
            f"focus {focus}\n"
            "keys arrows/n/p focus | +/- zoom | a all | g grid"
        )

    def save(self, path: str | Path | None = None) -> Path:
        if path is None:
            path = (
                f"{self.model.kind}_level{self.glevel}"
                f"_split{self.divide_face}.png"
            )
        output_path = Path(path)
        self.fig.savefig(output_path, dpi=180, facecolor=self.fig.get_facecolor())
        print(f"Saved {output_path}")
        return output_path

    def show(self) -> None:
        plt.show()


# ---------------------------------------------------------------------------
# Rendering and geometry helper functions
# ---------------------------------------------------------------------------

def _normalised_center(points: np.ndarray) -> np.ndarray:
    center = np.mean(points, axis=0)
    norm = np.linalg.norm(center)
    if norm < 1e-12:
        return center
    return center / norm


def _rgba(color: str, alpha: float) -> tuple[float, float, float, float]:
    """Return a Matplotlib RGBA color using the shared grid opacity value."""
    return to_rgba(color, alpha=float(np.clip(alpha, 0.0, 1.0)))


def _front_points_mask(
    points: np.ndarray,
    camera_direction: np.ndarray,
    cutoff: float = 0.0,
) -> np.ndarray:
    """Return True for points on the camera-facing side of the sphere."""
    if points.size == 0:
        return np.asarray([], dtype=bool)
    return np.dot(points, camera_direction) >= cutoff


def _visible_line_segments(
    lines: Iterable[np.ndarray],
    camera_direction: np.ndarray,
    cutoff: float = 0.0,
) -> list[np.ndarray]:
    """Split polylines into front-facing segments for hidden-line rendering.

    Matplotlib's 3D artists do not perform true depth-buffer occlusion between
    opaque polygon collections and line/scatter collections. When the surface
    is opaque, this function removes segments whose midpoint lies on the back
    hemisphere, so back-side grid lines do not show through the sphere.
    """
    visible_segments: list[np.ndarray] = []
    for line in lines:
        line = np.asarray(line, dtype=float)
        if line.shape[0] < 2:
            continue
        starts = line[:-1]
        ends = line[1:]
        midpoints = 0.5 * (starts + ends)
        mask = _front_points_mask(midpoints, camera_direction, cutoff=cutoff)
        for start, end, is_visible in zip(starts, ends, mask):
            if is_visible:
                visible_segments.append(np.vstack([start, end]))
    return visible_segments


def _domain_color(
    domain: GridDomain,
    model: GridModel,
    color_mode: str,
    alpha: float,
) -> tuple[float, float, float, float]:
    if color_mode == "Panels":
        cmap = plt.get_cmap("tab10")
        rgba = cmap(domain.panel_index % 10)
    elif color_mode == "Height":
        cmap = plt.get_cmap("viridis")
        center = _normalised_center(model.points[domain.point_grid.ravel()])
        rgba = cmap((center[2] + 1.0) * 0.5)
    else:
        cmap = plt.get_cmap("tab20")
        rgba = cmap(domain.index % 20)
    return (float(rgba[0]), float(rgba[1]), float(rgba[2]), alpha)


def _surface_step(domains: Iterable[GridDomain], max_quads: int) -> int:
    domains = list(domains)
    total_cells = sum(domain.cells_per_edge * domain.cells_per_edge for domain in domains)
    if total_cells <= max_quads:
        return 1
    return max(1, int(np.ceil(np.sqrt(total_cells / max_quads))))


def _domain_quads(
    domain: GridDomain,
    points: np.ndarray,
    step: int = 1,
) -> Iterable[np.ndarray]:
    grid = domain.point_grid
    edge = domain.cells_per_edge
    for row in range(0, edge, step):
        row2 = min(row + step, edge)
        for col in range(0, edge, step):
            col2 = min(col + step, edge)
            yield points[
                [
                    grid[row, col],
                    grid[row, col2],
                    grid[row2, col2],
                    grid[row2, col],
                ]
            ]


def _sample_axis_indices(size: int, max_count: int) -> list[int]:
    if size <= max_count:
        return list(range(size))
    indices = np.linspace(0, size - 1, max_count, dtype=int)
    return sorted(set(int(index) for index in indices))


def _panel_grid_lines(
    panel_grids: list[np.ndarray],
    points: np.ndarray,
    max_lines_per_axis: int = 36,
) -> Iterable[np.ndarray]:
    for panel_grid in panel_grids:
        rows = _sample_axis_indices(panel_grid.shape[0], max_lines_per_axis)
        cols = _sample_axis_indices(panel_grid.shape[1], max_lines_per_axis)
        for row in rows:
            yield points[panel_grid[row, :]]
        for col in cols:
            yield points[panel_grid[:, col]]


def _domain_grid_lines(domain: GridDomain, points: np.ndarray) -> Iterable[np.ndarray]:
    grid = domain.point_grid
    for row in range(grid.shape[0]):
        yield points[grid[row, :]]
    for col in range(grid.shape[1]):
        yield points[grid[:, col]]


def _domain_boundary_lines(domain: GridDomain, points: np.ndarray) -> list[np.ndarray]:
    grid = domain.point_grid
    return [
        points[grid[0, :]],
        points[grid[-1, :]],
        points[grid[:, 0]],
        points[grid[:, -1]],
    ]


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="oasis-grid",
        description="Explore icosahedral and cube-sphere spherical grids.",
    )
    parser.add_argument(
        "--grid",
        default="icosahedral",
        help="Grid type: icosahedral, ico, cube-sphere, or cube.",
    )
    parser.add_argument(
        "--glevel",
        type=int,
        default=4,
        help="Refinement level; each panel edge has 2**glevel cells.",
    )
    parser.add_argument(
        "--divide-face",
        type=int,
        default=2,
        help="Domain split level; each panel has 4**divide_face domains.",
    )
    parser.add_argument(
        "--projection",
        choices=("equiangular", "normalized"),
        default="equiangular",
        help="Cube-sphere projection.",
    )
    parser.add_argument(
        "--focus-domain",
        type=int,
        default=None,
        help="Start with a specific domain highlighted.",
    )
    parser.add_argument(
        "--surface-alpha",
        type=float,
        default=0.50,
        help="Grid opacity from 0.0 transparent to 1.0 opaque.",
    )
    parser.add_argument(
        "--zoom",
        type=float,
        default=1.0,
        help="Initial view zoom. Values above 1 zoom in; below 1 zoom out.",
    )
    parser.add_argument(
        "--save",
        default=None,
        help="Optional PNG path for the initial view.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open the interactive window; useful with --save.",
    )
    parser.add_argument(
        "--postprocess",
        action="store_true",
        help="Print grid summary and a concise MPI-domain table.",
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print grid/decomposition summary as JSON.",
    )
    parser.add_argument(
        "--list-domains",
        action="store_true",
        help="Print the MPI-domain/rank table.",
    )
    parser.add_argument(
        "--domain-info",
        type=int,
        default=None,
        help="Print detailed JSON metadata for one MPI domain/rank.",
    )
    parser.add_argument(
        "--export-npz",
        default=None,
        help="Write postprocessing arrays to a compressed NPZ file.",
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=40,
        help="Maximum domain rows to print in tabular CLI output.",
    )
    return parser.parse_args()


def main() -> GridModel:
    args = _parse_args()
    postprocessing_requested = any(
        [
            args.postprocess,
            args.summary,
            args.list_domains,
            args.domain_info is not None,
            args.export_npz is not None,
        ]
    )

    needs_explorer = args.save is not None or (
        not args.no_show and not postprocessing_requested
    )
    if needs_explorer:
        explorer = GridExplorer(
            grid=args.grid,
            glevel=args.glevel,
            divide_face=args.divide_face,
            projection=args.projection,
            focus_domain=args.focus_domain,
            surface_alpha=args.surface_alpha,
            zoom=args.zoom,
        )
        model = explorer.model
        if args.save:
            explorer.save(args.save)
        if not args.no_show and not postprocessing_requested:
            explorer.show()
    else:
        model = build_grid_model(
            args.grid,
            args.glevel,
            args.divide_face,
            projection=args.projection,
        )

    if args.summary or args.postprocess:
        print(json.dumps(_json_ready(postprocess_grid(model)), indent=2))
    if args.list_domains or args.postprocess:
        print(_format_domain_table(postprocess_domains(model), max_rows=args.max_rows))
    if args.domain_info is not None:
        print(json.dumps(_json_ready(postprocess_domain(model, args.domain_info)), indent=2))
    if args.export_npz:
        output_path = export_postprocessing_npz(model, args.export_npz)
        print(f"Wrote {output_path}")

    return model


def _format_domain_table(rows: list[dict], max_rows: int = 40) -> str:
    if not rows:
        return "No domains."

    headers = ["rank", "panel", "tx", "ty", "cells", "points", "neighbors"]
    lines = [
        f"{headers[0]:>5}  {headers[1]:>8}  {headers[2]:>3}  {headers[3]:>3}  "
        f"{headers[4]:>5}  {headers[5]:>6}  {headers[6]}"
    ]
    displayed = rows[:max(0, max_rows)]
    for row in displayed:
        neighbors = ",".join(str(value) for value in row.get("neighbors", []))
        lines.append(
            f"{row['rank']:5d}  {row['panel_name'][:8]:>8}  "
            f"{row['tile_x']:3d}  {row['tile_y']:3d}  "
            f"{row['cell_count']:5d}  {row['point_count']:6d}  {neighbors}"
        )
    remaining = len(rows) - len(displayed)
    if remaining > 0:
        lines.append(f"... {remaining} more domains. Use --max-rows to print more.")
    return "\n".join(lines)


def _json_ready(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    return value


if __name__ == "__main__":
    main()

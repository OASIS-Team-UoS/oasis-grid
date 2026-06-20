import numpy as np

from oasis_grid import (
    build_domain_adjacency,
    build_grid_model,
    build_mpi_domain_map,
    extract_domain_mesh,
    generate_oasis_icosahedral_points,
    postprocess_domain,
    postprocess_grid,
)


def test_icosahedral_grid_counts_and_unit_sphere():
    model = build_grid_model("icosahedral", glevel=3, divide_face=1)

    assert model.panel_count == 10
    assert model.domain_count == 40
    assert model.total_cells == 640
    assert np.allclose(np.linalg.norm(model.points, axis=1), 1.0)


def test_cube_sphere_grid_counts_and_unit_sphere():
    model = build_grid_model("cube-sphere", glevel=3, divide_face=1)

    assert model.panel_count == 6
    assert model.domain_count == 24
    assert model.total_cells == 384
    assert np.allclose(np.linalg.norm(model.points, axis=1), 1.0)


def test_mpi_rank_map_matches_domain_order():
    model = build_grid_model("icosahedral", glevel=3, divide_face=1)
    rank_map = build_mpi_domain_map(model)

    assert rank_map.shape == (10, 2, 2)
    assert rank_map[0, 0, 0] == 0
    assert rank_map[0, 0, 1] == 1
    assert rank_map[0, 1, 0] == 2
    assert rank_map[1, 0, 0] == 4


def test_domain_postprocessing_and_mesh_shapes():
    model = build_grid_model("cube-sphere", glevel=3, divide_face=1)
    domain = postprocess_domain(model, 7)
    mesh = extract_domain_mesh(model, 7)

    assert domain["rank"] == 7
    assert domain["cells_per_edge"] == 4
    assert domain["cell_count"] == 16
    assert len(domain["neighbors"]) == 4
    assert mesh["point_grid"].shape == (5, 5)
    assert mesh["cell_point_indices"].shape == (4, 4, 4)
    assert mesh["cell_centers_xyz"].shape == (4, 4, 3)


def test_domain_adjacency_is_symmetric():
    model = build_grid_model("cube-sphere", glevel=3, divide_face=1)
    adjacency = build_domain_adjacency(model)

    for rank, neighbors in adjacency.items():
        for neighbor in neighbors:
            assert rank in adjacency[neighbor]


def test_oasis_icosahedral_generator_compatibility():
    xyz, rhomb, panel_grids = generate_oasis_icosahedral_points(
        3,
        return_coarse=True,
    )

    assert xyz.shape == (782, 3)
    assert rhomb.shape == (10, 64)
    assert panel_grids.shape == (10, 9, 9)
    assert np.allclose(np.linalg.norm(xyz, axis=1), 1.0)


def test_postprocess_summary_is_serializable_and_consistent():
    model = build_grid_model("icosahedral", glevel=3, divide_face=1)
    summary = postprocess_grid(model)

    assert summary["grid"] == "icosahedral"
    assert summary["domain_count"] == model.domain_count
    assert summary["total_cells"] == model.total_cells
    assert summary["unit_sphere_max_error"] < 1.0e-12

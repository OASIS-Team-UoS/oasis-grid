"""Export a coarse native OASIS dual without building or running the model.

Compile only three unchanged, extracted geometry methods with a tiny class
declaration. Connectivity comes from oasis-grid's ten rhombi, checked against
the native coordinates. No MPI, CUDA, solver, spring integration, or model
initialization is invoked. The generated C++ excerpt retains its GPL notice.
"""

from __future__ import annotations

import hashlib
import io
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np

from oasis_grid import build_grid_model


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def method(source, name):
    """Extract one complete C++ definition, ignoring braces in comments/strings."""
    start = source.index("void Icogrid::" + name + "(")
    opening = source.index("{", start)
    depth = 0
    tokens = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|[{}]', re.S)
    for match in tokens.finditer(source, opening):
        if match.group() == "{":
            depth += 1
        elif match.group() == "}":
            depth -= 1
            if depth == 0:
                return source[start:match.end()]
    raise ValueError(f"Unclosed native method: {name}")


PRELUDE = r'''// Extracted OASIS geometry: SPDX-License-Identifier: GPL-3.0-or-later
// Generated mesh-only wrapper; original function bodies are unchanged.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using real_t = double;
using realvector = std::vector<double>;
using intvector = std::vector<int>;
class Icogrid {
public:
 void sphere_ico(realvector&, int, int, int, int, int, int, intvector&, int);
 void find_qpoints_thor(const intvector&, const realvector&, const intvector&, realvector&) const;
 void relocate_centres_thor(const realvector&, const intvector&, realvector&) const;
};
'''

MAIN = r'''
int main(int argc, char** argv) {
 if (argc != 3) return 2;
 const int level = std::stoi(argv[1]);
 const int n = 1 << level, count = 10*n*n+2;
 Icogrid grid;
 realvector xyz(3*count), q;
 intvector pent(12), rings(6*count, -1);
 grid.sphere_ico(xyz, level, n*n, n, n*n, 1, 1, pent, count);
 std::cout << std::setprecision(17) << std::scientific;
 if (std::string(argv[2]) == "centres") {
   for (int i=0; i<count; ++i)
     std::cout << xyz[3*i] << ' ' << xyz[3*i+1] << ' ' << xyz[3*i+2] << '\n';
   return 0;
 }
 for (auto& v : rings) if (!(std::cin >> v)) return 3;
 grid.find_qpoints_thor(rings, xyz, pent, q);
 realvector seed(xyz);
 grid.relocate_centres_thor(q, pent, xyz);
 for (int i=0; i<count; ++i) {
   for (int k=0; k<3; ++k) std::cout << seed[3*i+k] << ' ';
   for (int k=0; k<3; ++k) std::cout << xyz[3*i+k] << ' ';
   for (int k=0; k<18; ++k) std::cout << q[18*i+k] << ' ';
   std::cout << '\n';
 }
}
'''


def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def deduplicate(points):
    """Deduplicate panel seams while preserving an original, unrounded point."""
    _, first, inverse = np.unique(np.round(points, 12), axis=0,
                                   return_index=True, return_inverse=True)
    return points[first], inverse


def panel_triangles(model, inverse):
    # Native neighbor deltas include (+1,+1), so the rhombus diagonal is
    # (i,j)--(i+1,j+1), never the other diagonal.
    triangles = []
    for panel in model.panel_grids:
        g = inverse[panel]
        for i in range(len(g)-1):
            for j in range(len(g)-1):
                a, b, c, d = g[i,j], g[i,j+1], g[i+1,j+1], g[i+1,j]
                triangles.extend([(a,b,c), (a,c,d)])
    return np.asarray(triangles)


def area_centroid(triangle):
    """Independent spherical area-moment formula, used only for validation."""
    a = triangle
    b = np.roll(a, -1, axis=0)
    cross = np.cross(a, b)
    lengths = np.linalg.norm(cross, axis=1)
    angles = np.arctan2(lengths, np.sum(a*b, axis=1))
    moment = np.sum(cross * (angles/lengths)[:,None], axis=0)
    if np.dot(moment, a.sum(axis=0)) < 0:
        moment = -moment
    return unit(moment)


def export_native(checkout, level, compiler="c++"):
    checkout = Path(checkout).resolve()
    build_path = checkout / "src/grid/grid_build.cpp"
    q_path = checkout / "src/grid/grid_pointsq.cpp"
    build, qsource = build_path.read_text(), q_path.read_text()
    helpers = build[build.index("namespace {"):build.index("// Icosahedral grid construction")]
    excerpt = PRELUDE + helpers + "\n" + method(build, "sphere_ico") + "\n"
    excerpt += method(qsource, "find_qpoints_thor") + "\n"
    excerpt += method(qsource, "relocate_centres_thor") + "\n" + MAIN

    model = build_grid_model("icosahedral", level, 0)
    py_points, inverse = deduplicate(model.points)
    triangles_py = panel_triangles(model, inverse)
    with tempfile.TemporaryDirectory(prefix="oasis-figure-mesh-") as tmp:
        cpp = Path(tmp) / "mesh.cpp"
        executable = Path(tmp) / "mesh"
        cpp.write_text(excerpt)
        subprocess.run([compiler, "-std=c++17", "-O1", str(cpp), "-o", str(executable)],
                       check=True, capture_output=True, text=True)
        command = [str(executable), str(level)]
        seed = np.loadtxt(io.StringIO(subprocess.run(command+["centres"], check=True,
                             capture_output=True, text=True).stdout))
        # Match independently generated native and Python coordinates bijectively.
        distance = np.linalg.norm(py_points[:,None,:]-seed[None,:,:], axis=2)
        mapping = distance.argmin(axis=1)
        coordinate_error = float(distance[np.arange(len(mapping)), mapping].max())
        assert len(np.unique(mapping)) == len(seed) == 10*4**level+2
        assert coordinate_error < 2e-13, coordinate_error
        triangles = mapping[triangles_py]
        neighbors = [set() for _ in seed]
        for a,b,c in triangles:
            neighbors[a].update([b,c])
            neighbors[b].update([a,c])
            neighbors[c].update([a,b])
        degrees = np.array([len(n) for n in neighbors])
        assert dict(zip(*np.unique(degrees, return_counts=True))) == {5:12, 6:len(seed)-12}
        assert np.array_equal(np.flatnonzero(degrees == 5),
                              np.r_[np.arange(10)*4**level,len(seed)-2,len(seed)-1])
        # Native methods require a clockwise ring viewed from outside.
        rings = np.full((len(seed),6), -1, dtype=int)
        for i,(p,ns) in enumerate(zip(seed,neighbors)):
            axis = np.eye(3)[np.argmin(abs(p))]
            east = unit(np.cross(axis,p))
            north = np.cross(p,east)
            ns = np.array(sorted(ns))
            angle = np.arctan2(seed[ns]@north, seed[ns]@east)
            rings[i,:len(ns)] = ns[np.argsort(-angle)]
        ring_text = " ".join(map(str,rings.ravel()))
        raw = np.loadtxt(io.StringIO(subprocess.run(command+["corners"], input=ring_text,
                              check=True, capture_output=True, text=True).stdout))

    centres = raw[:,3:6]
    q = raw[:,6:].reshape(-1,6,3)
    # One native vertex per construction triangle. This key also validates
    # shared q-points across all three incident cells without rounding them.
    vertex_by_triangle, vertices, cells = {}, [], []
    sharing_error, formula_error, relocation_error = 0., 0., 0.
    triangle_keys = {tuple(sorted(t)) for t in triangles}
    for i,degree in enumerate(degrees):
        cell = []
        for k in range(degree):
            tri = (i, rings[i,k], rings[i,(k+1)%degree])
            key = tuple(sorted(tri))
            assert key in triangle_keys
            if key not in vertex_by_triangle:
                vertex_by_triangle[key] = len(vertices)
                vertices.append(q[i,k])
            vi = vertex_by_triangle[key]
            sharing_error = max(sharing_error, float(np.linalg.norm(q[i,k]-vertices[vi])))
            formula_error = max(formula_error, float(np.linalg.norm(q[i,k]-area_centroid(seed[list(tri)]))))
            cell.append(vi)
        assert np.all(q[i,degree:] == 0)
        relocation_error = max(relocation_error, float(np.linalg.norm(centres[i]-area_centroid(q[i,:degree]))))
        cells.append(cell[::-1])  # outward CCW for the renderer/topology checks
    assert max(sharing_error,formula_error,relocation_error) < 2e-12
    face_array = np.full((len(cells),6), -1, dtype=int)
    for i,c in enumerate(cells):
        face_array[i,:len(c)] = c
    data = dict(vertices=np.array(vertices), faces=face_array, degrees=degrees,
                centres=centres, seed_centres=seed, neighbors=rings,
                construction_triangles=triangles)
    provenance = dict(
        checkout=str(checkout), glevel=level, radius=1.0,
        spring_relaxation=False, centre_relocation="Icogrid::relocate_centres_thor",
        generator="Icogrid::sphere_ico + Icogrid::find_qpoints_thor (unchanged native C++ bodies)",
        compiler=subprocess.run([compiler,"--version"], check=True, capture_output=True,
                                text=True).stdout.splitlines()[0],
        compile_flags=["-std=c++17","-O1"],
        sources={str(p.relative_to(checkout)):sha256(p) for p in [build_path,q_path]},
        native_python_coordinate_max_error=coordinate_error,
        shared_qpoint_max_error=sharing_error,
        independent_area_centroid_max_error=formula_error,
        independent_relocated_centre_max_error=relocation_error,
        constructor_minimum_level=4,
        coarse_mesh_note=f"Level {level} uses the native mesh methods directly; it is not a supported full-model run configuration.")
    return data, provenance, excerpt

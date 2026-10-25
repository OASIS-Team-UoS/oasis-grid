#!/usr/bin/env python3
"""Generate and validate the review's vector geometry/staggering figure.

From the repository root: python3 tools/figures/grid_geometry_and_staggering/generate.py
Dependencies follow requirements.txt; a C++17 compiler and sibling OASIS
checkout are needed only to refresh the native mesh. See README.md for the
portable, checksum-verified mesh reuse mode. No simulation is performed.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.text import Text
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from oasis_grid import build_grid_model
from native_mesh import export_native, sha256, unit


# ---------------------------------------------------------------------------
# Mesh construction and global geometry checks
# ---------------------------------------------------------------------------

def ordered_faces(points, faces):
    """Preserve connectivity while putting every face in outward CCW order."""
    result = []
    for face in faces:
        f = list(map(int, face))
        p = points[f]
        if np.dot(np.cross(p, np.roll(p,-1,axis=0)).sum(axis=0),p.sum(axis=0)) < 0:
            f.reverse()
        result.append(f)
    return result


def reference_latlon(nlon, nlat):
    """Regular latitude/longitude cells, including single shared pole vertices.

    Small-circle latitude edges are rendered and integrated analytically.
    This reference generator is a figure utility, not an OASIS backend.
    """
    points = [[0.,0.,-1.]]
    for j in range(1,nlat):
        lat = -np.pi/2 + j*np.pi/nlat
        for k in range(nlon):
            lon = k*2*np.pi/nlon
            points.append([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    points.append([0.,0.,1.])
    points = np.array(points)
    def index(j,k):
        return 0 if j == 0 else len(points)-1 if j == nlat else 1+(j-1)*nlon+k%nlon
    faces, areas = [], []
    for j in range(nlat):
        for k in range(nlon):
            face = list(dict.fromkeys([index(j,k),index(j,k+1),index(j+1,k+1),index(j+1,k)]))
            faces.append(face)
            areas.append((2*np.pi/nlon)*(np.sin(-np.pi/2+(j+1)*np.pi/nlat)-np.sin(-np.pi/2+j*np.pi/nlat)))
    return points, ordered_faces(points,faces), np.array(areas)


def cube_mesh(level, projection):
    model = build_grid_model("cube-sphere",level,0,projection=projection)
    faces, seams = [], set()
    for g in model.panel_grids:
        n = len(g)-1
        for i in range(n):
            for j in range(n):
                faces.append([g[i,j],g[i,j+1],g[i+1,j+1],g[i+1,j]])
        for line in [g[0,:],g[-1,:],g[:,0],g[:,-1]]:
            seams.update(tuple(sorted((int(a),int(b)))) for a,b in zip(line[:-1],line[1:]))
    return model.points, ordered_faces(model.points,faces), seams


def spherical_triangle_area(a,b,c):
    return 2*np.arctan2(np.linalg.det(np.array([a,b,c])),
                        1+np.dot(a,b)+np.dot(b,c)+np.dot(c,a))


def polygon_areas(points, faces):
    areas = []
    for f in faces:
        p = points[f]
        centre = unit(p.sum(axis=0))
        area = [spherical_triangle_area(centre,a,b) for a,b in zip(p,np.roll(p,-1,axis=0))]
        assert min(area) > 0, "Inconsistent or non-convex face ordering"
        areas.append(sum(area))
    return np.array(areas)


def check_mesh(points, faces, config, areas=None):
    assert np.isfinite(points).all() and points.ndim == 2 and points.shape[1] == 3
    sphere_error = float(abs(np.linalg.norm(points,axis=1)-1).max())
    assert sphere_error < config["coordinate_tolerance"]
    assert len(np.unique(np.round(points,11),axis=0)) == len(points)
    directed, edges, adjacency, used = Counter(), Counter(), defaultdict(set), set()
    for face in faces:
        assert len(face) >= 3 and len(set(face)) == len(face)
        assert min(face) >= 0 and max(face) < len(points)
        used.update(face)
        for a,b in zip(face,face[1:]+face[:1]):
            directed[a,b] += 1
            edges[tuple(sorted((a,b)))] += 1
            adjacency[a].add(b)
            adjacency[b].add(a)
            assert np.dot(points[a],points[b]) > 0, "Ambiguous long spherical edge"
    assert used == set(range(len(points)))
    assert set(edges.values()) == {2}, "An edge does not have exactly two incident cells"
    assert all(directed[a,b] == directed[b,a] == 1 for a,b in edges)
    seen, pending = set(), [0]
    while pending:
        i = pending.pop()
        if i not in seen:
            seen.add(i)
            pending.extend(adjacency[i]-seen)
    assert len(seen) == len(points)
    euler = len(points)-len(edges)+len(faces)
    assert euler == 2
    # This also checks outward face winding for the latitude/longitude caps.
    gc_areas = polygon_areas(points,faces)
    areas = gc_areas if areas is None else areas
    assert (areas > 0).all()
    total_area = float(areas.sum())
    error = abs(total_area-4*np.pi)
    assert error < config["area_tolerance_steradians"], error
    checks = dict(vertices=len(points),edges=len(edges),cells=len(faces),
                  euler_characteristic=euler,connected=True,
                  edge_incidence_histogram=dict(Counter(edges.values())),
                  all_face_windings_outward=True,opposite_shared_edge_directions=True,
                  unit_radius_max_error=sphere_error,total_area_steradians=total_area,
                  total_area_error_steradians=error,
                  area_method="oriented great-circle triangle solid angles" if areas is gc_areas else
                              "analytic longitude_width * (sin north - sin south)",
                  polygon_degree_histogram=dict(Counter(map(len,faces))))
    return checks, sorted(edges)


def check_native_centres(native,faces,config):
    centres,seed,rings = [native[k] for k in ["centres","seed_centres","neighbors"]]
    assert len(centres) == len(faces) == len(seed) == len(rings)
    assert np.isfinite(centres).all() and np.isfinite(seed).all()
    assert max(abs(np.linalg.norm(centres,axis=1)-1)) < config["coordinate_tolerance"]
    inside_margin = 1.
    for i,(centre,face) in enumerate(zip(centres,faces)):
        p = native["vertices"][face]
        inside_margin = min(inside_margin,float((np.cross(p,np.roll(p,-1,axis=0))@centre).min()))
        neighbors = rings[i,:len(face)]
        assert len(set(neighbors)) == len(face)
        assert min(neighbors) >= 0 and max(neighbors) < len(centres)
        for j in neighbors:
            assert i in rings[j], "Neighbor connectivity must be reciprocal"
    assert inside_margin > 0
    return dict(relocated_centres_inside_native_polygons=True,
                minimum_inward_edge_plane_dot=inside_margin,reciprocal_neighbor_connectivity=True)


# ---------------------------------------------------------------------------
# Exact hemisphere cuts for both great circles and small circles
# ---------------------------------------------------------------------------

def camera(longitude, latitude):
    lon,lat = np.deg2rad([longitude,latitude])
    east = np.array([-np.sin(lon),np.cos(lon),0.])
    north = np.array([-np.sin(lat)*np.cos(lon),-np.sin(lat)*np.sin(lon),np.cos(lat)])
    view = np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    return np.stack([east,north,view])


def visible_circle(centre, cosine, sine, start, end, basis, step_deg):
    """Sample only analytically visible intervals of c+a*cos(t)+b*sin(t).

    All limb roots of A+B*cos(t)+C*sin(t)=0 are inserted before sampling.
    Thus a visible arc survives even if both original endpoints are hidden;
    no line is joined across a hidden interval and no limb chord is invented.
    """
    A,B,C = [float(v@basis[2]) for v in [centre,cosine,sine]]
    amplitude = np.hypot(B,C)
    cuts = [start,end]
    if amplitude > 1e-15 and abs(A) <= amplitude:
        phase = np.arctan2(C,B)
        offset = np.arccos(np.clip(-A/amplitude,-1,1))
        for root in [phase-offset,phase+offset]:
            for turn in range(int(np.floor((start-root)/(2*np.pi)))-1,
                              int(np.ceil((end-root)/(2*np.pi)))+2):
                t = root+turn*2*np.pi
                if start+1e-13 < t < end-1e-13:
                    cuts.append(t)
    cuts = sorted(set(cuts))
    result = []
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        mid = (lo+hi)/2
        if A+B*np.cos(mid)+C*np.sin(mid) < -1e-14:
            continue
        t = np.linspace(lo,hi,max(2,int(np.ceil(np.rad2deg(hi-lo)/step_deg))+1))
        xyz = centre+np.cos(t)[:,None]*cosine+np.sin(t)[:,None]*sine
        projected = xyz@basis.T
        assert projected[:,2].min() > -2e-12
        assert np.linalg.norm(projected[:,:2],axis=1).max() <= 1+2e-12
        result.append(projected[:,:2])
    return result


def great_circle(a,b,basis,step):
    cosine = float(np.clip(a@b,-1,1))
    angle = np.arccos(cosine)
    return visible_circle(np.zeros(3),a,unit(b-cosine*a),0,angle,basis,step)


def check_clipping(config):
    basis = camera(0,0)
    # Both hidden endpoints but a visible middle part of a latitude circle.
    arcs = visible_circle(np.zeros(3),np.array([1.,0.,0.]),np.array([0.,1.,0.]),
                          -2.,2.,basis,config["curve_step_deg"])
    assert len(arcs) == 1 and abs(arcs[0][0,0]+1) < 1e-12 and abs(arcs[0][-1,0]-1) < 1e-12
    # Hidden great-circle edge, straddling edge, and a true small circle.
    assert not great_circle(unit(np.array([-1.,1.,0.])),unit(np.array([-1.,0.,1.])),basis,.5)
    cross = great_circle(unit(np.array([1.,1.,0.])),unit(np.array([-1.,1.,0.])),basis,.5)
    assert len(cross) == 1 and abs(np.linalg.norm(cross[0][-1])-1) < 1e-12
    small = visible_circle(np.array([0.,0.,.5]),np.array([np.sqrt(.75),0.,0.]),
                           np.array([0.,np.sqrt(.75),0.]),0,2*np.pi,basis,.5)
    assert all(np.allclose(a[:,1],.5) for a in small)
    return dict(analytic_limb_intersections=True,hidden_endpoints_visible_middle=True,
                hidden_edges_omitted=True,straddling_edges_preserved=True,
                latitude_small_circles=True)


# ---------------------------------------------------------------------------
# Independent local staggering schematics and marker-location checks
# ---------------------------------------------------------------------------

def staggering_markers(n):
    centres = np.array([(i+.5,j+.5) for j in range(n) for i in range(n)])
    u = np.array([(i,j+.5) for j in range(n) for i in range(n+1)])
    v = np.array([(i+.5,j) for j in range(n+1) for i in range(n)])
    return {"A":dict(h=centres,u=centres.copy(),v=centres.copy()),
            "C":dict(h=centres,u=u,v=v),
            "Z":dict(h=centres,delta=centres.copy(),zeta=centres.copy())}


def check_markers(markers):
    a,c,z = [markers[k] for k in ["A","C","Z"]]
    assert np.array_equal(a["h"],a["u"]) and np.array_equal(a["h"],a["v"])
    assert np.array_equal(z["h"],z["delta"]) and np.array_equal(z["h"],z["zeta"])
    assert "u" not in z and "v" not in z
    assert (c["h"]%1 == .5).all()
    assert (c["u"][:,0]%1 == 0).all() and (c["u"][:,1]%1 == .5).all()
    assert (c["v"][:,0]%1 == .5).all() and (c["v"][:,1]%1 == 0).all()
    return dict(A_collocated=True,C_components_on_correct_faces=True,
                Z_collocated_no_stored_velocity=True,
                positions={k:{v:p.tolist() for v,p in values.items()} for k,values in markers.items()})


# ---------------------------------------------------------------------------
# Publication layout: all coordinates below are page millimetres
# ---------------------------------------------------------------------------

def page_text(fig,x,y,text,config,fontsize=None,**kwargs):
    return fig.text(x/config["width_mm"],y/config["height_mm"],text,
                    fontsize=fontsize or config["font_pt"],color=config["ink"],**kwargs)


def panel_axes(fig,x,y,size,config,limits):
    ax = fig.add_axes([x/config["width_mm"],y/config["height_mm"],
                       size/config["width_mm"],size/config["height_mm"]])
    ax.set(xlim=limits,ylim=limits,aspect="equal")
    ax.axis("off")
    return ax


def render_geometries(fig,meshes,config,basis):
    page_text(fig,8,142,"Horizontal geometry",config,config["group_font_pt"],weight="bold")
    labels = ["(a) Latitude-longitude","(b) Cubed sphere","(c) Icosahedral dual"]
    details = ["Regular latitude-longitude","Equiangular gnomonic","OASIS · pentagons outlined"]
    for k,(label,detail) in enumerate(zip(labels,details)):
        x = 7+60*k
        page_text(fig,x+25,134,label,config,config["panel_font_pt"],ha="center")
        ax = panel_axes(fig,x,80,config["sphere_axes_size_mm"],config,(-1.025,1.025))
        paths,weights = [],[]
        if k == 0:
            spec = config["latitude_longitude"]
            for lon in np.arange(spec["longitude_cells"])*2*np.pi/spec["longitude_cells"]:
                paths.extend(visible_circle(np.zeros(3),np.array([np.cos(lon),np.sin(lon),0.]),
                             np.array([0.,0.,1.]),-np.pi/2,np.pi/2,basis,config["curve_step_deg"]))
            for lat in np.linspace(-np.pi/2,np.pi/2,spec["latitude_cells"]+1)[1:-1]:
                paths.extend(visible_circle(np.array([0.,0.,np.sin(lat)]),
                             np.array([np.cos(lat),0.,0.]),np.array([0.,np.cos(lat),0.]),
                             0,2*np.pi,basis,config["curve_step_deg"]))
            weights = [config["mesh_line_pt"]]*len(paths)
        else:
            m = meshes[k]
            pentagon_edges = set()
            if k == 2:
                for f in m["faces"]:
                    if len(f) == 5:
                        pentagon_edges.update(tuple(sorted((a,b))) for a,b in zip(f,f[1:]+f[:1]))
            for a,b in m["edges"]:
                segments = great_circle(m["points"][a],m["points"][b],basis,config["curve_step_deg"])
                width = config["mesh_line_pt"]
                if (a,b) in m.get("seams",set()):
                    width = config["boundary_line_pt"]
                if (a,b) in pentagon_edges:
                    width = config["pentagon_line_pt"]
                paths.extend(segments)
                weights.extend([width]*len(segments))
        collection = LineCollection(paths,colors=config["mesh_ink"],linewidths=weights,
                                    capstyle="round",joinstyle="round",rasterized=False)
        ax.add_collection(collection)
        ax.add_patch(Circle((0,0),1,fill=False,lw=config["sphere_outline_pt"],color=config["ink"]))
        page_text(fig,x+25,75,detail,config,ha="center")
    # A fine horizontal separator and whitespace make the two groups distinct.
    fig.add_artist(Line2D([8/config["width_mm"],172/config["width_mm"]],
                          [66/config["height_mm"]]*2,transform=fig.transFigure,
                          lw=.45,color="#bdbdbd"))


def draw_arrow(ax,anchor,component,config):
    delta = np.array([.27,0.]) if component == "u" else np.array([0.,.27])
    arrow = FancyArrowPatch(anchor,np.asarray(anchor)+delta,arrowstyle="-|>",
                            mutation_scale=5.7,linewidth=config["arrow_line_pt"],
                            shrinkA=0,shrinkB=0,color=config["ink"],zorder=3)
    ax.add_patch(arrow)


def render_staggering(fig,markers,config):
    page_text(fig,8,58,"Variable staggering",config,config["group_font_pt"],weight="bold")
    page_text(fig,172,58,"Independent local schematics",config,ha="right")
    n = config["schematic_cells_per_axis"]
    annotations = ["Scalar and velocity\ncollocated","Normal velocity\nat cell faces",
                   "Scalar, divergence and\nvorticity collocated"]
    for k,kind in enumerate(["A","C","Z"]):
        x = 7+60*k
        page_text(fig,x+25,50,f"({chr(100+k)}) {kind} grid",config,config["panel_font_pt"],ha="center")
        ax = panel_axes(fig,x+8.5,14,config["schematic_axes_size_mm"],config,(-.24,n+.24))
        for i in range(n+1):
            ax.plot([0,n],[i,i],color=config["schematic_ink"],lw=config["schematic_line_pt"],zorder=1)
            ax.plot([i,i],[0,n],color=config["schematic_ink"],lw=config["schematic_line_pt"],zorder=1)
        m = markers[kind]
        ax.plot(*m["h"].T,"o",ms=config["dot_size_pt"],color=config["ink"],zorder=4)
        for component in ["u","v"]:
            for anchor in m.get(component,[]):
                draw_arrow(ax,anchor,component,config)
        # One exemplar label per placement/component; only labels are offset.
        if kind == "A":
            ax.annotate("$h$",(.5,.5),xytext=(-2,-2),textcoords="offset points",ha="right",va="top")
            ax.annotate("$u$",(.635,.5),xytext=(0,-2),textcoords="offset points",ha="center",va="top")
            ax.annotate("$v$",(.5,.70),xytext=(-2,0),textcoords="offset points",ha="right",va="center")
        elif kind == "C":
            ax.annotate("$h$",(.5,.5),xytext=(2,2),textcoords="offset points",ha="left",va="bottom")
            ax.annotate("$u$",(1.135,.5),xytext=(0,-2),textcoords="offset points",ha="center",va="top")
            ax.annotate("$v$",(.5,1.135),xytext=(-2,0),textcoords="offset points",ha="right",va="center")
        else:
            ax.text(.61,.58,r"$h,\,\delta,\,\zeta$",ha="left",va="bottom")
        page_text(fig,x+25,12,annotations[k],config,ha="center",va="top",linespacing=1.25)
    page_text(fig,152,0.4,"Velocity reconstructed",config,ha="center")


def shared_legend(fig,config):
    # One neutral shared legend for both A and C arrows and all scalar dots.
    y = -8
    page_text(fig,14,y,"●  Stored scalar(s)",config,va="center")
    fig.add_artist(FancyArrowPatch((67/config["width_mm"],y/config["height_mm"]),
                                   (73/config["width_mm"],y/config["height_mm"]),
                                   transform=fig.transFigure,arrowstyle="-|>",mutation_scale=6,
                                   lw=config["arrow_line_pt"],color=config["ink"]))
    page_text(fig,76,y,"Stored velocity component",config,va="center")


def text_check(fig, config):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    boxes = [(t,t.get_window_extent(renderer)) for t in texts]
    # Layout is cropped below. Test actual nonempty text boxes for overlaps.
    overlaps = []
    for i,(a,ba) in enumerate(boxes):
        for b,bb in boxes[i+1:]:
            if ba.overlaps(bb):
                overlaps.append([a.get_text(),b.get_text()])
    assert not overlaps, overlaps
    assert min(t.get_fontsize() for t in texts) >= 8
    return dict(text_overlap_pairs=overlaps,minimum_font_pt=min(t.get_fontsize() for t in texts))


# ---------------------------------------------------------------------------
# Provenance, protected exports and reproducible execution
# ---------------------------------------------------------------------------

def repository_state(path, inputs):
    def git(*args):
        return subprocess.run(["git","-C",str(path),*args],check=True,capture_output=True,text=True).stdout
    status = git("status","--porcelain=v1","--untracked-files=normal").splitlines()
    input_status = git("status","--porcelain=v1","--untracked-files=normal","--",*inputs).splitlines()
    return dict(path=str(Path(path).resolve()),revision=git("rev-parse","HEAD").strip(),
                dirty=bool(status),status_counts=dict(Counter(line[:2] for line in status)),
                input_status_porcelain=input_status,
                tracked_diff_sha256=hashlib.sha256(git("diff","HEAD","--binary").encode()).hexdigest())


def protected_outputs(output):
    manifest = output/"generated_files.json"
    previous = json.loads(manifest.read_text()) if manifest.exists() else {}
    for name,digest in previous.items():
        path = output/name
        if path.exists() and sha256(path) != digest:
            raise RuntimeError(f"Refusing to overwrite modified generated file: {path}")
    expected = [f"{stem}.{ext}" for stem in ["grid_geometry_and_staggering",
                "horizontal_grid_geometries","horizontal_variable_staggering"] for ext in ["pdf","svg","png"]]
    expected += ["metadata.json","native_mesh.npz","native_provenance.json","native_mesh_extract.cpp",
                 "generation_config.json"]
    for name in expected:
        if (output/name).exists() and name not in previous:
            raise RuntimeError(f"Refusing to overwrite unregistered output: {output/name}")
    return manifest


def export(fig,output,config):
    # Fixed physical crops preserve the exact same dimensions/fonts in all groups.
    # Composite includes the common legend below the schematic axes.
    from matplotlib.transforms import Bbox
    crops = {"grid_geometry_and_staggering":(0.,config["export_bottom_mm"],config["width_mm"],config["height_mm"]),
             "horizontal_grid_geometries":(0.,config["geometry_crop_bottom_mm"],config["width_mm"],config["height_mm"]),
             "horizontal_variable_staggering":(0.,config["export_bottom_mm"],config["width_mm"],config["staggering_crop_top_mm"])}
    saved,checks = [],{}
    for stem,(left,bottom,right,top) in crops.items():
        bbox = Bbox.from_extents(left/25.4,bottom/25.4,right/25.4,top/25.4)
        renderer = fig.canvas.get_renderer()
        for t in fig.findobj(Text):
            if t.get_visible() and t.get_text():
                box = t.get_window_extent(renderer).transformed(fig.dpi_scale_trans.inverted())
                if bbox.overlaps(box):
                    assert (box.x0 >= bbox.x0 and box.y0 >= bbox.y0 and
                            box.x1 <= bbox.x1 and box.y1 <= bbox.y1), (stem,t.get_text())
        for ext in ["pdf","svg","png"]:
            path = output/f"{stem}.{ext}"
            metadata = {"Creator":"oasis-grid review figure utility"}
            if ext == "pdf":
                metadata.update(CreationDate=None,ModDate=None)
            if ext == "svg":
                metadata.update(Date=None)
            fig.savefig(path,format=ext,dpi=config["preview_dpi"],bbox_inches=bbox,
                        facecolor="white",metadata=metadata)
            saved.append(path)
        svg = (output/f"{stem}.svg").read_text()
        assert "<image" not in svg and "<text" in svg
        checks[stem] = dict(width_mm=right-left,height_mm=top-bottom,
                           svg_has_editable_text=True,svg_raster_images=0,
                           no_partially_clipped_text=True,
                           preview_dpi=config["preview_dpi"])
    return saved,checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=HERE/"config.json")
    parser.add_argument("--output-dir",type=Path,default=HERE)
    parser.add_argument("--oasis-checkout",type=Path,default=ROOT.parent/"oasis")
    parser.add_argument("--compiler",default="c++")
    parser.add_argument("--reuse-native",action="store_true",
                        help="Reuse the checksum-verified native snapshot in the output directory")
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    assert config["icosahedral_dual"]["spring_relaxation"] is False
    assert 1 <= config["icosahedral_dual"]["glevel"] <= 4
    assert config["schematic_cells_per_axis"] == 2
    output = args.output_dir.resolve()
    output.mkdir(parents=True,exist_ok=True)
    manifest = protected_outputs(output)
    destination = output
    staging = tempfile.TemporaryDirectory(prefix="oasis-figure-exports-")
    output = Path(staging.name)
    grid_state = repository_state(ROOT,["oasis_grid.py","tools/figures/grid_geometry_and_staggering"])
    if args.reuse_native:
        native_provenance = json.loads((destination/"native_provenance.json").read_text())
        assert sha256(destination/"native_mesh.npz") == native_provenance["mesh_sha256"]
        assert sha256(destination/"native_mesh_extract.cpp") == native_provenance["excerpt_sha256"]
        assert native_provenance["glevel"] == config["icosahedral_dual"]["glevel"]
        native = dict(np.load(destination/"native_mesh.npz",allow_pickle=False))
        for name in ["native_mesh.npz","native_mesh_extract.cpp","native_provenance.json"]:
            shutil.copy2(destination/name,output/name)
        native_mode = "checksum-verified archived native mesh"
    else:
        native, native_provenance, excerpt = export_native(args.oasis_checkout,
                                            config["icosahedral_dual"]["glevel"],args.compiler)
        native_provenance["repository"] = repository_state(args.oasis_checkout,list(native_provenance["sources"]))
        np.savez_compressed(output/"native_mesh.npz",**native)
        (output/"native_mesh_extract.cpp").write_text(excerpt)
        native_provenance["mesh_sha256"] = sha256(output/"native_mesh.npz")
        native_provenance["excerpt_sha256"] = sha256(output/"native_mesh_extract.cpp")
        (output/"native_provenance.json").write_text(json.dumps(native_provenance,indent=2)+"\n")
        native_mode = "fresh mesh-only native C++ generation"
    spec = config["latitude_longitude"]
    ll_points,ll_faces,ll_areas = reference_latlon(spec["longitude_cells"],spec["latitude_cells"])
    cube_points,cube_faces,seams = cube_mesh(**dict(level=config["cubed_sphere"]["glevel"],
                                                   projection=config["cubed_sphere"]["projection"]))
    native_faces = [list(map(int,f[:degree])) for f,degree in zip(native["faces"],native["degrees"])]
    meshes = [dict(points=ll_points,faces=ll_faces),dict(points=cube_points,faces=cube_faces,seams=seams),
              dict(points=native["vertices"],faces=native_faces)]
    geometry_checks = {}
    for i,name in enumerate(["latitude_longitude","cubed_sphere","icosahedral_dual"]):
        m = meshes[i]
        checks,m["edges"] = check_mesh(m["points"],m["faces"],config,ll_areas if i == 0 else None)
        geometry_checks[name] = checks
    assert Counter(map(len,native_faces)) == {5:12,6:len(native_faces)-12}
    geometry_checks["icosahedral_dual"].update(check_native_centres(native,native_faces,config))
    basis = camera(config["camera_longitude_deg"],config["camera_latitude_deg"])
    visible_pentagons = 0
    for f in native_faces:
        if len(f) == 5 and any(great_circle(native["vertices"][a],native["vertices"][b],basis,.5)
                                  for a,b in zip(f,f[1:]+f[:1])):
            visible_pentagons += 1
    geometry_checks["icosahedral_dual"].update(pentagons_global=12,pentagons_with_visible_boundary=visible_pentagons)
    markers = staggering_markers(config["schematic_cells_per_axis"])
    marker_checks = check_markers(markers)
    clipping_checks = check_clipping(config)
    plt.rcParams.update({"font.family":config["font_family"],"font.size":config["font_pt"],
                         "mathtext.fontset":"dejavusans","pdf.fonttype":42,"ps.fonttype":42,
                         "svg.fonttype":"none","svg.hashsalt":"oasis-grid-review-figure-v1",
                         "savefig.facecolor":"white"})
    fig = plt.figure(figsize=(config["width_mm"]/25.4,config["height_mm"]/25.4),dpi=120,facecolor="white")
    render_geometries(fig,meshes,config,basis)
    render_staggering(fig,markers,config)
    shared_legend(fig,config)
    typography_checks = text_check(fig,config)
    saved,export_checks = export(fig,output,config)
    plt.close(fig)
    (output/"generation_config.json").write_text(json.dumps(config,indent=2)+"\n")
    metadata = dict(
        purpose="Horizontal grid geometry and variable staggering are separate modelling choices.",
        generation_command="python3 tools/figures/grid_geometry_and_staggering/generate.py",
        invocation=[sys.executable,*sys.argv],repository=grid_state,native_generation_mode=native_mode,
        native_provenance=native_provenance,
        source_sha256={str(p.resolve().relative_to(ROOT)) if p.resolve().is_relative_to(ROOT) else str(p.resolve()):sha256(p)
                       for p in [Path(__file__),HERE/"native_mesh.py",HERE/"audit_exports.py",args.config,
                                 HERE/"caption.tex",HERE/"README.md",ROOT/"oasis_grid.py"]},
        dependencies=dict(python=platform.python_version(),numpy=np.__version__,matplotlib=matplotlib.__version__),
        config=config,camera=dict(projection="orthographic",basis_rows_east_north_view=basis.tolist(),
                                 visible_hemisphere="dot(position,view) >= 0",sphere_diameter_mm=config["sphere_axes_size_mm"]/1.025),
        generators=dict(latitude_longitude="reference_latlon in generate.py; regular, not Gaussian",
                        cubed_sphere="oasis_grid.build_grid_model; equiangular gnomonic; explorer generator, not a runtime OASIS backend",
                        icosahedral_dual=native_provenance["generator"]),
        geometry_checks=geometry_checks,staggering_checks=marker_checks,clipping_checks=clipping_checks,
        typography_checks=typography_checks,export_checks=export_checks,
        limitations=["Mesh density is chosen for legibility; no claim of equal effective resolution.",
                     "No spring relaxation in the selected mesh; native polygon corners and centre relocation are retained.",
                     "Coarse native mesh-only methods bypass the full model's glevel>=4 constructor guard.",
                     "Checks validate this figure's construction, not dynamical accuracy or solver conservation.",
                     "SVG text is editable; PDF fonts are embedded; curves use vector polylines at the configured angular step."],
        bibliography="No manuscript .bib found in the inspected repositories; caption uses no invented citation keys.")
    (output/"metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    generated = saved+[output/p for p in ["metadata.json","native_mesh.npz","native_provenance.json",
                                         "native_mesh_extract.cpp","generation_config.json"]]
    for path in generated:
        shutil.copy2(path,destination/path.name)
    manifest.write_text(json.dumps({p.name:sha256(p) for p in generated},indent=2)+"\n")
    staging.cleanup()
    print(json.dumps({"output":str(destination),"cell_counts":{k:v["cells"] for k,v in geometry_checks.items()},
                      "visible_pentagons":visible_pentagons,"checks":"passed","exports":len(saved)},indent=2))


if __name__ == "__main__":
    main()

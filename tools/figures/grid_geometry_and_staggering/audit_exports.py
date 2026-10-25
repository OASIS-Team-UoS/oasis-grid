#!/usr/bin/env python3
"""Render PDF/SVG/PNG exports at 180 mm and check their vector/text content.

Optional dependencies are listed in requirements-audit.txt. Report JSON is
written to stdout; inspection images go to a new temporary directory unless
--render-dir is supplied. This audits the files; visual review is still needed.
"""

import argparse
import hashlib
import io
import json
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

import cairosvg
import numpy as np
from PIL import Image
import PIL
import pymupdf

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument("--render-dir",type=Path)
    args = parser.parse_args()
    output = args.output_dir
    renders = args.render_dir or Path(tempfile.mkdtemp(prefix="oasis-figure-inspection-"))
    renders.mkdir(parents=True,exist_ok=True)
    metadata = json.loads((output/"metadata.json").read_text())
    report = dict(export_sha256={},formats={},render_directory=str(renders),inspection_dpi=144,
                  versions=dict(pymupdf=pymupdf.VersionBind,cairosvg=cairosvg.__version__,pillow=PIL.__version__))
    namespace = {"s":"http://www.w3.org/2000/svg"}
    for stem,spec in metadata["export_checks"].items():
        results = {}
        for ext in ["pdf","svg","png"]:
            path = output/f"{stem}.{ext}"
            report["export_sha256"][path.name] = sha256(path)
            if ext == "pdf":
                doc = pymupdf.open(path)
                assert len(doc) == 1
                page = doc[0]
                assert len(page.get_images()) == 0
                assert abs(page.rect.width*25.4/72-spec["width_mm"]) < 1e-5
                assert abs(page.rect.height*25.4/72-spec["height_mm"]) < 1e-5
                fonts = page.get_fonts()
                assert fonts and all(doc.extract_font(f[0])[3] for f in fonts)
                spans = [s for b in page.get_text("dict")["blocks"] if "lines" in b
                         for line in b["lines"] for s in line["spans"] if s["text"].strip()]
                assert spans and all(s["size"] >= 8-1e-6 for s in spans)
                assert all(pymupdf.Rect(s["bbox"]) in page.rect for s in spans)
                pixmap = page.get_pixmap(matrix=pymupdf.Matrix(2,2),alpha=False)
                image = Image.frombytes("RGB",[pixmap.width,pixmap.height],pixmap.samples)
                results[ext] = dict(images=0,embedded_fonts=True,minimum_font_pt=min(s["size"] for s in spans),
                                    all_text_within_page=True)
            elif ext == "svg":
                tree = ET.parse(path)
                assert not tree.findall(".//s:image",namespace)
                assert tree.findall(".//s:text",namespace)
                assert "DejaVu Sans" in path.read_text()
                # CairoSVG understands the SVG CSS font shorthand; it avoids
                # the fallback serif fonts used by some limited SVG importers.
                image = Image.open(io.BytesIO(cairosvg.svg2png(url=str(path),dpi=144,
                                                background_color="white"))).convert("RGB")
                results[ext] = dict(images=0,editable_text=True,declared_font="DejaVu Sans")
            else:
                source = Image.open(path)
                assert all(abs(d-600) < .01 for d in source.info["dpi"])
                assert abs(source.width-spec["width_mm"]*600/25.4) < 1
                assert abs(source.height-spec["height_mm"]*600/25.4) < 1
                results[ext] = dict(pixels=list(source.size),dpi=list(source.info["dpi"]))
                image = source.convert("RGB").resize((round(spec["width_mm"]*144/25.4),
                            round(spec["height_mm"]*144/25.4)),Image.Resampling.LANCZOS)
            raster = np.asarray(image)
            assert np.max(abs(raster[:,:,0].astype(int)-raster[:,:,1])) <= 1
            assert np.max(abs(raster[:,:,1].astype(int)-raster[:,:,2])) <= 1
            assert (raster < 200).any()
            # No ink touches the page edges at the intended print dimensions.
            assert (raster[[0,-1],:,:] > 245).all() and (raster[:,[0,-1],:] > 245).all()
            image.save(renders/f"{stem}-{ext}.png")
        report["formats"][stem] = results
    report["checks"] = "passed"
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()

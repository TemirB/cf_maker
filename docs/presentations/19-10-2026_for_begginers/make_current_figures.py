"""Make PDF crops of the author's JPEG plots supplied on 2026-10-19.

The copied JPEG files are the primary sources. They are embedded byte for byte
and clipped in PDF; no points, curves, axes, uncertainties, or fit boxes are
redrawn. Run ``python make_current_figures.py`` from any working directory.

These current author-provided plots are independent of the archival 2025 PDF.
Unreported settings are deliberately left unknown in both manifests.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pymupdf


BASE = Path(__file__).resolve().parent
SOURCES = BASE / "sources" / "current_2026_10_19"
FIGURES = BASE / "figures"
POINTS_PER_PIXEL = 0.25
SUPPLIED_DATE = "2026-10-19"
SOURCE_NAMES = (
    "kt_c_all_graphs_neg.jpg",
    "c_all_graphs_neg.jpg",
    "kt_cfs_neg_0-10_0.15-0.25.jpg",
    "cfs_neg_0-10_0.0_0.2.jpg",
)
AUTHOR_ROOT = Path("/Users/temir-balinov")

# All rectangles are actual source-JPEG pixels, not resized preview pixels.
# The parameter panels retain their original axis labels and all plotted marks.
PARAMETER_RECTS = {
    "out": (60, 50, 1460, 770),
    "side": (1658, 50, 3058, 770),
    "long": (60, 836, 1460, 1556),
    "lambda": (1658, 836, 3058, 1556),
}
CF_RECTS = {
    "out": (0, 0, 799, 772),
    "side": (799, 0, 1598, 772),
    "long": (1598, 0, 2396, 772),
}
CF_HEADER_HEIGHT = 44
CF_AXIS_HEADER_RECTS = (
    (0, 0, 136, CF_HEADER_HEIGHT),
    (799, 0, 935, CF_HEADER_HEIGHT),
    (1598, 0, 1733, CF_HEADER_HEIGHT),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def copy_sources() -> dict[str, dict]:
    """Import once; subsequent builds use the frozen project copies only."""
    SOURCES.mkdir(parents=True, exist_ok=True)
    sources = {}
    for name in SOURCE_NAMES:
        original = AUTHOR_ROOT / name
        copied = SOURCES / name
        if not copied.exists():
            if not original.exists():
                raise FileNotFoundError(f"Missing source copy {copied}; import the author's file {original}")
            shutil.copyfile(original, copied)
            if sha256(original) != sha256(copied):
                raise RuntimeError(f"Source copy is not byte-identical: {name}")
        pixmap = pymupdf.Pixmap(str(copied))
        sources[name] = {
            "source": str(copied.relative_to(BASE)),
            "original_absolute_path": str(original),
            "sha256": sha256(copied),
            "bytes": copied.stat().st_size,
            "pixels": [pixmap.width, pixmap.height],
            "provenance_type": "current_author_provided",
            "supplied_date": SUPPLIED_DATE,
            "archival_2025_source": False,
        }
    return sources


def image_document(path: Path, size: list[int]) -> pymupdf.Document:
    document = pymupdf.open()
    page = document.new_page(width=size[0], height=size[1])
    page.insert_image(page.rect, stream=path.read_bytes())
    return document


def add_clip(
    page: pymupdf.Page,
    original: pymupdf.Document,
    clip: pymupdf.Rect,
    crop: pymupdf.Rect,
) -> None:
    target = pymupdf.Rect(
        (clip.x0 - crop.x0) * POINTS_PER_PIXEL,
        (clip.y0 - crop.y0) * POINTS_PER_PIXEL,
        (clip.x1 - crop.x0) * POINTS_PER_PIXEL,
        (clip.y1 - crop.y0) * POINTS_PER_PIXEL,
    )
    page.show_pdf_page(target, original, 0, clip=clip)


def save_crop(
    original: pymupdf.Document,
    source: dict,
    name: str,
    rect: tuple[int, int, int, int],
    description: str,
    conditions: dict,
    *,
    remove_cf_header: bool = False,
) -> dict:
    crop = pymupdf.Rect(rect)
    page_rect = original[0].rect
    if not page_rect.contains(crop):
        raise ValueError(f"Crop outside source image: {name}: {rect}")
    output = pymupdf.open()
    page = output.new_page(width=crop.width * POINTS_PER_PIXEL, height=crop.height * POINTS_PER_PIXEL)
    retained_rects = []
    if remove_cf_header:
        # Titles and the upper y-axis tick labels share a horizontal band.
        # Clip out only the technical-title portions, retaining the left axis
        # strips at their exact original coordinates. Nothing is painted over.
        body = crop & pymupdf.Rect(0, CF_HEADER_HEIGHT, page_rect.width, page_rect.height)
        retained_rects.append(list(body))
        add_clip(page, original, body, crop)
        for raw_strip in CF_AXIS_HEADER_RECTS:
            strip = crop & pymupdf.Rect(raw_strip)
            if not strip.is_empty:
                retained_rects.append(list(strip))
                add_clip(page, original, strip, crop)
    else:
        retained_rects.append(list(crop))
        add_clip(page, original, crop, crop)
    output.set_metadata({
        "title": description,
        "author": "Temir Balinov",
        "subject": "Current author-provided SMASH plots, supplied 2026-10-19. Original JPEG, PDF clipping only.",
        "creator": "make_current_figures.py / PyMuPDF",
    })
    destination = FIGURES / f"{name}.pdf"
    output.save(destination, garbage=4, deflate=True)
    output.close()
    # Verify that PDF creation retained the exact original JPEG bitstream.
    with pymupdf.open(destination) as check:
        image_hashes = {
            hashlib.sha256(check.extract_image(image[0])["image"]).hexdigest()
            for image in check[0].get_images(full=True)
        }
        if image_hashes != {source["sha256"]}:
            raise RuntimeError(f"Embedded JPEG changed in {destination}")
    return {
        "asset": str(destination.relative_to(BASE)),
        "source": source["source"],
        "source_sha256": source["sha256"],
        "provenance_type": "current_author_provided",
        "supplied_date": SUPPLIED_DATE,
        "archival_2025_source": False,
        "source_raster_pixels": source["pixels"],
        "crop_pixels": list(rect),
        "retained_source_rectangles_pixels": retained_rects,
        "technical_cf_header_removed": remove_cf_header,
        "header_band_pixels": [0, 0, source["pixels"][0], CF_HEADER_HEIGHT] if remove_cf_header else None,
        "pdf_points": [crop.width * POINTS_PER_PIXEL, crop.height * POINTS_PER_PIXEL],
        "description": description,
        "conditions": conditions,
        "method": "Exact original JPEG stream embedded; PDF clipping without resampling, retouching, redrawing, or error reconstruction",
        "embedded_jpeg_sha256_verified": True,
        "source_is_vector": False,
    }


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    sources = copy_sources()
    for name in SOURCE_NAMES[:2]:
        if sources[name]["pixels"] != [3196, 1572]:
            raise ValueError(f"Unexpected parameter-image dimensions: {name}")
    for name in SOURCE_NAMES[2:]:
        if sources[name]["pixels"] != [2396, 772]:
            raise ValueError(f"Unexpected correlation-image dimensions: {name}")
    common = {
        "model": "SMASH",
        "collision_system": "Au+Au",
        "sqrt_sNN_GeV": 3,
        "system_energy_model_evidence": "Explicitly confirmed by the author in the current session",
        "radii_reference_frame": "LCMS",
        "rapidity_reference_frame": "collision center-of-mass system",
        "reference_frame_evidence": "Explicit current-session author confirmation: Да: радиусы в LCMS, y в системе центра масс. This evidence is author confirmation, not independent verification of a code run.",
        "source_charge_tag": "neg",
        "charge_evidence": "Source filenames; the correlation-panel headers explicitly state charge=neg",
        "centrality_colors_percent": {"red": [0, 10], "blue": [10, 30], "pink": [30, 50], "green": [50, 80]},
        "centrality_color_evidence": "Explicitly confirmed by the author in the current session",
        "not_inferred_from_2025": True,
    }
    unknown = [
        "SMASH version and number of generated events",
        "Full single-particle momentum and pseudorapidity cuts",
        "A rapidity cut is not assigned to the kT plots",
        "kT integration range for rapidity-dependent panels",
        "Projection windows in the other q components",
        "Uncertainty estimator and meaning of the quoted fit errors",
        "Full fit implementation and final-state-interaction treatment",
    ]
    assets = []
    for variable, filename in (("kt", SOURCE_NAMES[0]), ("y", SOURCE_NAMES[1])):
        conditions = dict(common)
        conditions["horizontal_variable"] = "kT" if variable == "kt" else "rapidity"
        if variable == "y":
            conditions["rapidity_binning_confirmed"] = {"range": [-1, 1], "step": 0.2}
        else:
            conditions["kT_unit"] = "GeV/c is the usual project convention; the source axis prints only kt"
        source = sources[filename]
        with image_document(BASE / source["source"], source["pixels"]) as original:
            for panel, rect in PARAMETER_RECTS.items():
                assets.append(save_crop(original, source, f"current_2026_{variable}_{panel}", rect,
                                        f"Current author-provided {panel} parameter versus {conditions['horizontal_variable']}", conditions))
            assets.append(save_crop(original, source, f"current_2026_{variable}_all", (0, 0, 3196, 1572),
                                    f"Current author-provided four parameter panels versus {conditions['horizontal_variable']}", conditions))
    for variable, filename in (("kt", SOURCE_NAMES[2]), ("y", SOURCE_NAMES[3])):
        source = sources[filename]
        conditions = dict(common)
        conditions.update({
            "centrality_percent": [0, 10],
            "centrality_evidence": "Explicitly written in original correlation-panel headers",
            "source_charge_tag": "neg",
            "charge_evidence": "Explicitly written in original correlation-panel headers",
            "fit_box_retained_in_long_panel": True,
        })
        if variable == "kt":
            conditions["kT_GeV_per_c"] = [0.15, 0.25]
            conditions["kT_evidence"] = "Original CF header explicitly states kt=0.15-0.25; the header does not print a unit"
        else:
            conditions["rapidity"] = [0.0, 0.2]
            conditions["rapidity_evidence"] = "Original CF header explicitly states rapidity=[0.0,0.2]"
        with image_document(BASE / source["source"], source["pixels"]) as original:
            for panel, rect in CF_RECTS.items():
                assets.append(save_crop(original, source, f"current_2026_cf_{variable}_{panel}", rect,
                                        f"Current author-provided {panel} correlation projection, {variable} selection", conditions,
                                        remove_cf_header=True))
            assets.append(save_crop(original, source, f"current_2026_cf_{variable}_all", (0, 0, 2396, 772),
                                    f"Current author-provided three correlation projections, {variable} selection", conditions,
                                    remove_cf_header=True))
    source_manifest = {
        "description": "Current original author-provided JPEG plots; separate from the 2025 archival presentation",
        "supplied_date": SUPPLIED_DATE,
        "provenance_type": "current_author_provided",
        "byte_identical_source_copies": True,
        "source_files": list(sources.values()),
        "author_confirmed_common_conditions": common,
        "unknown_or_not_verified": unknown,
        "numerical_processing": "None; no digitization, new fits, uncertainty reconstruction, or transfer of 2025 settings",
    }
    write_json(SOURCES / "manifest.json", source_manifest)
    write_json(FIGURES / "current_2026_figure_manifest.json", {
        "description": "PDF clipping of current JPEG plots supplied by the author on 2026-10-19",
        "source_manifest": str((SOURCES / "manifest.json").relative_to(BASE)),
        "coordinate_system": "Original JPEG pixels, origin at upper left; rectangles [x0,y0,x1,y1]",
        "points_per_pixel": POINTS_PER_PIXEL,
        "unknown_or_not_verified": unknown,
        "assets": assets,
    })
    print(f"Prepared {len(sources)} byte-identical source JPEGs and {len(assets)} clipped PDF assets.")


if __name__ == "__main__":
    main()

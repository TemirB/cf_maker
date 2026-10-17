"""Crop research figures without redrawing or resampling their original data.

Run with the presentation's Python environment::

    python extract_research_figures.py

The source PDFs contain *raster* plots.  The resulting PDF assets preserve the
embedded original raster streams and use PDF clipping only; these are not newly
vectorized plots.  Composite crops contain separate, explicitly listed copies
of the original legend.  No data points, error bars, or axis labels are redrawn.
All pixel rectangles and conditions are recorded in research_figure_manifest.json.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pymupdf


BASE = Path(__file__).resolve().parent
FIGURES = BASE / "figures"
SOURCES = BASE / "sources"
POINTS_PER_PIXEL = 0.25

KRAEVA_URL = "https://inspirehep.net/files/6a9386d330d51f4723a2145d87465b26"
KRAEVA_DOI = "https://doi.org/10.3390/universe10040188"
BALINOV_URL = "http://web-static-data.msu.dubna.ru:8093/epp/kursov_4/2025/Balinov.pdf"


def image_document(source: pymupdf.Document, page_number: int, size: tuple[int, int], image_index: int):
    """Return the exact embedded bitmap placed on a pixel-coordinate PDF page."""
    page = source[page_number - 1]
    matches = [item for item in page.get_images(full=True) if item[2:4] == size]
    matches.sort(key=lambda item: page.get_image_rects(item[0])[0].x0)
    if image_index >= len(matches):
        raise ValueError(f"PDF page {page_number}: expected image {image_index} of size {size}, got {len(matches)}")
    image = source.extract_image(matches[image_index][0])
    document = pymupdf.open()
    image_page = document.new_page(width=size[0], height=size[1])
    image_page.insert_image(image_page.rect, stream=image["image"])
    return document


def save_crop(
    source: pymupdf.Document,
    metadata: dict,
    name: str,
    page_number: int,
    image_size: tuple[int, int],
    rect: tuple[int, int, int, int],
    description: str,
    conditions: str,
    *,
    image_index: int = 0,
    legend_rect: tuple[int, int, int, int] | None = None,
    legend_scale: float = 0.6,
    header_rect: tuple[int, int, int, int] | None = None,
    footer_rect: tuple[int, int, int, int] | None = None,
):
    """Create a PDF crop with optional original legend, header, or shared axis.

    Header/footer strips retain their original horizontal pixel coordinates and
    scale. They allow a panel of a shared-axis figure to remain self-contained.
    """
    original = image_document(source, page_number, image_size, image_index)
    clip = pymupdf.Rect(rect)
    gap = 12 if legend_rect else 0
    legend = pymupdf.Rect(legend_rect) if legend_rect else None
    header = pymupdf.Rect(header_rect) if header_rect else None
    footer = pymupdf.Rect(footer_rect) if footer_rect else None
    if any(strip and strip.width != clip.width for strip in (header, footer)):
        raise ValueError("Header/footer strips must have the same width as the main crop")
    width = clip.width + gap + (legend.width * legend_scale if legend else 0)
    header_height = header.height if header else 0
    footer_height = footer.height if footer else 0
    height = max(clip.height + header_height + footer_height,
                 legend.height * legend_scale if legend else 0)
    output = pymupdf.open()
    page = output.new_page(width=width * POINTS_PER_PIXEL, height=height * POINTS_PER_PIXEL)
    target = pymupdf.Rect(0, header_height * POINTS_PER_PIXEL,
                         clip.width * POINTS_PER_PIXEL,
                         (header_height + clip.height) * POINTS_PER_PIXEL)
    page.show_pdf_page(target, original, 0, clip=clip)
    if header:
        page.show_pdf_page(pymupdf.Rect(0, 0, header.width * POINTS_PER_PIXEL,
                                       header_height * POINTS_PER_PIXEL),
                           original, 0, clip=header)
    if footer:
        page.show_pdf_page(pymupdf.Rect(0, (header_height + clip.height) * POINTS_PER_PIXEL,
                                       footer.width * POINTS_PER_PIXEL,
                                       (header_height + clip.height + footer_height) * POINTS_PER_PIXEL),
                           original, 0, clip=footer)
    if legend:
        target = pymupdf.Rect(
            (clip.width + gap) * POINTS_PER_PIXEL,
            0,
            width * POINTS_PER_PIXEL,
            legend.height * legend_scale * POINTS_PER_PIXEL,
        )
        page.show_pdf_page(target, original, 0, clip=legend)
    output.set_metadata({
        "title": description,
        "subject": "Original embedded raster; PDF clipping only, no data redrawing. " + conditions,
        "creator": "extract_research_figures.py / PyMuPDF",
    })
    destination = FIGURES / f"{name}.pdf"
    output.save(destination, garbage=4, deflate=True)
    output.close()
    original.close()
    entry = metadata.copy()
    entry.update({
        "asset": str(destination.relative_to(BASE)),
        "source_pdf_page": page_number,
        "source_raster_pixels": image_size,
        "source_image_index_left_to_right": image_index,
        "crop_pixels": rect,
        "legend_crop_pixels": legend_rect,
        "legend_scale": legend_scale if legend else None,
        "header_crop_pixels": header_rect,
        "footer_crop_pixels": footer_rect,
        "pdf_points": [round(width * POINTS_PER_PIXEL, 2), round(height * POINTS_PER_PIXEL, 2)],
        "description": description,
        "conditions": conditions,
        "source_is_vector": False,
        "method": "Original embedded raster preserved; PDF clipping without resampling or redrawing",
    })
    return entry


def main():
    FIGURES.mkdir(exist_ok=True)
    kraeva_path = SOURCES / "kraeva_2024.pdf"
    balinov_path = SOURCES / "balinov_2025.pdf"
    if not kraeva_path.exists():
        raise FileNotFoundError(f"Download {KRAEVA_URL} to {kraeva_path} first")
    assets = []
    with pymupdf.open(kraeva_path) as kraeva, pymupdf.open(balinov_path) as balinov:
        kraeva_meta = {
            "source": "sources/kraeva_2024.pdf",
            "source_url": KRAEVA_URL,
            "doi": KRAEVA_DOI,
            "sha256": hashlib.sha256(kraeva_path.read_bytes()).hexdigest(),
            "credit": "A. Kraeva on behalf of STAR, Universe 10 (2024) 188, CC BY 4.0; original STAR Preliminary labels retained",
        }
        balinov_meta = {
            "source": "sources/balinov_2025.pdf",
            "source_url": BALINOV_URL,
            "sha256": hashlib.sha256(balinov_path.read_bytes()).hexdigest(),
            "credit": "T. M. Balinov, archival defense presentation (2025); original graphs and uncertainties unchanged",
            "limitations": "Historical model result. 2025 fit code and error algorithm unavailable; source formulas contain dimensional/sign ambiguities.",
        }

        def k(name, page, size, rect, description, conditions, **kwargs):
            assets.append(save_crop(kraeva, kraeva_meta, "kraeva_2024_" + name,
                                    page, size, rect, description, conditions, **kwargs))

        def b(name, page, size, rect, description, conditions, **kwargs):
            assets.append(save_crop(balinov, balinov_meta, "balinov_2025_" + name,
                                    page, size, rect, description, conditions, **kwargs))

        # Figure 2: detector particle identification.  Full bitmap preserves both
        # PID methods, signed p/q axis, particle labels, and STAR/energy headings.
        k("pid", 6, (598, 976), (0, 0, 598, 976),
          "Figure 2: TPC ionization and TOF mass-squared particle identification",
          "STAR Au+Au sqrt(sNN)=3 GeV, 0-80%; signed p/q; top dE/dx, bottom TOF m^2. "
          "Illustrates PID, not the femtoscopic source.")
        k("pid_tpc_dedx", 6, (598, 976), (0, 0, 598, 487),
          "Figure 2 top: TPC dE/dx PID with the original shared momentum axis",
          "STAR Au+Au sqrt(sNN)=3 GeV, 0-80%. dE/dx versus signed p/q; colored curves are species expectations. "
          "The common x-axis tick labels and p/q label are copied from the bottom of the same figure. "
          "This illustrates single-track PID, not a position measurement at the femtoscopic scale.",
          footer_rect=(0, 901, 598, 976))
        k("pid_tof_mass2", 6, (598, 976), (0, 487, 598, 976),
          "Figure 2 bottom: TOF mass-squared PID with the original STAR/energy header",
          "STAR Au+Au sqrt(sNN)=3 GeV, 0-80%. m^2 inferred from TOF and TPC momentum/path length versus signed p/q. "
          "The STAR/energy/centrality header is copied from the top of the same figure. "
          "Finite resolution allows negative reconstructed m^2 values; these are measurement estimates, not negative physical mass squared.",
          header_rect=(0, 0, 598, 72))

        # Figure 3: pair acceptance with the selected rapidity windows.
        k("acceptance", 7, (1798, 776), (0, 0, 1798, 776),
          "Figure 3: pion-pair acceptance and selected rapidity windows",
          "STAR Au+Au sqrt(sNN)=3 GeV, 0-80%; positive pairs left, negative right. "
          "Horizontal coordinate is pair rapidity in the CM; y_beam=1.05 is used for the lab-to-CM shift.")
        k("acceptance_positive", 7, (1798, 776), (0, 0, 900, 776),
          "Figure 3 left: positive-pion pair acceptance",
          "STAR Au+Au sqrt(sNN)=3 GeV, 0-80%; pi+pi+; pair CM rapidity versus kT. "
          "Dashed lines mark rapidity-analysis windows; this is detector acceptance.")

        correlation_conditions = (
            "STAR Preliminary, Au+Au sqrt(sNN)=3 GeV, 0-10%, 0.15<kT<0.25 GeV/c; "
            "no additional pair-rapidity constraint. Red pi+pi+, blue pi-pi-; curves are Bowler-Sinyukov fits. "
            "For each projection the other q components are integrated over +/-0.05 GeV/c. "
            "Error bars are statistical. q=p1-p2, LCMS; these are projections of a 3D correlation, not q_inv."
        )
        k("correlation_projections", 9, (1844, 676), (0, 35, 1844, 655),
          "Figure 4: measured out-side-long correlation projections and fits",
          correlation_conditions)
        k("correlation_out", 9, (1844, 676), (0, 35, 638, 655),
          "Figure 4 out projection, with a separate original charge/fit legend",
          correlation_conditions + " Right-hand legend is copied from the same original figure.",
          legend_rect=(1060, 85, 1205, 210), legend_scale=1.0)

        k("radii_kt", 9, (997, 969), (0, 0, 997, 969),
          "Figure 5: measured LCMS radii versus pair transverse momentum",
          "STAR Preliminary, Au+Au sqrt(sNN)=3 GeV, 0-10%; rapidity-integrated acceptance without additional y_pair cut. "
          "Red Rout, blue Rside, black Rlong; filled pi-pi-, open pi+pi+. "
          "Error bars statistical; shaded boxes systematic. Horizontal axis is kT, not mT.")

        rapidity_conditions = (
            "STAR Preliminary, Au+Au sqrt(sNN)=3 GeV, 0.15<kT<0.6 GeV/c; black 0-10%, red 10-30%, blue 30-50%. "
            "Filled pi-pi-, open pi+pi+; boxes systematic. Gray symbols are mirrored measured points, not independent data. "
            "R_out-long^2 is a signed coefficient in -2*q_out*q_long*R_out-long^2, with fm^2 units."
        )
        k("rapidity_all", 10, (1844, 1112), (0, 0, 1844, 1112),
          "Figure 6: radii, lambda, and signed out-long coefficient versus pair CM rapidity",
          rapidity_conditions)
        k("rapidity_out_long", 10, (1844, 1112), (615, 565, 1233, 1112),
          "Figure 6e: signed out-long coefficient, with the original legend",
          rapidity_conditions,
          legend_rect=(1255, 585, 1844, 870), legend_scale=0.65)
        k("rapidity_side", 10, (1844, 1112), (615, 0, 1233, 558),
          "Figure 6b: Rside versus pair CM rapidity, with the original legend",
          rapidity_conditions,
          legend_rect=(1255, 585, 1844, 870), legend_scale=0.65)

        model_conditions = (
            "SMASH model Au+Au sqrt(sNN)=3 GeV; selected charged pions, |p| and pT in [0.15,1.5] GeV/c, "
            "eta in [-1.5,1.5]. Historical 2025 result; no new simulation or refit."
        )
        # These larger single-panel PDFs retain complete original headings and
        # legends.  The 3D raster is split at empty gutters, never across axes.
        b("correlation_kt", 13, (896, 672), (0, 0, 896, 672),
          "Archival 1D correlations for four kT intervals at centrality 0-10%",
          model_conditions + " Fixed centrality 0-10%; kT=[0.15,0.25],[0.25,0.35],[0.35,0.45],[0.45,0.60] GeV/c. Horizontal axis q_inv.")
        b("correlation_centrality", 13, (896, 672), (0, 0, 896, 672),
          "Archival 1D correlations for four centralities at fixed kT",
          model_conditions + " Fixed kT=[0.15,0.25] GeV/c; centralities 0-10%,10-30%,30-50%,50-80%. Horizontal axis q_inv.",
          image_index=1)

        b("rinv", 14, (696, 472), (0, 0, 696, 472),
          "Archival R_inv versus kT, four centralities",
          model_conditions + " Effective 1D Gaussian parameter in fm; original fit errors retained, error algorithm unknown.")
        b("lambda_1d", 14, (696, 472), (0, 0, 696, 472),
          "Archival 1D lambda versus kT, four centralities",
          model_conditions + " Dimensionless 1D fit strength; original fit errors retained, error algorithm unknown.",
          image_index=1)
        for name, rect, description in [
            ("rout_3d", (0, 0, 437, 334), "Rout"),
            ("rside_3d", (443, 0, 896, 334), "Rside"),
            ("rlong_3d", (0, 342, 437, 672), "Rlong"),
            ("lambda_3d", (443, 342, 896, 672), "3D lambda"),
        ]:
            b(name, 17, (896, 672), rect,
              f"Archival {description} versus kT, individual 3D panel",
              model_conditions + " Centralities 0-10%,10-30%,30-50%,50-80%. Original axes and legend retained; error algorithm unknown.")

    (FIGURES / "research_figure_manifest.json").write_text(
        json.dumps(assets, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    for entry in assets:
        print(f"{entry['asset']}: {entry['pdf_points']} pt; source {entry['source_raster_pixels']} px")


if __name__ == "__main__":
    main()

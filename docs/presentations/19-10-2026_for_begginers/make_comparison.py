"""Digitize published marker centres for a qualitative STAR/SMASH comparison.

This is NOT an original numerical table. All coordinates below refer to the
unscaled embedded rasters in the two cited PDFs. Pixel-coordinate calibration
and conservative manual centre-location uncertainty are exported with the data.
The latter is digitization precision, NOT a statistical or systematic error.
Original error bars and systematic boxes are intentionally not reconstructed.

Run with the presentation's Python environment: python make_comparison.py
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/cf-beginners-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pymupdf


BASE = Path(__file__).resolve().parent
FIGURES = BASE / "figures"
CSV_NAME = "comparison_digitized_2025_star.csv"
JSON_NAME = "comparison_digitized_2025_star_provenance.json"
PLOT_NAME = "comparison_star_smash_2025"

SOURCES = {
    "SMASH_Balinov_2025": {
        "path": "sources/balinov_2025.pdf",
        "sha256": "8d822f7a2618fad46bb8c5376530f7c2622e8fd17952614a3b0217679febe4a4",
        "pdf_page": 17,
        "raster_size_pixels": [896, 672],
        "url": "http://web-static-data.msu.dubna.ru:8093/epp/kursov_4/2025/Balinov.pdf",
        "description": "Balinov 2025 defence, slide 17; red 0-10% points in the three 3D-radius panels.",
        "selection": {
            "system": "Au+Au", "sqrt_sNN_GeV": 3, "centrality_percent": [0, 10],
            "charge": "not_specified_on_source_panel",
            "single_particle_p_GeV_c": [0.15, 1.5],
            "single_particle_pT_GeV_c": [0.15, 1.5],
            "single_particle_eta": [-1.5, 1.5],
            "kt_bins_GeV_c": [[0.15, 0.25], [0.25, 0.35], [0.35, 0.45], [0.45, 0.60]],
        },
        "limitations": [
            "Historical numerical results; no new simulation or refit.",
            "The archive does not identify the charge of these radius panels.",
            "Centrality calibration, original fit implementation and error method are not independently established.",
            "Formula errors in the original slides prevent independently certifying the historical radius scale.",
        ],
    },
    "STAR_Kraeva_2024": {
        "path": "sources/kraeva_2024.pdf",
        "sha256": "9190fd23d5e4fdc30dd4c92fa53e4c0a0c69fc6a808468b47c69aa68cf7999d0",
        "pdf_page": 9,
        "raster_size_pixels": [997, 969],
        "url": "https://doi.org/10.3390/universe10040188",
        "description": "A. Kraeva for STAR, Universe 10 (2024) 188, Figure 5; STAR Preliminary, CC BY 4.0.",
        "selection": {
            "system": "Au+Au", "sqrt_sNN_GeV": 3, "centrality_percent": [0, 10],
            "single_particle_eta_lab": [-2, 0],
            "single_particle_p_GeV_c": [0.15, 1.5],
            "pair_rapidity": "Acceptance integrated; no additional pair-rapidity cut.",
            "frame": "LCMS",
            "fit": "Bowler-Sinyukov; Coulomb Gaussian source radius 5 fm; nonzero out-long term.",
        },
        "limitations": [
            "Only the four visible groups of Figure 5 are digitized. The article text lists five kT bins; an absent fifth group is not inferred.",
            "The small charge-dependent horizontal displacement is preserved as plotted. These x values are graphical positions, not a new determination of bin means.",
            "Marker centres only. Statistical bars and systematic boxes remain in the original figure and are not transferred to this comparison.",
        ],
    },
}

# Tick centres in the ORIGINAL rasters, in pixels. The fit uses every listed
# labelled tick, rather than the border of the plot or an assumed raster DPI.
MODEL_X = [[100, .2], [170, .3], [241, .4], [312, .5], [383, .6]]
MODEL_SIDE_X = [[548, .2], [618, .3], [689, .4], [760, .5], [831, .6]]
CALIBRATION_TICKS = {
    "model_out": {
        "x": MODEL_X,
        "y": [[34, 7], [60, 6.5], [87, 6], [114, 5.5], [141, 5],
              [168, 4.5], [194, 4], [221, 3.5], [248, 3], [275, 2.5], [301, 2]],
    },
    "model_side": {
        "x": MODEL_SIDE_X,
        "y": [[34, 5], [78, 4.5], [123, 4], [168, 3.5], [212, 3], [257, 2.5], [301, 2]],
    },
    "model_long": {
        "x": MODEL_X,
        "y": [[370, 5], [414, 4.5], [459, 4], [504, 3.5], [548, 3], [593, 2.5], [637, 2]],
    },
    "star_all": {
        "x": [[202, .2], [314, .25], [426, .3], [539, .35], [651, .4], [763, .45], [875, .5]],
        "y": [[103, 6], [272, 5], [441, 4], [609, 3], [778, 2]],
    },
}

# Red model points: centre of the horizontal error-bar line through each marker.
# The source's other centralities and its lambda panel are not used.
MODEL_CENTRES = {
    "out": [[100, 80], [170, 102], [241, 126], [330, 137]],
    "side": [[548, 67], [618, 110], [689, 141], [778, 152]],
    "long": [[100, 395], [170, 490], [241, 548], [330, 577]],
}

# STAR: filled = pi-pi-, open = pi+pi+. Original x displacements are retained.
# Triangular markers are located at their visual geometric centre; their
# conservative 4-pixel location budget includes marker-shape ambiguity.
STAR_CENTRES = {
    "out": {
        "pi-pi-": [[202, 221], [426, 310], [651, 405], [875, 480]],
        "pi+pi+": [[224, 236], [449, 308], [673, 394], [898, 468]],
    },
    "side": {
        "pi-pi-": [[202, 235], [426, 354], [651, 457], [875, 524]],
        "pi+pi+": [[224, 311], [449, 398], [673, 490], [898, 556]],
    },
    "long": {
        "pi-pi-": [[202, 421.5], [426, 549.5], [651, 650.5], [875, 725.5]],
        "pi+pi+": [[224, 480], [449, 586], [673, 668], [898, 746]],
    },
}


def linear_calibration(ticks: list[list[float]]) -> dict:
    pixels, values = np.asarray(ticks, dtype=float).T
    slope, intercept = np.polyfit(pixels, values, 1)
    residual_pixels = (values - (slope * pixels + intercept)) / abs(slope)
    if np.max(np.abs(residual_pixels)) > .7:
        raise ValueError("Tick calibration residual exceeds 0.7 source pixels")
    return {
        "value_per_pixel": float(slope), "value_at_pixel_zero": float(intercept),
        "ticks_pixel_value": ticks,
        "max_tick_residual_pixels": float(np.max(np.abs(residual_pixels))),
    }


def extract_source_metadata() -> dict:
    result = {}
    for name, source in SOURCES.items():
        path = BASE / source["path"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != source["sha256"]:
            raise ValueError(f"Source PDF changed: {path}")
        with pymupdf.open(path) as document:
            page = document[source["pdf_page"] - 1]
            matches = [item for item in page.get_images(full=True)
                       if list(item[2:4]) == source["raster_size_pixels"]]
            if len(matches) != 1:
                raise ValueError(f"Expected one exact-size source raster in {path}")
            raster = document.extract_image(matches[0][0])
        result[name] = {**source, "embedded_raster_sha256": hashlib.sha256(raster["image"]).hexdigest()}
    return result


def digitize() -> tuple[list[dict], dict]:
    calibration = {name: {axis: linear_calibration(ticks) for axis, ticks in axes.items()}
                   for name, axes in CALIBRATION_TICKS.items()}
    rows = []

    def add(source: str, charge: str, radius: str, points: list, key: str, pixel_uncertainty: float):
        cal = calibration[key]
        for index, (x, y) in enumerate(points):
            flags = ["digitized_marker_centre", "original_errors_not_transferred", "acceptance_differs"]
            if source == "SMASH_Balinov_2025":
                flags += ["historical_result", "charge_not_specified", "historical_fit_not_independently_validated"]
                if index == 3:
                    flags.append("last_kt_bin_differs_from_STAR")
            else:
                flags += ["STAR_Preliminary", "original_graphical_x_offset_preserved"]
            values = {axis: cal[axis]["value_per_pixel"] * pixel + cal[axis]["value_at_pixel_zero"]
                      for axis, pixel in [("x", x), ("y", y)]}
            rows.append({
                "source": source, "charge": charge, "radius": radius,
                "kt": round(values["x"], 6), "R": round(values["y"], 6),
                "estimatepixeluncertainty": pixel_uncertainty,
                "flags": ";".join(flags), "pixel_x": x, "pixel_y": y,
                "kt_digitization_precision_GeV_c": round(abs(cal["x"]["value_per_pixel"]) * pixel_uncertainty, 6),
                "R_digitization_precision_fm": round(abs(cal["y"]["value_per_pixel"]) * pixel_uncertainty, 6),
                "calibration": key, "visible_group_index": index + 1,
            })

    for radius, points in MODEL_CENTRES.items():
        add("SMASH_Balinov_2025", "not_specified", radius, points, f"model_{radius}", 2.0)
    for radius, charges in STAR_CENTRES.items():
        for charge, points in charges.items():
            add("STAR_Kraeva_2024", charge, radius, points, "star_all", 4.0 if radius == "side" else 3.0)
    if len(rows) != 36:
        raise ValueError("Expected 12 model and 24 STAR marker centres")
    return rows, calibration


def make_plot(rows: list[dict]):
    ink, teal, gold, muted = "#142d43", "#007f86", "#b06a16", "#617589"
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 18,
        "axes.labelsize": 22, "axes.titlesize": 24,
        "xtick.labelsize": 19, "ytick.labelsize": 19,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": muted, "axes.labelcolor": ink,
        "text.color": ink, "xtick.color": ink, "ytick.color": ink,
        "legend.frameon": False, "pdf.fonttype": 42,
    })
    figure, axes = plt.subplots(1, 3, figsize=(12.8, 5.2), sharex=True, sharey=True)
    figure.subplots_adjust(left=.072, right=.97, top=.80, bottom=.29, wspace=.16)
    styles = [
        ("SMASH_Balinov_2025", "not_specified", dict(color=teal, marker="D", markersize=7, linewidth=2.3, linestyle="-")),
        ("STAR_Kraeva_2024", "pi-pi-", dict(color=ink, marker="o", markersize=8, linestyle="none")),
        ("STAR_Kraeva_2024", "pi+pi+", dict(color=gold, marker="o", markersize=8, markerfacecolor="none", markeredgewidth=1.8, linestyle="none")),
    ]
    for axis, radius in zip(axes, ["out", "side", "long"]):
        for source, charge, style in styles:
            selected = [row for row in rows if row["source"] == source and row["charge"] == charge and row["radius"] == radius]
            axis.plot([row["kt"] for row in selected], [row["R"] for row in selected], **style)
        axis.set(xlim=(.15, .60), ylim=(2, 8), xticks=[.2, .3, .4, .5, .6], yticks=[2, 3, 4, 5, 6, 7, 8],
                 xlabel=r"$k_T$, ГэВ/$c$", title=rf"$R_{{\rm {radius}}}$")
        axis.grid(axis="y", color=muted, alpha=.16, linewidth=.8)
        axis.tick_params(length=4, pad=5)
    axes[0].set_ylabel("Радиус, фм", labelpad=10)
    handles = [
        Line2D([], [], color=teal, marker="D", linewidth=2.3, markersize=7, label="Расчёт SMASH (2025)"),
        Line2D([], [], color=ink, marker="o", linestyle="none", markersize=8, label=r"STAR $\pi^-\pi^-$"),
        Line2D([], [], color=gold, marker="o", markerfacecolor="none", markeredgewidth=1.8,
               linestyle="none", markersize=8, label=r"STAR $\pi^+\pi^+$"),
    ]
    figure.legend(handles=handles, loc="upper center", bbox_to_anchor=(.53, 1.02), ncol=3,
                  fontsize=17, handlelength=1.5, columnspacing=1.5)
    figure.text(.5, .075, "Центры точек оцифрованы; погрешности — на исходных графиках.",
                ha="center", va="center", fontsize=14.5, color=muted)
    figure.text(.5, .027, "STAR Preliminary. Условия отбора в расчёте и эксперименте различаются.",
                ha="center", va="center", fontsize=14.5, color=muted)
    figure.savefig(FIGURES / f"{PLOT_NAME}.pdf", transparent=True, metadata={
        "Title": "Qualitative comparison: Balinov SMASH 2025 and STAR Preliminary 2024",
        "Author": "Digitized from T. Balinov (2025) and A. Kraeva for STAR (2024)",
        "Subject": "Au+Au 3 GeV, 0-10%, LCMS radii. Digitized marker centres only; original errors not transferred; differing selections. Not original numerical tables.",
    })
    figure.savefig(FIGURES / f"{PLOT_NAME}.png", dpi=180, transparent=True)
    plt.close(figure)


def main():
    FIGURES.mkdir(exist_ok=True)
    sources = extract_source_metadata()
    rows, calibration = digitize()
    with (FIGURES / CSV_NAME).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "kind": "MANUAL DIGITIZATION; NOT ORIGINAL NUMERICAL TABLES",
        "method": "Manual marker-centre pixel coordinates on original embedded rasters; linear least-squares calibration from labelled tick centres.",
        "pixel_origin": "Top-left of unscaled original embedded raster; x rightwards, y downwards.",
        "units": {"kt": "GeV/c", "R": "fm", "estimatepixeluncertainty": "source-raster pixels"},
        "precision_warning": "estimatepixeluncertainty and derived digitization_precision fields are conservative manual location budgets, not confidence intervals, statistical errors or systematic errors. CSV decimal places are computational output, not measurement accuracy.",
        "error_bars": "Only centres are digitized. Original statistical errors and systematic boxes are not transferred and must be consulted in the original figures. They are not set to zero.",
        "purpose": "Compare radius scale and visible kT trends qualitatively; no goodness-of-fit, ratio, chi-square or statistical agreement is inferred.",
        "common_conditions": "Au+Au, sqrt(sNN)=3 GeV, 0-10%, out-side-long radii in LCMS.",
        "comparison_limitations": [
            "Model and experiment have different single-particle eta and resulting pair-rapidity acceptances.",
            "Model charge is not specified on the archival source panels; STAR signs are kept separate.",
            "Model and experiment use different final kT bins and incompletely matched fit/FSI procedures.",
            "This is a comparison with the displayed 2025 model archive, not the new 2026 analysis.",
            "The historical model implementation and absolute radius scale are not independently certified.",
            "The plotted experimental points retain the STAR Preliminary status of the source.",
        ],
        "plot": {"asset": f"figures/{PLOT_NAME}.pdf", "x_limits_GeV_c": [.15, .60], "y_limits_fm": [2, 8],
                 "model_connections": "Straight segments guiding the eye through the four digitized centres; not a fit.",
                 "caption_ru": "Центры точек оцифрованы; погрешности — на исходных графиках. STAR Preliminary. Условия отбора в расчёте и эксперименте различаются.",
                 "attribution": "T. Balinov, defence 2025, slide 17; A. Kraeva for STAR, Universe 10 (2024) 188, Fig. 5, doi:10.3390/universe10040188, CC BY 4.0."},
        "sources": sources, "calibrations": calibration, "points": rows,
    }
    (FIGURES / JSON_NAME).write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    make_plot(rows)
    print(f"Wrote {len(rows)} digitized marker centres and {PLOT_NAME}.pdf/.png")
    for name, axes in calibration.items():
        print(name, "max calibration residual [pixels]", max(axis["max_tick_residual_pixels"] for axis in axes.values()))


if __name__ == "__main__":
    main()

"""Reproduce a qualitative comparison of current SMASH and STAR negative pions.

The model centres are manually digitized from the original user-supplied JPEG.
STAR centres and their original pixel calibration are reused verbatim from the
previous, separately preserved digitization. Neither table is original fit data.
No error bars are reconstructed, set to zero, or treated as negligible.
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
from PIL import Image


BASE = Path(__file__).resolve().parent
FIGURES = BASE / "figures"
CURRENT_DIR = BASE / "sources/current_2026_10_19"
CURRENT_IMAGE = CURRENT_DIR / "kt_c_all_graphs_neg.jpg"
OLD_CSV = FIGURES / "comparison_digitized_2025_star.csv"
OLD_PROVENANCE = FIGURES / "comparison_digitized_2025_star_provenance.json"
CSV_NAME = "comparison_digitized_current_2026_star.csv"
JSON_NAME = "comparison_digitized_current_2026_star_provenance.json"
PLOT_NAME = "comparison_star_smash_current_2026"
CURRENT_SOURCE = "SMASH_current_images_2026_10_19"
STAR_SOURCE = "STAR_Kraeva_2024"
IMAGE_SHA256 = "3410a3b9e125bc73449b2c3369cd58af2acfd7d57505df2eedd6fab2875dd82a"
OLD_CSV_SHA256 = "a6557a447ec4430d234b417a198cc94681a11d639a779b145993fd67a343c133"

# These are original 3196 x 1572 JPEG pixel coordinates, not coordinates on a
# resized screenshot. Tick centres were checked against the actual black lines.
X_TICKS = [[287, .2], [444, .25], [602, .3], [759, .35], [917, .4],
           [1074, .45], [1231, .5], [1389, .55]]
Y_TICKS_UPPER = [[91, 10], [151, 9], [211, 8], [272, 7], [332, 6], [393, 5],
                 [453, 4], [513, 3], [574, 2], [634, 1], [694, 0]]
Y_TICKS_LONG = [[877, 10], [937, 9], [997, 8], [1058, 7], [1118, 6], [1179, 5],
               [1239, 4], [1299, 3], [1360, 2], [1420, 1], [1480, 0]]
CURRENT_PANELS = {
    "out": {"axes_bounds_pixels": [185, 91, 1412, 694],
            "x_ticks": X_TICKS, "y_ticks": Y_TICKS_UPPER,
            "centres": [[287, 321], [602, 342], [917, 369], [1310, 407]]},
    "side": {"axes_bounds_pixels": [1783, 91, 3010, 694],
             "x_ticks": [[x + 1598, value] for x, value in X_TICKS],
             "y_ticks": Y_TICKS_UPPER,
             "centres": [[1885, 403], [2200, 438], [2515, 462], [2908, 471]]},
    "long": {"axes_bounds_pixels": [185, 877, 1412, 1480],
             "x_ticks": X_TICKS, "y_ticks": Y_TICKS_LONG,
             "centres": [[287, 1190], [602, 1258], [917, 1302], [1310, 1322]]},
}

# Source annotations read from the separate first-kT-bin CF image, used only to
# cross-check marker location. The printed fit errors are NOT the pixel budgets.
FIT_BOX = {
    "out": {"R_fm": 6.185, "printed_fit_error_fm": .006},
    "side": {"R_fm": 4.832, "printed_fit_error_fm": .004},
    "long": {"R_fm": 4.804, "printed_fit_error_fm": .005},
}
SUPPORTING_IMAGES = {
    "kt_cfs_neg_0-10_0.15-0.25.jpg": "6ebd20159b97762c39677d7065f763fe329df4a48ea2039607782c85b3ca6809",
    "c_all_graphs_neg.jpg": "16051d82e8659f64125ee7eebaba069204a744503dbb14c4d66620d4b96d90a8",
    "cfs_neg_0-10_0.0_0.2.jpg": "395b91a05da0e1ba1e9f8db91e53721c3d7b20cb0e385facaa96b302246f3acd",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def calibrate(ticks: list[list[float]]) -> dict:
    pixels, values = np.asarray(ticks, dtype=float).T
    slope, intercept = np.polyfit(pixels, values, 1)
    residual = (values - slope * pixels - intercept) / abs(slope)
    maximum = float(np.max(np.abs(residual)))
    if maximum > .7:
        raise ValueError("Calibration residual exceeds 0.7 original-image pixels")
    return {"value_per_pixel": float(slope), "value_at_pixel_zero": float(intercept),
            "ticks_pixel_value": ticks, "max_tick_residual_pixels": maximum}


def current_rows() -> tuple[list[dict], dict, list[dict]]:
    if sha256(CURRENT_IMAGE) != IMAGE_SHA256:
        raise ValueError("Current source JPEG changed; recheck all manual pixel coordinates")
    if Image.open(CURRENT_IMAGE).size != (3196, 1572):
        raise ValueError("Unexpected current JPEG size")
    rows, calibrations, checks = [], {}, []
    pixel_budget = 2.0
    for radius, panel in CURRENT_PANELS.items():
        calibration = {axis: calibrate(panel[f"{axis}_ticks"]) for axis in ["x", "y"]}
        key = f"current_{radius}"
        calibrations[key] = {"axes_bounds_pixels": panel["axes_bounds_pixels"], **calibration}
        for index, (x, y) in enumerate(panel["centres"]):
            kt = calibration["x"]["value_per_pixel"] * x + calibration["x"]["value_at_pixel_zero"]
            radius_value = calibration["y"]["value_per_pixel"] * y + calibration["y"]["value_at_pixel_zero"]
            radius_budget = abs(calibration["y"]["value_per_pixel"]) * pixel_budget
            flags = ["digitized_marker_centre", "JPEG_location_budget_not_physical_error",
                     "original_fit_errors_not_reconstructed", "selections_may_differ"]
            if index == 3:
                flags.append("last_kt_group_position_differs_from_STAR")
            rows.append({
                "source": CURRENT_SOURCE, "charge": "pi-pi-", "radius": radius,
                "kt": round(float(kt), 6), "R": round(float(radius_value), 6),
                "estimatepixeluncertainty": pixel_budget, "flags": ";".join(flags),
                "pixel_x": x, "pixel_y": y,
                "kt_digitization_precision_GeV_c": round(abs(calibration["x"]["value_per_pixel"]) * pixel_budget, 6),
                "R_digitization_precision_fm": round(radius_budget, 6),
                "calibration": key, "visible_group_index": index + 1,
                "series_colour_in_source": "red", "source_digitization_flags": "",
            })
            if index == 0:
                difference = radius_value - FIT_BOX[radius]["R_fm"]
                checks.append({
                    "radius": radius,
                    "source_image": "sources/current_2026_10_19/kt_cfs_neg_0-10_0.15-0.25.jpg",
                    **FIT_BOX[radius], "digitized_R_fm": round(radius_value, 6),
                    "difference_fm": round(difference, 6),
                    "within_digitization_location_budget": abs(difference) <= radius_budget,
                })
                if abs(difference) > radius_budget:
                    raise ValueError(f"First-point fit-box cross-check failed for {radius}")
    return rows, calibrations, checks


def star_rows() -> tuple[list[dict], dict, dict]:
    if sha256(OLD_CSV) != OLD_CSV_SHA256:
        raise ValueError("Previous digitization CSV changed; explicitly review the STAR import")
    original = json.loads(OLD_PROVENANCE.read_text(encoding="utf-8"))
    with OLD_CSV.open(encoding="utf-8", newline="") as handle:
        imported = [row for row in csv.DictReader(handle)
                    if row["source"] == STAR_SOURCE and row["charge"] == "pi-pi-"]
    if len(imported) != 12:
        raise ValueError("Expected exactly 12 STAR negative-pion centres")
    rows = []
    for row in imported:
        converted = row.copy()
        for name in ["kt", "R", "estimatepixeluncertainty", "pixel_x", "pixel_y",
                     "kt_digitization_precision_GeV_c", "R_digitization_precision_fm"]:
            converted[name] = float(row[name])
        converted["visible_group_index"] = int(row["visible_group_index"])
        converted["source_digitization_flags"] = row["flags"]
        # The old acceptance_differs flag concerned the 2025 archive. For the
        # current batch the selections are not yet matched, so only "may differ"
        # is claimed. All original flags remain in source_digitization_flags.
        converted["flags"] = "digitized_marker_centre;original_errors_not_transferred;STAR_Preliminary;selections_may_differ;imported_STAR_coordinates_unchanged"
        converted["series_colour_in_source"] = {"out": "red", "side": "blue", "long": "black"}[row["radius"]]
        rows.append(converted)
    source = original["sources"][STAR_SOURCE]
    provenance = {
        "source_pdf": source["path"], "source_pdf_sha256": source["sha256"],
        "source_pdf_page": source["pdf_page"], "source_figure": 5,
        "source_url": source["url"], "source_raster_pixels": source["raster_size_pixels"],
        "embedded_raster_sha256": source["embedded_raster_sha256"],
        "credit": source["description"],
        "import_csv": str(OLD_CSV.relative_to(BASE)), "import_csv_sha256": sha256(OLD_CSV),
        "import_provenance": str(OLD_PROVENANCE.relative_to(BASE)),
        "import_provenance_sha256": sha256(OLD_PROVENANCE),
        "original_imported_csv_records": imported,
        "note": "STAR pixel coordinates, calibrated values and digitization budgets are unchanged. Only the negative-pion series is selected; no new STAR digitization is performed.",
    }
    return rows, original["calibrations"]["star_all"], provenance


def source_provenance() -> dict:
    supporting = []
    for name, expected_digest in SUPPORTING_IMAGES.items():
        path = CURRENT_DIR / name
        if sha256(path) != expected_digest:
            raise ValueError(f"Supporting source JPEG changed: {name}")
        supporting.append({"path": str(path.relative_to(BASE)), "sha256": expected_digest,
                           "image_pixels": list(Image.open(path).size)})
    return {
        "path": str(CURRENT_IMAGE.relative_to(BASE)), "sha256": IMAGE_SHA256,
        "image_pixels": [3196, 1572], "origin": "User-provided current calculation images, 19 October 2026.",
        "selected_series": "Red points; negative pions; centrality 0-10% confirmed by the user.",
        "confirmed_by_user": {"model": "SMASH", "system": "Au+Au", "sqrt_sNN_GeV": 3,
                              "centrality_percent": [0, 10], "radii_frame": "LCMS",
                              "rapidity_frame": "collision_center_of_mass"},
        "readable_in_fit_header": {"charge": "negative", "first_kt_interval_GeV_c": [.15, .25]},
        "supporting_images": supporting,
        "unconfirmed_current_conditions": [
            "Single-particle cuts and pair-rapidity integration for the kT analysis.",
            "Fit/FSI settings, centrality calibration and interpretation of printed fit errors.",
        ],
    }


def plot(rows: list[dict]):
    ink, teal, muted = "#142d43", "#007f86", "#617589"
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 18,
        "axes.labelsize": 22, "axes.titlesize": 24, "xtick.labelsize": 19, "ytick.labelsize": 19,
        "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": muted,
        "axes.labelcolor": ink, "text.color": ink, "xtick.color": ink, "ytick.color": ink,
        "legend.frameon": False, "pdf.fonttype": 42,
    })
    figure, axes = plt.subplots(1, 3, figsize=(12.8, 5.2), sharex=True, sharey=True)
    figure.subplots_adjust(left=.072, right=.97, top=.80, bottom=.29, wspace=.16)
    for axis, radius in zip(axes, ["out", "side", "long"]):
        model = [row for row in rows if row["source"] == CURRENT_SOURCE and row["radius"] == radius]
        star = [row for row in rows if row["source"] == STAR_SOURCE and row["radius"] == radius]
        axis.plot([row["kt"] for row in model], [row["R"] for row in model],
                  color=teal, marker="D", markersize=7.5, linewidth=2.3)
        axis.plot([row["kt"] for row in star], [row["R"] for row in star],
                  color=ink, marker="o", markersize=8.5, linestyle="none")
        axis.set(xlim=(.15, .60), ylim=(2, 8), xticks=[.2, .3, .4, .5, .6], yticks=[2, 3, 4, 5, 6, 7, 8],
                 xlabel=r"$k_T$, ГэВ/$c$", title=rf"$R_{{\rm {radius}}}$")
        axis.grid(axis="y", color=muted, alpha=.16, linewidth=.8)
        axis.tick_params(length=4, pad=5)
    axes[0].set_ylabel("Радиус, фм", labelpad=10)
    handles = [
        Line2D([], [], color=teal, marker="D", linewidth=2.3, markersize=7.5,
               label=r"Мой текущий расчёт SMASH: $\pi^-\pi^-$"),
        Line2D([], [], color=ink, marker="o", linestyle="none", markersize=8.5,
               label=r"STAR: $\pi^-\pi^-$"),
    ]
    figure.legend(handles=handles, loc="upper center", bbox_to_anchor=(.53, 1.02), ncol=2,
                  fontsize=20, handlelength=1.5, columnspacing=2.0)
    figure.text(.5, .075, "Оцифрованы центры точек; погрешности не восстановлены.",
                ha="center", va="center", fontsize=14.5, color=muted)
    figure.text(.5, .027, "STAR Preliminary. Отборы могут различаться.",
                ha="center", va="center", fontsize=14.5, color=muted)
    figure.savefig(FIGURES / f"{PLOT_NAME}.pdf", transparent=True, metadata={
        "Title": "Current SMASH calculation and STAR Preliminary: digitized negative-pion radii",
        "Author": "T. Balinov; STAR points from A. Kraeva for STAR (2024)",
        "Subject": "Au+Au 3 GeV, 0-10%, LCMS radii. Qualitative comparison of digitized centres only; uncertainties not reconstructed or assumed zero; selection matching unverified.",
    })
    figure.savefig(FIGURES / f"{PLOT_NAME}.png", dpi=180, transparent=True)
    plt.close(figure)


def main():
    model, model_calibrations, checks = current_rows()
    star, star_calibration, star_provenance = star_rows()
    rows = model + star
    if len(rows) != 24:
        raise ValueError("Expected 12 current-model and 12 STAR negative-pion centres")
    provenance = {
        "kind": "DIGITIZED MARKER CENTRES; NOT ORIGINAL NUMERICAL TABLES",
        "purpose": "Qualitative comparison of visible radius scales and kT trends; no chi-square or statistical agreement is inferred.",
        "common_confirmed_conditions": {"system": "Au+Au", "sqrt_sNN_GeV": 3,
                                         "centrality_percent": [0, 10], "pion_pair_charge": "negative",
                                         "radii_frame": "LCMS"},
        "method": "Current red-series centres located manually on the unscaled original JPEG; least-squares linear tick calibration. STAR negative-pion centres imported without recalibration or coordinate changes.",
        "pixel_origin": "Top-left of each original source raster; x rightwards, y downwards.",
        "units": {"kt": "GeV/c", "R": "fm", "estimatepixeluncertainty": "source-raster pixels"},
        "precision_warning": "Pixel location budgets and derived digitization_precision columns are not statistical or systematic uncertainties, confidence intervals, or fit errors. CSV decimal places express computation only.",
        "physical_errors": "Original physical uncertainties are not reconstructed and are not assumed zero. They are not represented by error bars in this comparison; consult source figures and original fit results.",
        "comparison_limitations": [
            "This compares the newly provided current model images with STAR; it is not the 2025 model archive.",
            "Detailed current track/pair selections and fit treatment have not yet been matched to STAR.",
            "The last plotted kT group differs between model and experiment; no rebinning or interpolation onto common bin centres is performed.",
            "Model line segments only guide the eye through the four centres; they are not an additional fit.",
            "STAR Preliminary status from the source is retained.",
        ],
        "sources": {CURRENT_SOURCE: source_provenance(), STAR_SOURCE: star_provenance},
        "calibrations": {**model_calibrations, "star_all": star_calibration},
        "first_point_fit_box_checks": checks,
        "points": rows,
        "assets": {"csv": f"figures/{CSV_NAME}", "pdf": f"figures/{PLOT_NAME}.pdf", "png": f"figures/{PLOT_NAME}.png"},
    }
    with (FIGURES / CSV_NAME).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (FIGURES / JSON_NAME).write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    plot(rows)
    print(f"Wrote {len(rows)} centres, current comparison CSV/JSON and PDF/PNG")
    for check in checks:
        print(check["radius"], "first point minus printed fit value [fm]", check["difference_fm"])


if __name__ == "__main__":
    main()

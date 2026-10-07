"""Bounded-memory plots for the CSV produced by statistics_diff."""
import collections
import math

from summarize_statistics import rows, number, key, usable


def histogram(counter, value, minimum=-6., maximum=6., bins=120):
    if not math.isfinite(value):
        counter["invalid"] += 1
    elif value < minimum:
        counter["below"] += 1
    elif value >= maximum:
        counter["above"] += 1
    else:
        counter[int((value-minimum)/(maximum-minimum)*bins)] += 1


def plot_directory(directory):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise SystemExit("Для графиков установите matplotlib: python3 -m pip install matplotlib") from error
    output = directory / "plots"
    output.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "figure.dpi": 120})

    def save(figure, name):
        figure.tight_layout()
        figure.savefig(output / (name + ".png"), dpi=180)
        figure.savefig(output / (name + ".pdf"))
        plt.close(figure)

    fits = {(key(row), row["variant"]): row for row in rows(directory / "fits.csv")}
    groups = sorted({group for group, _ in fits}, key=lambda group: tuple(map(int, group)))
    variants = ("program", "B", "independent")
    colors = ("black", "tab:orange", "tab:blue")
    positions = {group: index for index, group in enumerate(groups)}
    figure, axes = plt.subplots(2, 2, figsize=(12, 8))
    for axis, field in zip(axes.flat, ("r_out", "r_side", "r_long", "lambda")):
        for variant, color in zip(variants, colors):
            selected = [(group, row) for (group, name), row in fits.items()
                        if name == variant and usable(row)]
            selected.sort(key=lambda item: positions[item[0]])
            axis.errorbar([positions[group] for group, _ in selected],
                          [number(row, field) for _, row in selected],
                          yerr=[number(row, "e_"+field) if math.isfinite(number(row, "e_"+field)) else 0
                                for _, row in selected], fmt=".", capsize=2,
                          color=color, label=variant, alpha=.8)
        axis.set_title(field + (" [fm]" if field.startswith("r_") else ""))
        axis.set_xlabel("Fit index (mapping: fit_index.csv)")
        axis.grid(alpha=.2)
    axes.flat[0].legend()
    figure.suptitle("Parameters and fit errors; only converged fits away from limits")
    save(figure, "parameters")
    import csv
    with (output / "fit_index.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["index", "charge", "centrality", "bin"])
        writer.writerows((index, *group) for group, index in positions.items())
    figure, axes = plt.subplots(2, 2, figsize=(12, 8))
    for axis, field in zip(axes.flat, ("r_out", "r_side", "r_long", "lambda")):
        for variant, color in zip(variants[1:], colors[1:]):
            x, shifts = [], []
            for group in groups:
                base, other = fits.get((group, "program")), fits.get((group, variant))
                if not base or not other or not usable(base) or not usable(other):
                    continue
                error = number(base, "e_" + field)
                if error > 0 and math.isfinite(error):
                    x.append(positions[group])
                    shifts.append((number(other, field)-number(base, field))/error)
            axis.plot(x, shifts, ".", label=variant, color=color)
        axis.axhline(0, color="gray", linewidth=1)
        axis.set_title(field)
        axis.set_ylabel("(alternative - program) / program error")
        axis.set_xlabel("Fit index")
        axis.grid(alpha=.2)
    axes.flat[0].legend()
    figure.suptitle("Sensitivity to error model; these are NOT independent-fit significances")
    save(figure, "parameter_shifts")
    figure, axes = plt.subplots(1, 3, figsize=(15, 4))
    for variant, color in zip(variants, colors):
        selected = [(group, row) for (group, name), row in fits.items() if name == variant and usable(row)]
        axes[0].scatter([positions[group] for group, _ in selected],
                        [number(row, "chi2_ndf") for _, row in selected], label=variant, color=color, s=12)
    axes[0].set_yscale("symlog", linthresh=1e-5)
    axes[0].set_ylabel("chi2 / ndf")
    axes[0].set_xlabel("Fit index")
    axes[0].legend()
    categories = collections.Counter((variant, "usable" if usable(row) else "failed / limit")
                                      for (_, variant), row in fits.items())
    for index, label in enumerate(("usable", "failed / limit")):
        axes[1].bar([i + index*.35 for i in range(3)],
                    [categories[variant, label] for variant in variants], width=.35, label=label)
    axes[1].set_xticks([i+.175 for i in range(3)], variants)
    axes[1].set_ylabel("Fit count")
    axes[1].legend()
    checks = list(rows(directory / "chi2_check.csv"))
    for variant, color in zip(variants, colors):
        selected = [row for row in checks if row["variant"] == variant]
        axes[2].scatter([positions[key(row)] for row in selected],
                        [number(row, "used_cells") for row in selected],
                        label=variant, color=color, s=12)
    axes[2].set_xlabel("Fit index")
    axes[2].set_ylabel("Cells used in chi2")
    axes[2].legend()
    figure.suptitle("Different error models can select different cells: chi2 alone is not a model ranking")
    save(figure, "fit_quality")
    b_ratios = collections.Counter()
    b_zero = compared = 0
    for row in rows(directory / "cells.csv"):
        program, b = number(row, "sigma_program"), number(row, "sigma_B")
        if program > 0 and math.isfinite(b):
            compared += 1
            if b == 0:
                b_zero += 1
            elif b > 0:
                histogram(b_ratios, math.log10(b/program))
    figure, axis = plt.subplots(figsize=(9, 4))
    axis.bar([-6+(i+.5)*.1 for i in range(120)], [b_ratios[i] for i in range(120)], width=.1)
    axis.set_yscale("symlog", linthresh=1)
    axis.set_xlabel("log10(sigma_B / sigma_program); all occupied cells, program error > 0")
    axis.set_ylabel("Cell count")
    axis.set_title(f"Compared: {compared}; B=0: {b_zero}; outside plotted range: "
                   f"{b_ratios['below']+b_ratios['above']}")
    save(figure, "B_errors")
    rounded = collections.defaultdict(collections.Counter)
    for row in rows(directory / "rounding.csv"):
        stats = rounded[row["digits"]]
        original, changed = number(row, "sigma_raw_original"), number(row, "sigma_raw_rounded")
        stats["negative"] += number(row, "residual_rounded") < 0
        if original > 0 and math.isfinite(original) and math.isfinite(changed):
            delta = abs(changed/original-1)
            if delta > 0:
                histogram(stats, math.log10(delta), minimum=-12., maximum=4., bins=80)
            else:
                stats["unchanged"] += 1
    figure, axis = plt.subplots(figsize=(9, 4))
    for digits, stats in sorted(rounded.items(), key=lambda item: int(item[0])):
        axis.step([-12+(i+.5)*.2 for i in range(80)], [stats[i] for i in range(80)],
                  where="mid", label=f"{digits} digits; negative residual: {stats['negative']}; "
                  f"unchanged: {stats['unchanged']}; outside: {stats['below']+stats['above']}")
    axis.axvline(-2, color="gray", linestyle="--", label="1% change")
    axis.set_yscale("symlog", linthresh=1)
    axis.set_xlabel("log10(abs(sigma_rounded / sigma_original - 1)); raw pair-weight formula")
    axis.set_ylabel("Cell count")
    axis.set_title("Stored S rounded, N and Q fixed; all occupied cells")
    axis.legend(fontsize=7)
    save(figure, "rounding")
    # Only a bounded number of representative projections: lowest indices
    # with actual output. Other projections are still included in summary.txt.
    projection_rows = collections.defaultdict(list)
    chosen = []
    for row in rows(directory / "fit_over_cf.csv"):
        group = key(row)
        if group not in chosen:
            if len(chosen) >= 3:
                continue
            chosen.append(group)
        projection_rows[group, row["axis"]].append(row)
    for group in chosen:
        figure, axes = plt.subplots(2, 3, figsize=(12, 7))
        for index, axis_name in enumerate(("out", "side", "long")):
            selected = sorted(projection_rows[group, axis_name], key=lambda row: number(row, "q"))
            q = [number(row, "q") for row in selected]
            for field, label in (("sigma_old", "old"), ("sigma_new", "new")):
                axes[0, index].plot(q, [number(row, field) for row in selected], label=label)
            axes[1, index].plot(q, [number(row, "delta_sigma") for row in selected])
            axes[0, index].set_title(axis_name)
            axes[0, index].legend()
            axes[1, index].set_xlabel("q [GeV/c]")
            axes[1, index].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        axes[0, 0].set_ylabel("Error of fit/CF")
        axes[1, 0].set_ylabel("new error - old error")
        figure.suptitle(f"Fixed fitted curve; charge, centrality, bin = {group}")
        save(figure, "fit_over_cf_" + "_".join(group))
    print(f"Графики: {output}")

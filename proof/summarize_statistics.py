#!/usr/bin/env python3
"""Stream statistics_diff CSV files using only Python's standard library."""
import argparse
import collections
import csv
import heapq
import json
import math
from pathlib import Path


def rows(path):
    if not path.is_file():
        raise FileNotFoundError(f"Отсутствует {path}")
    with path.open(newline="") as source:
        yield from csv.DictReader(source)


def number(row, field):
    try:
        return float(row[field])
    except (KeyError, ValueError, TypeError):
        return math.nan


def key(row):
    return tuple(row[field] for field in ("charge", "centrality", "bin"))


def usable(row):
    return (number(row, "status") == 0 and number(row, "cov_status") in (1, 3)
            and number(row, "at_limit") == 0 and number(row, "ndf") > 0
            and all(math.isfinite(number(row, field)) for field in
                    ("chi2", "r_out", "r_side", "r_long", "lambda")))


def largest(items, score, detail):
    heapq.heappush(items, (score, detail))
    if len(items) > 5:
        heapq.heappop(items)


def contributions(directory, fits):
    """Keep a few outliers per fit, distinguishing optimizer status."""
    totals = collections.defaultdict(lambda: [0, 0., []])
    by_fit = collections.defaultdict(lambda: [0, 0., []])
    for row in rows(directory / "chi2_cells.csv"):
        value = number(row, "chi2_contribution")
        group, variant = key(row), row["variant"]
        state = "usable" if usable(fits.get((group, variant), {})) else "failed_or_limit"
        if not math.isfinite(value):
            raise ValueError(f"Нефинитный вклад χ²: {group}, {variant}, {row['cell']}")
        stats = totals[variant, state]
        stats[0] += 1
        stats[1] += value
        largest(stats[2], value, f"{group} cell={row['cell']}: вклад={value:.6g}")
        fit_stats = by_fit[group, variant]
        fit_stats[0] += 1
        fit_stats[1] += value
        # Integer cell IDs make ties deterministic and do not store whole rows.
        largest(fit_stats[2], value, int(row["cell"]))
    wanted = collections.defaultdict(dict)
    for (group, variant), (_, total, top) in by_fit.items():
        for value, cell in top:
            wanted[group + (str(cell),)][variant] = (value, value/total if total > 0 else 0)
    lines = ["Слагаемые χ² отдельно по статусу фита (usable означает сходимость, не качество модели):"]
    for (variant, state), (count, total, top) in sorted(totals.items()):
        lines.append(f"  {variant}/{state}: ячеек {count}, χ²={total:.6g}")
        lines.extend("    " + detail for _, detail in sorted(top, reverse=True)[:3])
    for variant in ("program", "B", "independent"):
        candidates = [(number(row, "chi2_ndf"), group) for (group, name), row in fits.items()
                      if name == variant and usable(row)]
        lines.append(f"Наибольшие χ²/ndf сходящихся {variant}:")
        for ratio, group in sorted(candidates, reverse=True)[:5]:
            count, total, top = by_fit.get((group, variant), (0, 0., []))
            maximum = max((value for value, _ in top), default=0.)
            top_sum = sum(value for value, _ in top)
            lines.append(f"  {group}: χ²/ndf={ratio:.6g}, ячеек={count}, "
                         f"max вклад/χ²={maximum/total if total > 0 else 0:.2%}, "
                         f"5 крупнейших/χ²={top_sum/total if total > 0 else 0:.2%}")
    return lines, wanted


def summarize(directory):
    lines = [f"\n=== {directory.name} ==="]
    notes = directory / "README.txt"
    pair_weights = True
    if notes.exists():
        with notes.open() as source:
            for line in source:
                if line.startswith(("Input:", "ROOT:", "Statistics:")):
                    lines.append(line.rstrip())
                if line.startswith("Statistics: fixed_reference"):
                    pair_weights = False
    if not pair_weights:
        lines.append("fixed_reference: сырая дисперсия среднего весов и её округление не являются ошибкой этой модели.")
    errors = collections.Counter(row["stage"] for row in rows(directory / "errors.csv"))
    lines.append(f"Ошибки проверки: {dict(errors)}")
    fits = { (key(row), row["variant"]): row for row in rows(directory / "fits.csv") }
    contribution_lines, wanted = contributions(directory, fits)
    for variant in ("program", "B", "independent"):
        selected = [row for (_, name), row in fits.items() if name == variant]
        good = [row for row in selected if usable(row)]
        statuses = collections.Counter((row["status"], row["cov_status"], row["at_limit"])
                                       for row in selected)
        lines.append(f"Фиты {variant}: всего {len(selected)}, допустимых {len(good)}; "
                     f"(status,cov,limit): {dict(statuses)}")
        values = [number(row, "chi2_ndf") for row in good]
        if values:
            lines.append(f"  χ²/ndf допустимых: min={min(values):.6g}, max={max(values):.6g}")
        bad = [(key(row), row["status"], row["at_limit"]) for row in selected if not usable(row)]
        if bad:
            lines.append(f"  Проблемные (charge,centrality,bin),status,limit: {bad}")
    for variant in ("B", "independent"):
        shifts = []
        by_parameter = collections.defaultdict(list)
        comparable = 0
        for (group, name), base in fits.items():
            other = fits.get((group, variant))
            if name != "program" or not other or not usable(base) or not usable(other):
                continue
            comparable += 1
            for field in ("r_out", "r_side", "r_long", "lambda"):
                original, changed = number(base, field), number(other, field)
                delta = changed - original
                error = number(base, "e_" + field)
                scale = abs(delta) / error if error > 0 else math.nan
                percent = 100 * delta / original if original != 0 else math.nan
                detail = (f"{group} {field}: program={original:.6g}, {variant}={changed:.6g}, "
                          f"Δ={delta:.6g} ({percent:.3g}%), |Δ|/σ_program={scale:.3g}")
                largest(shifts, scale if math.isfinite(scale) else abs(percent), detail)
                score = scale if math.isfinite(scale) else abs(percent)
                if math.isfinite(score):
                    largest(by_parameter[field], score, detail)
        lines.append(f"Сравнение program/{variant}: {comparable} пар допустимых фитов; наибольшие сдвиги:")
        lines.extend("  " + detail for _, detail in sorted(shifts, reverse=True))
        for field, top in by_parameter.items():
            lines.append(f"  Отдельно {field}:")
            lines.extend("    " + detail for _, detail in sorted(top, reverse=True)[:2])
    lines.append("|Δ|/σ_program — масштаб чувствительности, не значимость: фиты используют одни данные.")
    checks = 0
    mismatches = 0
    maximum = 0
    max_relative = max_usable_absolute = 0
    check_outliers = []
    for row in rows(directory / "chi2_check.csv"):
        checks += 1
        difference = abs(number(row, "difference"))
        maximum = max(maximum, difference) if math.isfinite(difference) else math.inf
        tolerance = 1e-8 * max(1, abs(number(row, "program_chi2")))
        relative = difference/max(1., abs(number(row, "program_chi2")))
        max_relative = max(max_relative, relative)
        if usable(fits.get((key(row), row["variant"]), {})):
            max_usable_absolute = max(max_usable_absolute, difference)
        largest(check_outliers, relative, f"{key(row)} {row['variant']}: "
                f"χ²={number(row, 'program_chi2'):.6g}, |Δ|={difference:.6g}, "
                f"|Δ|/max(1,χ²)={relative:.3g}")
        if (not math.isfinite(difference) or difference > tolerance
                or row["manual_ndf"] != row["program_ndf"]):
            mismatches += 1
    lines.append(f"Ручной χ²/ndf: проверок {checks}, несовпадений {mismatches}, max |Δχ²|={maximum:.6g}.")
    lines.append(f"  max |Δ|/max(1,χ²)={max_relative:.6g}; max |Δχ²| сходящихся={max_usable_absolute:.6g}")
    lines.extend("  " + detail for _, detail in sorted(check_outliers, reverse=True)[:2])
    count = changed = zero_old = 0
    max_absolute = max_relative = 0
    for row in rows(directory / "fit_over_cf.csv"):
        old, new = number(row, "sigma_old"), number(row, "sigma_new")
        if not (math.isfinite(old) and math.isfinite(new)):
            continue
        count += 1
        delta = abs(new-old)
        changed += delta > 1e-10 * max(abs(new), abs(old), 1e-300)
        max_absolute = max(max_absolute, delta)
        if old > 0:
            max_relative = max(max_relative, delta/old)
        elif new > 0:
            zero_old += 1
    lines.append(f"fit/CF: {count} бинов, изменённых {changed}; max |Δσ|={max_absolute:.6g}, "
                 f"max |Δσ|/σ_old={max_relative:.6g}; σ_old=0→σ_new>0: {zero_old}.")
    moments = collections.Counter()
    fit_moments = collections.Counter()
    config_file = directory / "config.json"
    if config_file.is_file():
        with config_file.open() as source:
            q_max = json.load(source).get("fit", {}).get("q_max", .20)
    else:
        raise FileNotFoundError(f"Нет {config_file}; область фита не определяется")
    ratios = []
    diagnostic_rows = []
    for row in rows(directory / "cells.csv"):
        moments["occupied"] += 1
        raw, program, b = (number(row, field) for field in
                           ("sigma_raw", "sigma_program", "sigma_B"))
        moments["negative_residual"] += number(row, "residual") < 0
        moments["sigma_program_zero"] += program == 0
        moments["sigma_program_nan"] += not math.isfinite(program)
        moments["B_zero_program_positive"] += b == 0 and program > 0
        moments["raw_positive_program_zero"] += raw > 0 and program == 0
        inside = all(abs(number(row, field)) <= q_max for field in ("q_out", "q_side", "q_long"))
        if inside:
            fit_moments["occupied"] += 1
            fit_moments["N_le_1"] += number(row, "N") <= 1
            fit_moments["negative_residual"] += number(row, "residual") < 0
            fit_moments["program_zero"] += program == 0
            fit_moments["B_positive_program_zero"] += b > 0 and program == 0
            fit_moments["raw_positive_program_zero"] += raw > 0 and program == 0
        candidate = wanted.get(key(row) + (row["cell"],))
        if candidate:
            for variant, (value, share) in candidate.items():
                record = dict(row)
                record.update(variant=variant, chi2_contribution=value, fraction_of_fit_chi2=share,
                              usable_fit=int(usable(fits.get((key(row), variant), {}))),
                              inside_fit=int(inside))
                diagnostic_rows.append(record)
        if program > 0 and math.isfinite(b):
            ratio = b/program
            largest(ratios, abs(ratio-1), f"{key(row)} cell={row['cell']}: σ_B/σ_program={ratio:.6g}")
    lines.append(f"Моменты по всем занятым ячейкам (включая вне области фита): {dict(moments)}")
    lines.append(f"Моменты только в области фита |q_i|≤{q_max}: {dict(fit_moments)}")
    if diagnostic_rows:
        path = directory / "diagnostic_cells.csv"
        with path.open("w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=list(diagnostic_rows[0]))
            writer.writeheader()
            writer.writerows(diagnostic_rows)
        lines.append(f"Исходные моменты до 5 ведущих ячеек каждого фита: {path.name} ({len(diagnostic_rows)} строк)")
        for variant in ("program", "B", "independent"):
            selected = sorted((row for row in diagnostic_rows if row["variant"] == variant),
                              key=lambda row: row["chi2_contribution"], reverse=True)[:3]
            for row in selected:
                lines.append(f"  {variant} {key(row)} cell={row['cell']}, usable={row['usable_fit']}: "
                             f"N={number(row, 'N'):.6g}, C={number(row, 'C'):.6g}, "
                             f"Q={number(row, 'Q'):.6g}, residual={number(row, 'residual'):.6g}, "
                             f"σ_program={number(row, 'sigma_program'):.6g}, "
                             f"σ_B={number(row, 'sigma_B'):.6g}, доля χ² фита={row['fraction_of_fit_chi2']:.2%}")
    lines.extend("  " + detail for _, detail in sorted(ratios, reverse=True))
    rounding = collections.defaultdict(collections.Counter)
    extremes = collections.defaultdict(list)
    for row in rows(directory / "rounding.csv"):
        digits = row["digits"]
        stats = rounding[digits]
        stats["total"] += 1
        original, rounded = number(row, "sigma_raw_original"), number(row, "sigma_raw_rounded")
        stats["negative_residual_after"] += number(row, "residual_rounded") < 0
        stats["positive_to_zero"] += original > 0 and rounded == 0
        stats["zero_to_positive"] += original == 0 and rounded > 0
        stats["positive_to_unavailable"] += original > 0 and not math.isfinite(rounded)
        if original > 0 and math.isfinite(original) and math.isfinite(rounded):
            stats["comparable"] += 1
            relative = abs(rounded/original-1)
            stats["over_1pct"] += relative > .01
            stats["over_10pct"] += relative > .1
            largest(extremes[digits], relative,
                    f"{key(row)} cell={row['cell']}: |Δσ|/σ={relative:.6g}, ΔC={number(row, 'delta_C'):.6g}")
    lines.append("Округление S при фиксированных N,Q; сырые ошибки по всем занятым ячейкам:")
    for digits in sorted(rounding, key=int):
        lines.append(f"  {digits} знаков: {dict(rounding[digits])}")
        stats = rounding[digits]
        denominator = stats["comparable"]
        if denominator:
            lines.append(f"    Среди сравнимых: >1% у {stats['over_1pct']/denominator:.4%}, "
                         f">10% у {stats['over_10pct']/denominator:.4%}")
        lines.extend("    " + detail for _, detail in sorted(extremes[digits], reverse=True)[:2])
    lines.extend(contribution_lines)
    if not checks or not count:
        lines.append("Проверьте полноту расчёта: некоторые таблицы могут быть пустыми.")
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plots", action="store_true", help="Сохранить PNG/PDF; нужен matplotlib")
    args = parser.parse_args()
    directories = ([args.results] if (args.results / "fits.csv").exists() else
                   sorted(path for path in args.results.iterdir() if path.is_dir()))
    if not directories:
        parser.error("Нет папок результатов")
    report = ["Сводка statistics_diff", "Статусы проверяются до интерпретации параметров."]
    for directory in directories:
        try:
            report.extend(summarize(directory))
        except (OSError, KeyError, ValueError) as error:
            report.append(f"НЕПОЛНАЯ СВОДКА {directory}: {error}")
    output = args.output or args.results / "summary.txt"
    output.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Сводка: {output}")
    if args.plots:
        from plot_statistics import plot_directory
        for directory in directories:
            plot_directory(directory)


if __name__ == "__main__":
    main()

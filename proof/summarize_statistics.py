#!/usr/bin/env python3
"""Stream statistics_diff CSV files using only Python's standard library."""
import argparse
import collections
import csv
import heapq
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
    for variant in ("B", "independent"):
        shifts = []
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
        lines.append(f"Сравнение program/{variant}: {comparable} пар допустимых фитов; наибольшие сдвиги:")
        lines.extend("  " + detail for _, detail in sorted(shifts, reverse=True))
    lines.append("|Δ|/σ_program — масштаб чувствительности, не значимость: фиты используют одни данные.")
    checks = 0
    mismatches = 0
    maximum = 0
    for row in rows(directory / "chi2_check.csv"):
        checks += 1
        difference = abs(number(row, "difference"))
        maximum = max(maximum, difference) if math.isfinite(difference) else math.inf
        tolerance = 1e-8 * max(1, abs(number(row, "program_chi2")))
        if (not math.isfinite(difference) or difference > tolerance
                or row["manual_ndf"] != row["program_ndf"]):
            mismatches += 1
    lines.append(f"Ручной χ²/ndf: проверок {checks}, несовпадений {mismatches}, max |Δχ²|={maximum:.6g}.")
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
    ratios = []
    for row in rows(directory / "cells.csv"):
        moments["occupied"] += 1
        raw, program, b = (number(row, field) for field in
                           ("sigma_raw", "sigma_program", "sigma_B"))
        moments["negative_residual"] += number(row, "residual") < 0
        moments["sigma_program_zero"] += program == 0
        moments["sigma_program_nan"] += not math.isfinite(program)
        moments["B_zero_program_positive"] += b == 0 and program > 0
        moments["raw_positive_program_zero"] += raw > 0 and program == 0
        if program > 0 and math.isfinite(b):
            ratio = b/program
            largest(ratios, abs(ratio-1), f"{key(row)} cell={row['cell']}: σ_B/σ_program={ratio:.6g}")
    lines.append(f"Моменты по всем занятым ячейкам (включая вне области фита): {dict(moments)}")
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
        lines.extend("    " + detail for _, detail in sorted(extremes[digits], reverse=True)[:2])
    contributions = collections.defaultdict(lambda: [0, 0., []])
    for row in rows(directory / "chi2_cells.csv"):
        value = number(row, "chi2_contribution")
        stats = contributions[row["variant"]]
        stats[0] += 1
        stats[1] += value
        if math.isfinite(value):
            largest(stats[2], value, f"{key(row)} cell={row['cell']}: вклад={value:.6g}")
    lines.append("Крупнейшие слагаемые χ² (сумма по всем фитам каждого варианта):")
    for variant, (used_count, total, top) in contributions.items():
        lines.append(f"  {variant}: ячеек {used_count}, суммарный χ²={total:.6g}")
        lines.extend("    " + detail for _, detail in sorted(top, reverse=True))
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

#!/usr/bin/env python3
"""Calculate diagonal chi-square from one fit's CSV, without ROOT or refitting."""

import argparse
import csv
import math
from pathlib import Path


def calculate(rows, free_parameters, details=False):
    """Use the supplied model values and errors; exclude bins with sigma=0."""
    used = skipped = 0
    group = None

    def contributions():
        nonlocal used, skipped, group
        for index, row in enumerate(rows, 1):
            # A proof/chi2_cells.csv can contain many different fits.
            current_group = tuple(row.get(field) for field in
                                  ("charge", "centrality", "bin", "variant"))
            if group is None:
                group = current_group
            elif group != current_group:
                raise ValueError("CSV содержит несколько фитов; оставьте один фит")
            try:
                value, model, sigma = (float(row[field]) for field in ("C", "model", "sigma"))
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"Строка {index}: нужны числа C, model, sigma") from error
            if not all(math.isfinite(number) for number in (value, model, sigma)):
                raise ValueError(f"Строка {index}: C, model, sigma должны быть конечными")
            if sigma < 0:
                raise ValueError(f"Строка {index}: sigma должна быть >= 0")
            if sigma == 0:
                skipped += 1
                continue
            pull = (value - model) / sigma
            contribution = pull * pull
            if not math.isfinite(contribution):
                raise ValueError(f"Строка {index}: переполнение вклада chi2")
            if details:
                print(f"{index:6d} {value:.12g} {model:.12g} {sigma:.12g} "
                      f"{pull:.12g} {contribution:.12g}")
            used += 1
            yield contribution

    if details:
        print("row C model sigma pull chi2_contribution")
    chi2 = math.fsum(contributions())
    ndf = used - free_parameters
    print(f"used_cells = {used}; skipped_zero_error = {skipped}")
    print(f"free_parameters = {free_parameters}")
    print(f"chi2 = {chi2:.12g}")
    print(f"ndf = {used} - {free_parameters} = {ndf}")
    if ndf > 0:
        print(f"chi2/ndf = {chi2 / ndf:.12g}")
    else:
        print("chi2/ndf не определено: ndf <= 0")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="?", type=Path,
                        help="CSV одного фита с колонками C,model,sigma; без файла — пример из отчёта")
    parser.add_argument("--free-parameters", type=int,
                        help="число свободных параметров фита (замороженные не считать)")
    parser.add_argument("--details", action="store_true", help="показать вклад каждой строки")
    args = parser.parse_args()
    if args.free_parameters is not None and args.free_parameters < 0:
        parser.error("--free-parameters должно быть >= 0")
    if args.csv is not None and args.free_parameters is None:
        parser.error("для CSV укажите --free-parameters")
    try:
        if args.csv is None:
            rows = ({"C": value, "model": 1.2, "sigma": 0.1}
                    for value in (1.0, 1.1, 1.2, 1.3, 1.4))
            free_parameters = 1 if args.free_parameters is None else args.free_parameters
            calculate(rows, free_parameters, details=True)
        else:
            with args.csv.open(newline="", encoding="utf-8-sig") as source:
                rows = csv.DictReader(source)
                if not {"C", "model", "sigma"}.issubset(rows.fieldnames or []):
                    raise ValueError("CSV должен содержать заголовок C,model,sigma")
                calculate(rows, args.free_parameters, args.details)
    except (OSError, ValueError, OverflowError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()

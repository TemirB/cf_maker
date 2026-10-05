#!/usr/bin/env python3
"""Build a presentation chapter from the exact artifacts of one cf_maker run."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


REPO_ROOT = Path(__file__).resolve().parents[1]
SUPPORTED_IMAGES = {".pdf", ".png", ".jpg", ".jpeg"}
CHARGES = {0: ("pos", r"$\pi^+$"), 1: ("neg", r"$\pi^-$")}
CENTRALITIES = {0: "0-10", 1: "10-30", 2: "30-50", 3: "50-80"}


def repo_path(value):
    path = Path(value).expanduser()
    return (path if path.is_absolute() else REPO_ROOT / path).resolve()


def read_json(path):
    with path.open(encoding="utf-8") as source:
        result = json.load(source)
    if not isinstance(result, dict):
        raise ValueError(f"Ожидался JSON-объект: {path}")
    return result


def output_from_config(config_path):
    config = read_json(config_path)
    base = config["machine"]["base_output"]
    directory = config["vars"]["output"].get("dir", "results")
    if not isinstance(base, str) or not base or not isinstance(directory, str):
        raise ValueError("Неверный путь результатов в конфигурации")
    # Match the C++ string concatenation, including a leading slash in directory.
    return repo_path(base + "/" + directory)


def atomic_write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.read_text(encoding="utf-8") == content:
        return
    descriptor, temporary = tempfile.mkstemp(prefix="." + path.name, dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output:
            output.write(content)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def tex_escape(value):
    replacements = {
        "\\": r"\textbackslash{}", "{": r"\{", "}": r"\}",
        "$": r"\$", "&": r"\&", "#": r"\#", "%": r"\%",
        "_": r"\_", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
        "\n": " ", "\r": " ",
    }
    return "".join(replacements.get(char, char) for char in str(value))


def short_text(value, limit=64):
    text = str(value)
    return tex_escape(text if len(text) <= limit else text[:limit - 1] + "…")


def load_manifest(results):
    path = results / "run_manifest.json"
    raw = path.read_bytes()
    manifest = json.loads(raw)
    if not isinstance(manifest, dict) or manifest.get("schema_version") != 1:
        raise ValueError(f"Неподдерживаемый manifest: {path}")
    status = manifest.get("status")
    if status not in ("completed", "incomplete"):
        raise ValueError("Расчёт не завершён. Перезапустите программу перед сборкой результатов.")
    expected_code = 0 if status == "completed" else 2
    if manifest.get("exit_code") != expected_code or not manifest.get("run_id"):
        raise ValueError("Статус расчёта не согласован с кодом завершения")
    config = manifest["config"]
    if repo_path(config["vars"]["output"]["dir"]) != results:
        raise ValueError("Manifest относится к другой директории результатов")
    fits = manifest["fits"]
    for key in ("requested", "usable", "unusable", "retried", "at_limit"):
        if type(fits.get(key)) is not int or fits[key] < 0:
            raise ValueError(f"Неверный счётчик фитов: {key}")
    if fits["requested"] != fits["usable"] + fits["unusable"]:
        raise ValueError("Счётчики фитов не согласованы")
    if status == "completed" and fits["unusable"]:
        raise ValueError("Успешный статус содержит непригодные фиты")
    if not isinstance(manifest.get("images"), list):
        raise ValueError("В manifest отсутствует список рисунков")
    return manifest, hashlib.sha256(raw).hexdigest()


def plot_labels(config):
    """Map exported filenames to captions without inferring scientific results."""
    labels = {}
    selection = config["selection"]
    names = config["binning"]["names"]
    file_names = config["binning"]["file_names"]
    if len(names) != len(file_names):
        raise ValueError("Подписи и имена бинов не согласованы")
    mode = r"$k_t$" if config["vars"]["input"]["type"] == "kt" else "$y$"
    for ch in selection["charges"]:
        slug, charge = CHARGES[ch]
        for stem, category, caption in (
            ("c_all_graphs", "radii", "Радиусы и сила корреляции"),
            ("c_fit_quality", "quality", "Качество аппроксимации"),
            ("c_pvalues", "pvalues", "Вероятности согласия"),
            ("c_all_cross_graphs", "cross", "Перекрёстные параметры"),
        ):
            labels[f"dependency/{stem}_{slug}"] = (
                category, caption, charge + ", зависимость от " + mode)
        for centr in selection["centralities"]:
            centrality = CENTRALITIES[centr]
            sample = charge + ", центральность " + centrality + r"\%"
            labels[f"all_2d_histos/all_out-long_2d_histos_centr_{centrality}_{slug}"] = (
                "2d", "Двумерные проекции КФ", sample + ", плоскость out–long")
            for name, filename in zip(names, file_names):
                caption = sample + ", " + mode + ": " + short_text(name, 35)
                labels[f"all_1d_histos/cfs_{slug}_{centrality}_{filename}"] = (
                    "1d", "Одномерные проекции КФ", caption)
                labels[f"dependency/fit_over_cf/fit_over_cf_{slug}_{centrality}_{filename}"] = (
                    "fit_ratio", "Отношение аппроксимации к КФ", caption)
    return labels


def validate_images(results, manifest):
    labels = plot_labels(manifest["config"])
    images = []
    seen = set()
    for entry in manifest["images"]:
        relative = Path(entry["path"])
        image = (results / relative).resolve()
        if relative.is_absolute() or image == results or results not in image.parents:
            raise ValueError(f"Рисунок находится вне результатов: {relative}")
        if relative.as_posix() in seen:
            raise ValueError(f"Рисунок перечислен дважды: {relative}")
        seen.add(relative.as_posix())
        if image.suffix.lower() not in SUPPORTED_IMAGES:
            raise ValueError(f"Формат {image.suffix} не поддерживается. Задайте images.format=pdf.")
        if not image.is_file() or image.stat().st_size <= 0:
            raise ValueError(f"Нет рисунка текущего расчёта: {image}")
        if image.stat().st_size != entry["size_bytes"]:
            raise ValueError(f"Рисунок изменён после расчёта: {image}")
        stem = relative.with_suffix("").as_posix()
        category, title, caption = labels.get(stem, (
            "other", "Результат расчёта", short_text(relative.name)))
        images.append({"source": image, "path": relative.as_posix(), "category": category,
                       "title": title, "caption": caption})
    return images


def select_plots(images, maximum):
    if maximum == 0 or len(images) <= maximum:
        return images
    # Give each kind of graph room before adding more per-bin projections.
    categories = ("radii", "quality", "1d", "2d", "pvalues", "fit_ratio", "cross", "other")
    groups = [[image for image in images if image["category"] == category]
              for category in categories]
    selected = []
    while len(selected) < maximum:
        added = False
        for group in groups:
            if group and len(selected) < maximum:
                selected.append(group.pop(0))
                added = True
        if not added:
            break
    return selected


def render_chapter(manifest, images, total):
    config = manifest["config"]
    fits = manifest["fits"]
    diagnostic = manifest["status"] == "incomplete"
    status = (r"\alert{Диагностический расчёт: есть непригодные фиты.}"
              if diagnostic else "Запрошенные стадии завершены.")
    mode = r"$k_t$" if config["vars"]["input"]["type"] == "kt" else "$y$"
    charges = ", ".join(CHARGES[ch][1] for ch in config["selection"]["charges"])
    centralities = ", ".join(CENTRALITIES[c] + r"\%"
                            for c in config["selection"]["centralities"])
    fit_line = (f"Пригодные фиты: {fits['usable']} из {fits['requested']}; "
                f"непригодные: {fits['unusable']}." if fits["requested"]
                else "Аппроксимация в выбранных стадиях не запрашивалась.")
    lines = [r"\section{Результаты текущего расчёта}",
             r"\begin{frame}{Параметры и статус расчёта}", r"\small", status,
             r"\begin{itemize}",
             r"\item Вход: " + short_text(Path(config["vars"]["input"]["file"]).name),
             r"\item Завершение (UTC): " + tex_escape(manifest["completed_at"]),
             r"\item Биннинг: " + mode + f"; интервалов: {len(config['binning']['names'])}.",
             r"\item Заряды: " + charges + "; центральности: " + centralities + ".",
             r"\item " + fit_line,
             r"\item Повторные попытки: " + str(fits["retried"]) +
             "; параметры на границах: " + str(fits["at_limit"]) + ".",
             r"\item Рисунков в презентации: " + str(len(images)) +
             " из " + str(total) + " сохранённых в этом расчёте.",
             r"\end{itemize}", r"\end{frame}"]
    if not images:
        lines += [r"\begin{frame}{Рисунки расчёта}",
                  "В этом запуске рисунки не экспортировались.",
                  r"\par\medskip Для экспорта включите \texttt{images.need=true} "
                  r"и выберите формат \texttt{pdf}.", r"\end{frame}"]
    for image in images:
        title = image["title"] + (" (диагностика)" if diagnostic else "")
        lines += [r"\begin{frame}{" + tex_escape(title) + "}", r"\centering",
                  r"\includegraphics[width=.96\linewidth,height=.72\textheight,keepaspectratio]{" +
                  image["asset"] + "}", r"\par\smallskip",
                  r"{\footnotesize " + image["caption"] + "}", r"\end{frame}"]
    return "\n".join(lines) + "\n"


def prepare(args):
    output = repo_path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    tex_path = output / "presentation-results.tex"
    meta_path = output / "presentation-results.json"
    try:
        results = repo_path(args.results) if args.results else output_from_config(repo_path(args.config))
        if not (results / "run_manifest.json").is_file():
            if args.results:
                raise ValueError("В выбранных результатах нет run_manifest.json. Перезапустите cf_maker.")
            atomic_write(tex_path, "")
            atomic_write(meta_path, "{}\n")
            print("Результаты ещё не подготовлены; собираются исходные слайды. "
                  "Для обновления результатов выполните make run-presentation.")
            return 0
        manifest, digest = load_manifest(results)
        images = validate_images(results, manifest)
        selected = select_plots(images, args.max_plots)
        previous = read_json(meta_path) if meta_path.is_file() else {}
        if (previous.get("manifest_sha256") == digest and
                previous.get("max_plots") == args.max_plots and
                previous.get("results") == str(results) and tex_path.is_file() and
                previous.get("tex_sha256") == hashlib.sha256(tex_path.read_bytes()).hexdigest() and
                len(previous.get("assets", [])) == len(selected) and
                all((output / path).is_file() and
                    (output / path).stat().st_size == image["source"].stat().st_size
                    for path, image in zip(previous.get("assets", []), selected))):
            print(f"Слайды результатов актуальны: {len(selected)} рисунков.")
            return 0
        assets_root = output / "current-results"
        assets_root.mkdir(exist_ok=True)
        assets = Path(tempfile.mkdtemp(prefix="run-", dir=assets_root))
        for i, image in enumerate(selected, 1):
            destination = assets / f"plot-{i:04d}{image['source'].suffix.lower()}"
            shutil.copyfile(image["source"], destination)
            image["asset"] = destination.relative_to(output).as_posix()
        if load_manifest(results)[1] != digest:
            raise ValueError("Результаты обновились во время подготовки. Повторите сборку.")
        tex = render_chapter(manifest, selected, len(images))
        metadata = {"results": str(results), "run_id": manifest["run_id"],
                    "status": manifest["status"], "manifest_sha256": digest,
                    "max_plots": args.max_plots, "total_images": len(images),
                    "assets": [image["asset"] for image in selected],
                    "tex_sha256": hashlib.sha256(tex.encode()).hexdigest()}
        atomic_write(tex_path, tex)
        atomic_write(meta_path, json.dumps(metadata, ensure_ascii=False, indent=2) + "\n")
        print(f"Подготовлены слайды: {len(selected)} из {len(images)} рисунков, "
              f"статус {manifest['status']}.")
        return 0
    except Exception:
        # A later standalone build must never pick up an earlier run's chapter.
        atomic_write(tex_path, "")
        atomic_write(meta_path, "{}\n")
        raise


def publish(args):
    output = repo_path(args.output_dir)
    pdf = output / "presentation.pdf"
    if not pdf.is_file() or pdf.stat().st_size == 0:
        raise ValueError("Сборка не создала presentation.pdf")
    metadata = read_json(output / "presentation-results.json")
    if not metadata:
        print(f"PDF: {pdf}")
        return 0
    results = Path(metadata["results"])
    manifest, digest = load_manifest(results)
    if digest != metadata["manifest_sha256"]:
        raise ValueError("Расчёт изменился во время сборки; PDF не опубликован в его результатах.")
    destination = results / "presentation.pdf"
    descriptor, temporary = tempfile.mkstemp(prefix=".presentation-", suffix=".pdf", dir=results)
    os.close(descriptor)
    try:
        shutil.copyfile(pdf, temporary)
        if load_manifest(results)[1] != digest:
            raise ValueError("Расчёт изменился во время копирования; PDF не опубликован.")
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    if manifest["status"] == "incomplete":
        print("ДИАГНОСТИЧЕСКИЙ PDF: анализ вернул код 2; непригодные фиты отмечены в слайдах.")
    print(f"PDF: {destination}")
    return 0


def run(args):
    config = repo_path(args.config)
    results = output_from_config(config)
    previous_id = None
    try:
        previous_id = read_json(results / "run_manifest.json").get("run_id")
    except (OSError, ValueError):
        pass
    result = subprocess.run([str(repo_path(args.executable)), str(config)], cwd=REPO_ROOT)
    if result.returncode not in (0, 2):
        print("Расчёт завершился ошибкой; сборка результатов остановлена.", file=sys.stderr)
        return result.returncode if result.returncode > 0 else 1
    manifest, _ = load_manifest(results)
    if manifest["run_id"] == previous_id:
        raise ValueError("Программа не сохранила результаты нового запуска; сборка остановлена.")
    if manifest["exit_code"] != result.returncode:
        raise ValueError("Код программы не совпадает со статусом сохранённых результатов")
    if result.returncode == 2:
        print("Анализ вернул код 2. Будет собран PDF с явной отметкой диагностики.")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--config", default="config/kt.json")
    prepare_parser.add_argument("--results")
    prepare_parser.add_argument("--output-dir", required=True)
    prepare_parser.add_argument("--max-plots", type=int, default=12)
    publish_parser = commands.add_parser("publish")
    publish_parser.add_argument("--output-dir", required=True)
    run_parser = commands.add_parser("run")
    run_parser.add_argument("--config", required=True)
    run_parser.add_argument("--executable", required=True)
    args = parser.parse_args()
    if args.command == "prepare" and args.max_plots < 0:
        parser.error("--max-plots должен быть неотрицательным; 0 включает все рисунки")
    try:
        return {"prepare": prepare, "publish": publish, "run": run}[args.command](args)
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"Презентация: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

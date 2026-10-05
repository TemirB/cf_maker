#!/usr/bin/env python3
"""Regression checks for generating presentation chapters from run artifacts."""

import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_SPEC = importlib.util.spec_from_file_location(
    "prepare_presentation", REPO_ROOT / "tools" / "prepare_presentation.py"
)
PRESENTATION = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(PRESENTATION)


class PresentationTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="cf-maker-presentation-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name).resolve()
        self.results = self.directory / "results"
        self.results.mkdir()
        self.output = self.directory / "slides"
        self.config_path = self.directory / "config.json"
        self.config = {
            "machine": {
                "base_input": str(self.directory),
                "base_output": str(self.directory),
            },
            "vars": {
                "input": {"type": "kt", "file": "sample_&%{run}.root"},
                "output": {"dir": "results"},
            },
            "selection": {"charges": [0, 1], "centralities": [0]},
            "binning": {
                "names": ["k_t&50%{A}", "second", "third", "fourth"],
                "file_names": ["special_bin", "second", "third", "fourth"],
            },
        }
        self.write_json(self.config_path, self.config)
        snapshot = json.loads(json.dumps(self.config))
        snapshot["vars"]["input"]["file"] = str(
            self.directory / self.config["vars"]["input"]["file"]
        )
        snapshot["vars"]["output"]["dir"] = str(self.results)
        self.manifest = {
            "schema_version": 1,
            "run_id": "test-run",
            "status": "completed",
            "started_at": "2026-10-06T10:00:00.000Z",
            "completed_at": "2026-10-06T10:01:00.000Z",
            "exit_code": 0,
            "config": snapshot,
            "fits": {
                "requested": 8,
                "usable": 8,
                "unusable": 0,
                "retried": 1,
                "at_limit": 0,
            },
            "images": [],
        }

    @staticmethod
    def write_json(path, value):
        path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")

    def write_manifest(self):
        self.write_json(self.results / "run_manifest.json", self.manifest)

    def add_image(self, relative, content=b"test image content"):
        image = self.results / relative
        image.parent.mkdir(parents=True, exist_ok=True)
        image.write_bytes(content)
        self.manifest["images"].append(
            {"path": relative, "size_bytes": len(content)}
        )
        return image

    def prepare_args(self, maximum=12, explicit=True):
        return SimpleNamespace(
            config=str(self.config_path),
            results=str(self.results) if explicit else None,
            output_dir=str(self.output),
            max_plots=maximum,
        )

    def prepare(self, maximum=12, explicit=True):
        with contextlib.redirect_stdout(io.StringIO()):
            return PRESENTATION.prepare(self.prepare_args(maximum, explicit))

    def metadata(self):
        return json.loads(
            (self.output / "presentation-results.json").read_text(encoding="utf-8")
        )

    def chapter(self):
        return (self.output / "presentation-results.tex").read_text(encoding="utf-8")

    def assert_invalidated(self):
        self.assertEqual(self.chapter(), "")
        self.assertEqual(self.metadata(), {})

    def test_current_images_only_and_unchanged_prepare_cache(self):
        current = self.add_image("dependency/c_all_graphs_pos.png", b"current plot")
        stale = self.results / "stale.png"
        stale.write_bytes(b"old plot")
        self.write_manifest()
        self.assertEqual(self.prepare(), 0)
        metadata = self.metadata()
        self.assertEqual(metadata["total_images"], 1)
        self.assertEqual(len(metadata["assets"]), 1)
        asset = self.output / metadata["assets"][0]
        self.assertEqual(asset.read_bytes(), current.read_bytes())
        self.assertNotIn("stale", self.chapter())
        self.assertEqual(stale.read_bytes(), b"old plot")
        paths = [self.output / "presentation-results.tex",
                 self.output / "presentation-results.json", asset]
        before = {path: path.stat().st_mtime_ns for path in paths}
        asset_directories = list((self.output / "current-results").iterdir())
        self.assertEqual(self.prepare(), 0)
        self.assertEqual(self.metadata(), metadata)
        self.assertEqual(before, {path: path.stat().st_mtime_ns for path in paths})
        self.assertEqual(list((self.output / "current-results").iterdir()), asset_directories)

    def test_plot_limit_balances_categories_and_zero_includes_all(self):
        for charge in ("pos", "neg"):
            for filename in self.config["binning"]["file_names"]:
                self.add_image(f"all_1d_histos/cfs_{charge}_0-10_{filename}.png")
        for name in ("c_all_graphs_pos", "c_fit_quality_pos", "c_pvalues_pos"):
            self.add_image(f"dependency/{name}.png")
        self.add_image("all_2d_histos/all_out-long_2d_histos_centr_0-10_pos.png")
        self.write_manifest()
        self.assertEqual(self.prepare(maximum=4), 0)
        self.assertEqual(len(self.metadata()["assets"]), 4)
        for title in ("Радиусы и сила корреляции", "Качество аппроксимации",
                      "Одномерные проекции КФ", "Двумерные проекции КФ"):
            self.assertIn(title, self.chapter())
        self.assertEqual(self.prepare(maximum=0), 0)
        self.assertEqual(self.metadata()["total_images"], 12)
        self.assertEqual(len(self.metadata()["assets"]), 12)
        self.assertEqual(self.chapter().count(r"\includegraphics"), 12)

    def test_prepare_rebuilds_damaged_cached_asset(self):
        source = self.add_image("dependency/c_all_graphs_pos.png", b"complete current plot")
        self.write_manifest()
        self.prepare()
        old_asset = self.metadata()["assets"][0]
        (self.output / old_asset).write_bytes(b"truncated")
        self.assertEqual(self.prepare(), 0)
        new_asset = self.metadata()["assets"][0]
        self.assertNotEqual(new_asset, old_asset)
        self.assertEqual((self.output / new_asset).read_bytes(), source.read_bytes())
        self.assertIn(new_asset, self.chapter())
        self.assertNotIn(old_asset, self.chapter())

    def test_generated_input_and_bin_labels_escape_tex(self):
        self.add_image("all_1d_histos/cfs_pos_0-10_special_bin.png")
        self.write_manifest()
        self.prepare()
        self.assertIn(r"sample\_\&\%\{run\}.root", self.chapter())
        self.assertIn(r"k\_t\&50\%\{A\}", self.chapter())
        self.assertNotIn("sample_&%{run}.root", self.chapter())
        self.assertNotIn("k_t&50%{A}", self.chapter())
        self.assertEqual(
            PRESENTATION.tex_escape(r"\{}$&#%_~^"),
            r"\textbackslash{}\{\}\$\&\#\%\_\textasciitilde{}\textasciicircum{}",
        )

    def test_running_manifest_invalidates_previous_chapter(self):
        self.write_manifest()
        self.prepare()
        self.assertTrue(self.chapter())
        self.manifest.update(status="running", exit_code=None, completed_at=None)
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, "не завершён"):
            self.prepare()
        self.assert_invalidated()

    def test_explicit_missing_manifest_invalidates_previous_chapter(self):
        self.write_manifest()
        self.prepare()
        (self.results / "run_manifest.json").unlink()
        with self.assertRaisesRegex(ValueError, "нет run_manifest.json"):
            self.prepare()
        self.assert_invalidated()

    def test_missing_listed_image_invalidates_previous_chapter(self):
        image = self.add_image("dependency/c_all_graphs_pos.png")
        self.write_manifest()
        self.prepare()
        image.unlink()
        with self.assertRaisesRegex(ValueError, "Нет рисунка текущего расчёта"):
            self.prepare()
        self.assert_invalidated()

    def test_missing_config_invalidates_previous_chapter(self):
        self.write_manifest()
        self.prepare(explicit=False)
        self.assertTrue(self.chapter())
        self.config_path.unlink()
        with self.assertRaises(FileNotFoundError):
            self.prepare(explicit=False)
        self.assert_invalidated()

    def test_malformed_config_invalidates_previous_chapter(self):
        self.write_manifest()
        self.prepare(explicit=False)
        self.assertTrue(self.chapter())
        self.config_path.write_text("{not valid json", encoding="utf-8")
        with self.assertRaises(json.JSONDecodeError):
            self.prepare(explicit=False)
        self.assert_invalidated()

    def test_automatic_missing_manifest_prepares_empty_chapter(self):
        self.write_manifest()
        self.prepare(explicit=False)
        self.assertTrue(self.chapter())
        (self.results / "run_manifest.json").unlink()
        self.assertEqual(self.prepare(explicit=False), 0)
        self.assert_invalidated()

    def test_publish_copies_pdf_and_rejects_changed_run(self):
        self.write_manifest()
        self.prepare()
        pdf = self.output / "presentation.pdf"
        pdf.write_bytes(b"%PDF test presentation")
        args = SimpleNamespace(output_dir=str(self.output))
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(PRESENTATION.publish(args), 0)
        published = self.results / "presentation.pdf"
        self.assertEqual(published.read_bytes(), pdf.read_bytes())
        self.manifest["run_id"] = "next-run"
        self.write_manifest()
        pdf.write_bytes(b"%PDF changed presentation")
        with self.assertRaisesRegex(ValueError, "Расчёт изменился"):
            PRESENTATION.publish(args)
        self.assertEqual(published.read_bytes(), b"%PDF test presentation")

    def test_publish_rejects_run_changed_while_copying(self):
        self.write_manifest()
        self.prepare()
        (self.output / "presentation.pdf").write_bytes(b"%PDF current presentation")
        published = self.results / "presentation.pdf"
        published.write_bytes(b"%PDF preserve published presentation")
        copy_file = PRESENTATION.shutil.copyfile

        def copy_and_change_manifest(source, destination):
            copy_file(source, destination)
            self.manifest["run_id"] = "next-run-during-copy"
            self.write_manifest()

        with mock.patch.object(PRESENTATION.shutil, "copyfile", copy_and_change_manifest):
            with self.assertRaisesRegex(ValueError, "изменился во время копирования"):
                PRESENTATION.publish(SimpleNamespace(output_dir=str(self.output)))
        self.assertEqual(published.read_bytes(), b"%PDF preserve published presentation")
        self.assertEqual(list(self.results.glob(".presentation-*")), [])

    def fake_executable(self, exit_code, update_manifest=False):
        executable = self.directory / f"fake-analysis-{exit_code}"
        body = ""
        if update_manifest:
            body = (
                "import json\nfrom pathlib import Path\n"
                f"path = Path({str(self.results / 'run_manifest.json')!r})\n"
                "value = json.loads(path.read_text(encoding='utf-8'))\n"
                "value['run_id'] += '-new'\n"
                "path.write_text(json.dumps(value), encoding='utf-8')\n"
            )
        executable.write_text(
            f"#!{sys.executable}\nimport sys\n{body}sys.exit({exit_code})\n",
            encoding="utf-8",
        )
        executable.chmod(0o755)
        return executable

    def run_args(self, exit_code, update_manifest=False):
        return SimpleNamespace(config=str(self.config_path),
                               executable=str(self.fake_executable(exit_code, update_manifest)))

    def test_run_wrapper_accepts_incomplete_and_marks_diagnostics(self):
        self.manifest.update(status="incomplete", exit_code=2)
        self.manifest["fits"].update(usable=6, unusable=2, at_limit=1)
        self.write_manifest()
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertEqual(PRESENTATION.run(self.run_args(2, update_manifest=True)), 0)
        self.assertIn("код 2", output.getvalue())
        self.prepare()
        self.assertIn("Диагностический расчёт", self.chapter())
        self.assertIn("Пригодные фиты: 6 из 8; непригодные: 2.", self.chapter())

    def test_run_wrapper_stops_on_fatal_error(self):
        # A fatal CLI error must stop before inspecting a missing/invalid manifest.
        errors = io.StringIO()
        with contextlib.redirect_stderr(errors):
            self.assertEqual(PRESENTATION.run(self.run_args(1)), 1)
        self.assertIn("сборка результатов остановлена", errors.getvalue())
        self.assertFalse(self.output.exists())

    def test_run_wrapper_rejects_exit_code_manifest_mismatch(self):
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, "Код программы не совпадает"):
            PRESENTATION.run(self.run_args(2, update_manifest=True))

    def test_run_wrapper_rejects_unchanged_previous_run(self):
        for exit_code, status in ((0, "completed"), (2, "incomplete")):
            with self.subTest(exit_code=exit_code):
                self.manifest.update(status=status, exit_code=exit_code)
                self.write_manifest()
                with self.assertRaisesRegex(ValueError, "результаты нового запуска"):
                    PRESENTATION.run(self.run_args(exit_code))
                self.assertEqual(
                    json.loads((self.results / "run_manifest.json").read_text())["run_id"],
                    "test-run",
                )


if __name__ == "__main__":
    unittest.main()

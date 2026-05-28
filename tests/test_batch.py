import json
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.core.batch import batch_summary, read_manifest, run_batch
from biocompare.io.report_writers import batch_to_tsv


class BatchTests(TestCase):
    def test_manifest_paths_are_resolved_relative_to_manifest(self):
        items = read_manifest("tests/fixtures/batch_manifest.tsv")

        self.assertEqual(len(items), 3)
        self.assertTrue(items[0].file_a.endswith("tests/fixtures/table_a.tsv"))
        self.assertEqual(items[1].file_type, "deg")

    def test_run_batch_returns_reports_and_summary(self):
        results = run_batch("tests/fixtures/batch_manifest.tsv")
        summary = batch_summary(results)

        self.assertEqual(summary["total"], 3)
        self.assertEqual(summary["ok"], 3)
        self.assertEqual(summary["error"], 0)
        self.assertIn("overall_concordance", results[0].report.to_dict())

    def test_batch_tsv_writer_includes_one_row_per_result(self):
        results = run_batch("tests/fixtures/batch_manifest.tsv")
        output = batch_to_tsv(results)

        self.assertIn("label\tstatus\tcomparator", output)
        self.assertEqual(len(output.splitlines()), 4)


class BatchCliTests(TestCase):
    def test_cli_batch_outputs_json_summary(self):
        result = subprocess.run(
            [sys.executable, "-m", "biocompare", "batch", "tests/fixtures/batch_manifest.tsv", "--format", "json"],
            check=True,
            capture_output=True,
            text=True,
        )

        payload = json.loads(result.stdout)
        self.assertEqual(payload["summary"]["total"], 3)
        self.assertEqual(payload["summary"]["ok"], 3)

    def test_cli_batch_writes_output_file(self):
        with TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "batch.tsv"
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "biocompare",
                    "batch",
                    "tests/fixtures/batch_manifest.tsv",
                    "--output",
                    str(output_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            output = output_path.read_text(encoding="utf-8")

        self.assertIn("table\tok", output)

    def test_cli_batch_min_concordance_controls_exit_code(self):
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "biocompare",
                "batch",
                "tests/fixtures/batch_manifest.tsv",
                "--min-concordance",
                "0.95",
            ],
            capture_output=True,
            text=True,
        )

        self.assertEqual(result.returncode, 1)

    def test_cli_batch_outputs_html_summary(self):
        result = subprocess.run(
            [sys.executable, "-m", "biocompare", "batch", "tests/fixtures/batch_manifest.tsv", "--format", "html"],
            check=True,
            capture_output=True,
            text=True,
        )

        self.assertIn("<!doctype html>", result.stdout)
        self.assertIn("biocompare batch report", result.stdout)

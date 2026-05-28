import json
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase


class CliTests(TestCase):
    def test_cli_outputs_json_report(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text("id\tvalue\nA\t1\n", encoding="utf-8")
            right.write_text("id\tvalue\nA\t1\n", encoding="utf-8")
            result = subprocess.run(
                [sys.executable, "-m", "biocompare", "compare", str(left), str(right), "--key", "id"],
                check=True,
                capture_output=True,
                text=True,
            )

        report = json.loads(result.stdout)
        self.assertEqual(report["comparator"], "TableComparator")
        self.assertEqual(report["overall_concordance"], 1.0)

    def test_legacy_pair_invocation_still_works(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text("id\tvalue\nA\t1\n", encoding="utf-8")
            right.write_text("id\tvalue\nA\t1\n", encoding="utf-8")
            result = subprocess.run(
                [sys.executable, "-m", "biocompare", str(left), str(right), "--key", "id"],
                check=True,
                capture_output=True,
                text=True,
            )

        report = json.loads(result.stdout)
        self.assertEqual(report["comparator"], "TableComparator")

    def test_cli_outputs_html_report(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text("id\tvalue\nA\t1\n", encoding="utf-8")
            right.write_text("id\tvalue\nA\t1\n", encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "biocompare",
                    "compare",
                    str(left),
                    str(right),
                    "--key",
                    "id",
                    "--format",
                    "html",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

        self.assertIn("<!doctype html>", result.stdout)
        self.assertIn("TableComparator", result.stdout)

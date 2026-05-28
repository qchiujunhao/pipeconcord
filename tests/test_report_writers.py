from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.core.batch import run_batch
from biocompare.core.engine import ComparisonEngine
from biocompare.io.report_writers import batch_to_html, report_to_html, write_batch, write_report


class ReportWriterTests(TestCase):
    def test_single_report_html_contains_metric_table(self):
        report = ComparisonEngine().compare(
            "tests/fixtures/table_a.tsv",
            "tests/fixtures/table_b.tsv",
            file_type="table",
            key_column="gene_id",
        )

        html = report_to_html(report)

        self.assertIn("<!doctype html>", html)
        self.assertIn("TableComparator", html)
        self.assertIn("row_overlap", html)

    def test_write_report_uses_html_extension(self):
        report = ComparisonEngine().compare(
            "tests/fixtures/table_a.tsv",
            "tests/fixtures/table_b.tsv",
            file_type="table",
            key_column="gene_id",
        )
        with TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "report.html"
            write_report(report, str(output_path))
            html = output_path.read_text(encoding="utf-8")

        self.assertIn("Overall concordance", html)

    def test_batch_html_contains_summary_and_results(self):
        results = run_batch("tests/fixtures/batch_manifest.tsv")

        html = batch_to_html(results)

        self.assertIn("biocompare batch report", html)
        self.assertIn("Mean concordance", html.replace("_", " "))
        self.assertIn("TableComparator", html)

    def test_write_batch_uses_html_extension(self):
        results = run_batch("tests/fixtures/batch_manifest.tsv")
        with TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "batch.html"
            write_batch(results, str(output_path))
            html = output_path.read_text(encoding="utf-8")

        self.assertIn("Batch Comparison", html)


from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.comparators.bed import BEDComparator
from biocompare.core.engine import ComparisonEngine


class BEDComparatorTests(TestCase):
    def test_identical_bed_files_are_fully_concordant(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.bed"
            right = Path(temp_dir) / "right.bed"
            content = "chr1\t0\t100\nchr1\t200\t300\nchr2\t0\t50\n"
            left.write_text(content, encoding="utf-8")
            right.write_text(content, encoding="utf-8")

            report = BEDComparator().compare(str(left), str(right))

        self.assertEqual(report.overall_concordance, 1.0)
        self.assertEqual(report.metrics["bp_jaccard"], 1.0)
        self.assertEqual(report.metrics["interval_f1_file_b_vs_a"], 1.0)
        self.assertEqual(report.details["coordinate_system"], "BED 0-based half-open")

    def test_bed_metrics_capture_partial_overlap(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.bed"
            right = Path(temp_dir) / "right.bed"
            left.write_text("chr1\t0\t100\nchr1\t200\t300\nchr2\t0\t50\n", encoding="utf-8")
            right.write_text("chr1\t50\t150\nchr1\t200\t300\nchr2\t50\t100\n", encoding="utf-8")

            report = BEDComparator().compare(str(left), str(right))

        self.assertAlmostEqual(report.metrics["bp_jaccard"], 150 / 350)
        self.assertAlmostEqual(report.metrics["interval_precision_file_b_vs_a"], 2 / 3)
        self.assertAlmostEqual(report.metrics["interval_recall_file_b_vs_a"], 2 / 3)
        self.assertLess(report.overall_concordance, 1.0)

    def test_reciprocal_overlap_threshold_changes_interval_matching(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.bed"
            right = Path(temp_dir) / "right.bed"
            left.write_text("chr1\t0\t100\nchr1\t200\t300\nchr2\t0\t50\n", encoding="utf-8")
            right.write_text("chr1\t50\t150\nchr1\t200\t300\nchr2\t50\t100\n", encoding="utf-8")

            report = BEDComparator().compare(str(left), str(right), min_reciprocal_overlap=0.6)

        self.assertAlmostEqual(report.metrics["interval_precision_file_b_vs_a"], 1 / 3)
        self.assertAlmostEqual(report.metrics["interval_recall_file_b_vs_a"], 1 / 3)

    def test_engine_auto_selects_bed_for_bed_extension(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.bed"
            right = Path(temp_dir) / "right.bed"
            left.write_text("chr1\t0\t100\n", encoding="utf-8")
            right.write_text("chr1\t50\t150\n", encoding="utf-8")

            report = ComparisonEngine().compare(str(left), str(right))

        self.assertEqual(report.comparator, "BEDComparator")
        self.assertIn("bp_jaccard", report.metrics)


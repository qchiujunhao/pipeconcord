from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.comparators.counts import CountsComparator
from biocompare.core.engine import ComparisonEngine


class CountsComparatorTests(TestCase):
    def test_identical_count_matrices_are_fully_concordant(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            content = "gene_id\ts1\ts2\nG1\t10\t5\nG2\t0\t20\nG3\t30\t15\n"
            left.write_text(content, encoding="utf-8")
            right.write_text(content, encoding="utf-8")

            report = CountsComparator().compare(str(left), str(right))

        self.assertEqual(report.overall_concordance, 1.0)
        self.assertEqual(report.metrics["gene_overlap"], 1.0)
        self.assertEqual(report.metrics["sample_overlap"], 1.0)
        self.assertEqual(report.metrics["mean_sample_similarity"], 1.0)
        self.assertEqual(report.metrics["zero_pattern_jaccard"], 1.0)

    def test_count_matrix_metrics_capture_gene_and_value_differences(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text(
                "gene_id\ts1\ts2\ts3\nG1\t100\t50\t0\nG2\t20\t80\t5\nG3\t0\t10\t30\nG4\t300\t200\t100\n",
                encoding="utf-8",
            )
            right.write_text(
                "gene_id\ts1\ts2\ts3\nG1\t110\t55\t0\nG2\t18\t75\t5\nG3\t0\t12\t28\nG5\t400\t250\t120\n",
                encoding="utf-8",
            )

            report = CountsComparator().compare(str(left), str(right))

        self.assertAlmostEqual(report.metrics["gene_overlap"], 3 / 5)
        self.assertEqual(report.details["shared_genes"], 3)
        self.assertEqual(report.details["shared_samples"], ["s1", "s2", "s3"])
        self.assertLess(report.overall_concordance, 1.0)
        self.assertIn("sample.s1.pearson", report.metrics)

    def test_engine_auto_selects_counts_before_generic_table(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.csv"
            right = Path(temp_dir) / "right.csv"
            left.write_text("gene_id,s1,s2\nG1,10,5\nG2,0,20\n", encoding="utf-8")
            right.write_text("gene_id,s1,s2\nG1,11,5\nG2,0,19\n", encoding="utf-8")

            report = ComparisonEngine().compare(str(left), str(right))

        self.assertEqual(report.comparator, "CountsComparator")
        self.assertIn("mean_sample_pearson", report.metrics)

    def test_sample_column_override(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text("feature\ts1\ts2\tlength\nG1\t10\t5\t1000\nG2\t0\t20\t800\n", encoding="utf-8")
            right.write_text("feature\ts1\ts2\tlength\nG1\t10\t5\t1100\nG2\t0\t20\t900\n", encoding="utf-8")

            report = CountsComparator().compare(str(left), str(right), sample_columns="s1,s2")

        self.assertEqual(report.details["shared_samples"], ["s1", "s2"])
        self.assertEqual(report.overall_concordance, 1.0)

    def test_all_zero_shared_sample_has_full_library_size_ratio(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text("gene_id\ts1\ts2\nG1\t0\t5\nG2\t0\t10\n", encoding="utf-8")
            right.write_text("gene_id\ts1\ts2\nG1\t0\t5\nG2\t0\t10\n", encoding="utf-8")

            report = CountsComparator().compare(str(left), str(right))

        self.assertEqual(report.metrics["sample.s1.library_size_ratio"], 1.0)

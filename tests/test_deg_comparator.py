from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.comparators.deg import DEGComparator
from biocompare.core.engine import ComparisonEngine


class DEGComparatorTests(TestCase):
    def test_identical_deg_tables_are_fully_concordant(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            content = "gene_id\tlog2FoldChange\tpadj\nG1\t2.0\t0.001\nG2\t-1.0\t0.020\nG3\t0.1\t0.800\n"
            left.write_text(content, encoding="utf-8")
            right.write_text(content, encoding="utf-8")

            report = DEGComparator().compare(str(left), str(right), alpha=0.05)

        self.assertEqual(report.overall_concordance, 1.0)
        self.assertEqual(report.metrics["gene_overlap"], 1.0)
        self.assertEqual(report.metrics["significant_jaccard"], 1.0)
        self.assertEqual(report.metrics["significant_direction_agreement"], 1.0)
        self.assertEqual(report.details["file_a_columns"]["uses_adjusted_pvalue"], True)

    def test_deg_metrics_capture_significant_set_and_direction_changes(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text(
                "gene_id\tlog2FoldChange\tpadj\nG1\t2.0\t0.001\nG2\t-1.5\t0.010\nG3\t0.2\t0.700\nG4\t1.0\t0.040\n",
                encoding="utf-8",
            )
            right.write_text(
                "gene_id\tlog2FoldChange\tpadj\nG1\t1.8\t0.002\nG2\t1.2\t0.020\nG3\t0.1\t0.800\nG5\t3.0\t0.001\n",
                encoding="utf-8",
            )

            report = DEGComparator().compare(str(left), str(right), alpha=0.05)

        self.assertAlmostEqual(report.metrics["gene_overlap"], 3 / 5)
        self.assertAlmostEqual(report.metrics["significant_jaccard"], 2 / 4)
        self.assertAlmostEqual(report.metrics["significant_direction_agreement"], 0.5)
        self.assertLess(report.overall_concordance, 1.0)

    def test_engine_auto_selects_deg_before_generic_table(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.csv"
            right = Path(temp_dir) / "right.csv"
            left.write_text("symbol,logFC,FDR\nG1,1.0,0.01\nG2,-1.0,0.02\n", encoding="utf-8")
            right.write_text("symbol,logFC,FDR\nG1,1.1,0.01\nG2,-1.1,0.02\n", encoding="utf-8")

            report = ComparisonEngine().compare(str(left), str(right))

        self.assertEqual(report.comparator, "DEGComparator")
        self.assertIn("logfc_spearman", report.metrics)

    def test_uses_raw_pvalue_when_adjusted_pvalue_is_absent(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text("gene\tlogFC\tpvalue\nG1\t1.0\t0.01\n", encoding="utf-8")
            right.write_text("gene\tlogFC\tpvalue\nG1\t1.0\t0.01\n", encoding="utf-8")

            report = DEGComparator().compare(str(left), str(right))

        self.assertFalse(report.details["file_a_columns"]["uses_adjusted_pvalue"])
        self.assertTrue(report.warnings)


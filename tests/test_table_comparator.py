from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.comparators.table import TableComparator


class TableComparatorTests(TestCase):
    def test_identical_tables_are_fully_concordant(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            content = "gene_id\tlogFC\tstatus\nG1\t1.0\tup\nG2\t-2.0\tdown\n"
            left.write_text(content, encoding="utf-8")
            right.write_text(content, encoding="utf-8")

            report = TableComparator().compare(str(left), str(right), key_column="gene_id")

        self.assertEqual(report.overall_concordance, 1.0)
        self.assertEqual(report.metrics["row_overlap"], 1.0)
        self.assertEqual(report.metrics["numeric.logFC.similarity"], 1.0)
        self.assertEqual(report.metrics["categorical.status.agreement"], 1.0)

    def test_keyed_tables_report_row_and_value_differences(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text("gene_id\tlogFC\tstatus\nG1\t1.0\tup\nG2\t-2.0\tdown\nG3\t0.1\tns\n", encoding="utf-8")
            right.write_text("gene_id\tlogFC\tstatus\nG1\t1.5\tup\nG2\t-1.0\tup\nG4\t3.0\tup\n", encoding="utf-8")

            report = TableComparator().compare(str(left), str(right), key_column="gene_id")

        self.assertAlmostEqual(report.metrics["row_overlap"], 0.5)
        self.assertEqual(report.details["compared_rows"], 2)
        self.assertLess(report.overall_concordance, 1.0)
        self.assertIn("numeric.logFC.mae", report.metrics)
        self.assertIn("categorical.status.cohen_kappa", report.metrics)

    def test_auto_detects_key_column(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.csv"
            right = Path(temp_dir) / "right.csv"
            left.write_text("gene_id,value\nG1,1\nG2,2\n", encoding="utf-8")
            right.write_text("gene_id,value\nG2,2\nG1,1\n", encoding="utf-8")

            report = TableComparator().compare(str(left), str(right))

        self.assertEqual(report.details["key_column"], "gene_id")
        self.assertEqual(report.overall_concordance, 1.0)


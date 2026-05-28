from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.comparators.fasta import FASTAComparator
from biocompare.core.engine import ComparisonEngine


class FASTAComparatorTests(TestCase):
    def test_identical_fasta_files_are_fully_concordant(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.fa"
            right = Path(temp_dir) / "right.fa"
            content = ">seq1\nACGT\n>seq2\nGGGG\n"
            left.write_text(content, encoding="utf-8")
            right.write_text(content, encoding="utf-8")

            report = FASTAComparator().compare(str(left), str(right))

        self.assertEqual(report.overall_concordance, 1.0)
        self.assertEqual(report.metrics["record_overlap"], 1.0)
        self.assertEqual(report.metrics["exact_sequence_agreement"], 1.0)

    def test_fasta_metrics_capture_record_and_sequence_differences(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.fa"
            right = Path(temp_dir) / "right.fa"
            left.write_text(">seq1\nACGTACGT\n>seq2\nGGGG\n>seq3\nTTAA\n", encoding="utf-8")
            right.write_text(">seq1\nACGTACGT\n>seq2\nGGGA\n>seq4\nCCCC\n", encoding="utf-8")

            report = FASTAComparator().compare(str(left), str(right))

        self.assertAlmostEqual(report.metrics["record_overlap"], 2 / 4)
        self.assertAlmostEqual(report.metrics["exact_sequence_agreement"], 1 / 2)
        self.assertLess(report.overall_concordance, 1.0)
        self.assertEqual(report.details["shared_records"], 2)

    def test_engine_auto_selects_fasta_for_fasta_extension(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.fasta"
            right = Path(temp_dir) / "right.fasta"
            left.write_text(">seq1\nACGT\n", encoding="utf-8")
            right.write_text(">seq1\nACGA\n", encoding="utf-8")

            report = ComparisonEngine().compare(str(left), str(right))

        self.assertEqual(report.comparator, "FASTAComparator")
        self.assertIn("gc_content_similarity", report.metrics)

    def test_fastq_sequence_comparison_warns_quality_scores_not_compared(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.fastq"
            right = Path(temp_dir) / "right.fastq"
            left.write_text("@seq1\nACGT\n+\n!!!!\n", encoding="utf-8")
            right.write_text("@seq1\nACGT\n+\n####\n", encoding="utf-8")

            report = FASTAComparator().compare(str(left), str(right))

        self.assertEqual(report.overall_concordance, 1.0)
        self.assertTrue(any("quality" in warning for warning in report.warnings))


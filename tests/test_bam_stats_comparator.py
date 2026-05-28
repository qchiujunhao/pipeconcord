from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.comparators.bam_stats import BAMStatsComparator
from biocompare.core.engine import ComparisonEngine


class BAMStatsComparatorTests(TestCase):
    def test_identical_flagstat_outputs_are_fully_concordant(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.flagstat.txt"
            right = Path(temp_dir) / "right.flagstat.txt"
            content = (
                "1000 + 0 in total (QC-passed reads + QC-failed reads)\n"
                "900 + 0 mapped (90.00% : N/A)\n"
                "100 + 0 duplicates\n"
                "800 + 0 paired in sequencing\n"
                "700 + 0 properly paired (87.50% : N/A)\n"
            )
            left.write_text(content, encoding="utf-8")
            right.write_text(content, encoding="utf-8")

            report = BAMStatsComparator().compare(str(left), str(right))

        self.assertEqual(report.overall_concordance, 1.0)
        self.assertEqual(report.metrics["alignment_rate_similarity"], 1.0)
        self.assertEqual(report.metrics["total_reads_ratio"], 1.0)

    def test_flagstat_metrics_capture_alignment_summary_differences(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.txt"
            right = Path(temp_dir) / "right.txt"
            left.write_text(
                "1000 + 0 in total (QC-passed reads + QC-failed reads)\n"
                "900 + 0 mapped (90.00% : N/A)\n"
                "100 + 0 duplicates\n"
                "800 + 0 paired in sequencing\n"
                "700 + 0 properly paired (87.50% : N/A)\n",
                encoding="utf-8",
            )
            right.write_text(
                "1100 + 0 in total (QC-passed reads + QC-failed reads)\n"
                "990 + 0 mapped (90.00% : N/A)\n"
                "120 + 0 duplicates\n"
                "900 + 0 paired in sequencing\n"
                "760 + 0 properly paired (84.44% : N/A)\n",
                encoding="utf-8",
            )

            report = BAMStatsComparator().compare(str(left), str(right))

        self.assertAlmostEqual(report.metrics["alignment_rate_file_a"], 0.9)
        self.assertAlmostEqual(report.metrics["alignment_rate_file_b"], 0.9)
        self.assertLess(report.overall_concordance, 1.0)
        self.assertIn("proper_pair_rate_similarity", report.metrics)

    def test_samtools_stats_summary_lines_are_supported(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.stats"
            right = Path(temp_dir) / "right.stats"
            left.write_text(
                "SN\traw total sequences:\t1000\n"
                "SN\treads mapped:\t900\n"
                "SN\treads duplicated:\t100\n"
                "SN\tinsert size average:\t250.0\n"
                "SN\taverage length:\t150\n"
                "SN\terror rate:\t0.001\n",
                encoding="utf-8",
            )
            right.write_text(
                "SN\traw total sequences:\t1000\n"
                "SN\treads mapped:\t890\n"
                "SN\treads duplicated:\t110\n"
                "SN\tinsert size average:\t260.0\n"
                "SN\taverage length:\t150\n"
                "SN\terror rate:\t0.002\n",
                encoding="utf-8",
            )

            report = BAMStatsComparator().compare(str(left), str(right))

        self.assertIn("insert_size_average_ratio", report.metrics)
        self.assertIn("error_rate_similarity", report.metrics)

    def test_engine_auto_selects_bam_stats_for_flagstat_text(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.txt"
            right = Path(temp_dir) / "right.txt"
            content = "10 + 0 in total (QC-passed reads + QC-failed reads)\n9 + 0 mapped (90.00% : N/A)\n"
            left.write_text(content, encoding="utf-8")
            right.write_text(content, encoding="utf-8")

            report = ComparisonEngine().compare(str(left), str(right))

        self.assertEqual(report.comparator, "BAMStatsComparator")


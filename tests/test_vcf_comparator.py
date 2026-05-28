from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.comparators.vcf import VCFComparator
from biocompare.core.engine import ComparisonEngine


class VCFComparatorTests(TestCase):
    def test_identical_vcf_files_are_fully_concordant(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.vcf"
            right = Path(temp_dir) / "right.vcf"
            content = (
                "##fileformat=VCFv4.3\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n"
                "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t1/1\t0/1\n"
            )
            left.write_text(content, encoding="utf-8")
            right.write_text(content, encoding="utf-8")

            report = VCFComparator().compare(str(left), str(right))

        self.assertEqual(report.overall_concordance, 1.0)
        self.assertEqual(report.metrics["variant_jaccard"], 1.0)
        self.assertEqual(report.metrics["genotype_concordance"], 1.0)

    def test_vcf_metrics_capture_variant_and_genotype_differences(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.vcf"
            right = Path(temp_dir) / "right.vcf"
            left.write_text(
                "##fileformat=VCFv4.3\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n"
                "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t1/1\t0/1\n"
                "chr1\t300\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\t0/1\n",
                encoding="utf-8",
            )
            right.write_text(
                "##fileformat=VCFv4.3\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t0/0\n"
                "chr1\t200\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/1\n"
                "chr1\t400\t.\tT\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n",
                encoding="utf-8",
            )

            report = VCFComparator().compare(str(left), str(right))

        self.assertAlmostEqual(report.metrics["variant_jaccard"], 2 / 4)
        self.assertAlmostEqual(report.metrics["genotype_concordance"], 3 / 4)
        self.assertLess(report.overall_concordance, 1.0)
        self.assertEqual(report.details["genotype_calls_compared"], 4)

    def test_engine_auto_selects_vcf_for_vcf_extension(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.vcf"
            right = Path(temp_dir) / "right.vcf"
            left.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t.\tPASS\t.\n", encoding="utf-8")
            right.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t.\tPASS\t.\n", encoding="utf-8")

            report = ComparisonEngine().compare(str(left), str(right))

        self.assertEqual(report.comparator, "VCFComparator")
        self.assertIn("variant_jaccard", report.metrics)


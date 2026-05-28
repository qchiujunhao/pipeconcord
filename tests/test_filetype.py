from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.detection.filetype import detect_file_type, sniff_delimiter


class FileTypeTests(TestCase):
    def test_detects_tsv_by_extension(self):
        with TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "counts.tsv"
            path.write_text("gene\tsample_a\nG1\t10\n", encoding="utf-8")
            detected = detect_file_type(str(path))
        self.assertEqual(detected.kind, "tsv")
        self.assertEqual(detected.delimiter, "\t")

    def test_detects_vcf_by_header(self):
        with TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "calls.txt"
            path.write_text("##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", encoding="utf-8")
            detected = detect_file_type(str(path))
        self.assertEqual(detected.kind, "vcf")

    def test_sniffs_csv_delimiter(self):
        self.assertEqual(sniff_delimiter("id,value\na,1\nb,2\n"), ",")


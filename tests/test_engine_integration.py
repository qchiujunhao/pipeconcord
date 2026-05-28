from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from biocompare.core.engine import ComparisonEngine
from biocompare.io.report_writers import report_to_json


class EngineIntegrationTests(TestCase):
    def test_engine_finds_table_comparator(self):
        with TemporaryDirectory() as temp_dir:
            left = Path(temp_dir) / "left.tsv"
            right = Path(temp_dir) / "right.tsv"
            left.write_text("id\tvalue\nA\t1\nB\t2\n", encoding="utf-8")
            right.write_text("id\tvalue\nA\t1\nB\t3\n", encoding="utf-8")

            report = ComparisonEngine().compare(str(left), str(right), key_column="id")

        self.assertEqual(report.comparator, "TableComparator")
        self.assertIn("numeric.value.similarity", report.metrics)
        self.assertIn('"overall_concordance"', report_to_json(report))


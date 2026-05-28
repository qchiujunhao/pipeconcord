"""Semantic comparison of bioinformatics pipeline outputs."""

from biocompare._version import __version__
from biocompare.core.engine import ComparisonEngine
from biocompare.core.report import ConcordanceReport

__all__ = ["ComparisonEngine", "ConcordanceReport", "__version__"]


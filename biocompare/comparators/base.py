from __future__ import annotations

from abc import ABC, abstractmethod

from biocompare.core.report import ConcordanceReport


class Comparator(ABC):
    """Base class for semantic output comparators."""

    name = "base"
    supported_types: tuple[str, ...] = ()

    @abstractmethod
    def can_handle(self, file_a: str, file_b: str, **kwargs: object) -> bool:
        """Return True when this comparator can compare the two inputs."""

    @abstractmethod
    def compare(self, file_a: str, file_b: str, **kwargs: object) -> ConcordanceReport:
        """Run comparison and return a unified concordance report."""


from __future__ import annotations

from collections.abc import Iterable

from biocompare.comparators import register_builtin_comparators
from biocompare.core.registry import ComparatorRegistry
from biocompare.core.report import ConcordanceReport


class ComparisonEngine:
    """Orchestrates comparator selection and execution."""

    def __init__(self, registry: type[ComparatorRegistry] = ComparatorRegistry) -> None:
        self.registry = registry
        register_builtin_comparators(self.registry)

    def compare(
        self,
        file_a: str,
        file_b: str,
        *,
        file_type: str | None = None,
        **kwargs: object,
    ) -> ConcordanceReport:
        comparator = self.registry.find(file_a, file_b, file_type=file_type, **kwargs)
        return comparator.compare(file_a, file_b, file_type=file_type, **kwargs)

    def compare_many(
        self,
        pairs: Iterable[tuple[str, str]],
        *,
        file_type: str | None = None,
        **kwargs: object,
    ) -> list[ConcordanceReport]:
        return [self.compare(file_a, file_b, file_type=file_type, **kwargs) for file_a, file_b in pairs]


from __future__ import annotations

from biocompare.core.registry import ComparatorRegistry


def select_comparator(
    file_a: str,
    file_b: str,
    *,
    file_type: str | None = None,
    registry: type[ComparatorRegistry] = ComparatorRegistry,
    **kwargs: object,
):
    return registry.find(file_a, file_b, file_type=file_type, **kwargs)


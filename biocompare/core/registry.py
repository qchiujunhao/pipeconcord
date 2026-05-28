from __future__ import annotations

from importlib.metadata import entry_points
from typing import ClassVar

from biocompare.comparators.base import Comparator


class ComparatorRegistry:
    """Registry for built-in and third-party comparators."""

    _comparators: ClassVar[list[type[Comparator]]] = []
    _entry_points_loaded: ClassVar[bool] = False

    @classmethod
    def register(cls, comparator_cls: type[Comparator]) -> type[Comparator]:
        if not issubclass(comparator_cls, Comparator):
            raise TypeError("registered comparators must subclass Comparator")
        if comparator_cls not in cls._comparators:
            cls._comparators.append(comparator_cls)
        return comparator_cls

    @classmethod
    def clear(cls) -> None:
        cls._comparators = []
        cls._entry_points_loaded = False

    @classmethod
    def comparators(cls) -> list[type[Comparator]]:
        return list(cls._comparators)

    @classmethod
    def load_entry_points(cls) -> None:
        if cls._entry_points_loaded:
            return
        discovered = entry_points()
        if hasattr(discovered, "select"):
            group = discovered.select(group="biocompare.comparators")
        else:
            group = discovered.get("biocompare.comparators", [])

        for entry_point in group:
            comparator_cls = entry_point.load()
            cls.register(comparator_cls)
        cls._entry_points_loaded = True

    @classmethod
    def find(
        cls,
        file_a: str,
        file_b: str,
        *,
        file_type: str | None = None,
        **kwargs: object,
    ) -> Comparator:
        cls.load_entry_points()
        for comparator_cls in cls._comparators:
            comparator = comparator_cls()
            if comparator.can_handle(file_a, file_b, file_type=file_type, **kwargs):
                return comparator
        label = file_type or "auto-detected type"
        raise ValueError(f"No comparator found for {file_a!r}, {file_b!r} using {label}")


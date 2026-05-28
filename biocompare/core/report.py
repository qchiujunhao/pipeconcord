from __future__ import annotations

from dataclasses import asdict, dataclass, field
from math import isfinite
from typing import Any


@dataclass(slots=True)
class ConcordanceReport:
    """Unified output from a comparator."""

    comparator: str
    file_a: str
    file_b: str
    overall_concordance: float
    metrics: dict[str, float]
    details: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.overall_concordance = float(self.overall_concordance)
        if not isfinite(self.overall_concordance):
            raise ValueError("overall_concordance must be finite")
        if not 0.0 <= self.overall_concordance <= 1.0:
            raise ValueError("overall_concordance must be between 0.0 and 1.0")

        cleaned_metrics: dict[str, float] = {}
        for name, value in self.metrics.items():
            metric_value = float(value)
            if not isfinite(metric_value):
                raise ValueError(f"metric {name!r} must be finite")
            cleaned_metrics[str(name)] = metric_value
        self.metrics = cleaned_metrics
        self.details = dict(self.details)
        self.warnings = list(self.warnings)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


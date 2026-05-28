from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from biocompare.comparators.base import Comparator
from biocompare.core.report import ConcordanceReport
from biocompare.core.utils import clamp01


FLAGSTAT_COUNT_RE = re.compile(r"^(\d+)\s+\+\s+\d+\s+(.+)$")


@dataclass(slots=True)
class AlignmentStats:
    path: str
    total_reads: float | None = None
    mapped_reads: float | None = None
    duplicate_reads: float | None = None
    paired_reads: float | None = None
    properly_paired_reads: float | None = None
    insert_size_average: float | None = None
    average_length: float | None = None
    error_rate: float | None = None
    recognized_metrics: int = 0

    @property
    def alignment_rate(self) -> float | None:
        return fraction(self.mapped_reads, self.total_reads)

    @property
    def duplicate_rate(self) -> float | None:
        return fraction(self.duplicate_reads, self.total_reads)

    @property
    def proper_pair_rate(self) -> float | None:
        return fraction(self.properly_paired_reads, self.paired_reads)


class BAMStatsComparator(Comparator):
    """Comparator for samtools flagstat/stats alignment summaries."""

    name = "bam_stats"
    supported_types = ("bam_stats", "bam-stats", "flagstat", "samtools-stats")

    def can_handle(self, file_a: str, file_b: str, **kwargs: object) -> bool:
        requested_type = kwargs.get("file_type")
        if requested_type in self.supported_types:
            return True
        if requested_type is not None:
            return False
        return looks_like_bam_stats(file_a) and looks_like_bam_stats(file_b)

    def compare(self, file_a: str, file_b: str, **kwargs: object) -> ConcordanceReport:
        stats_a = parse_alignment_stats(file_a)
        stats_b = parse_alignment_stats(file_b)
        metrics: dict[str, float] = {}
        scores: list[float] = []

        add_ratio_metric(metrics, scores, "total_reads_ratio", stats_a.total_reads, stats_b.total_reads)
        add_ratio_metric(metrics, scores, "mapped_reads_ratio", stats_a.mapped_reads, stats_b.mapped_reads)
        add_rate_similarity(metrics, scores, "alignment_rate_similarity", stats_a.alignment_rate, stats_b.alignment_rate)
        add_rate_similarity(metrics, scores, "duplicate_rate_similarity", stats_a.duplicate_rate, stats_b.duplicate_rate)
        add_rate_similarity(metrics, scores, "proper_pair_rate_similarity", stats_a.proper_pair_rate, stats_b.proper_pair_rate)
        add_ratio_metric(metrics, scores, "insert_size_average_ratio", stats_a.insert_size_average, stats_b.insert_size_average)
        add_ratio_metric(metrics, scores, "average_length_ratio", stats_a.average_length, stats_b.average_length)
        add_rate_similarity(metrics, scores, "error_rate_similarity", stats_a.error_rate, stats_b.error_rate)

        if stats_a.alignment_rate is not None:
            metrics["alignment_rate_file_a"] = stats_a.alignment_rate
        if stats_b.alignment_rate is not None:
            metrics["alignment_rate_file_b"] = stats_b.alignment_rate
        if stats_a.duplicate_rate is not None:
            metrics["duplicate_rate_file_a"] = stats_a.duplicate_rate
        if stats_b.duplicate_rate is not None:
            metrics["duplicate_rate_file_b"] = stats_b.duplicate_rate

        details = {
            "file_a": stats_to_dict(stats_a),
            "file_b": stats_to_dict(stats_b),
            "scored_metric_count": len(scores),
        }
        return ConcordanceReport(
            comparator=self.__class__.__name__,
            file_a=str(file_a),
            file_b=str(file_b),
            overall_concordance=clamp01(mean(scores)),
            metrics=metrics,
            details=details,
            warnings=[],
        )


def looks_like_bam_stats(path: str) -> bool:
    text = Path(path).read_text(encoding="utf-8", errors="replace")[:8192]
    if text.startswith("SN\t") or "\nSN\t" in text:
        return True
    lowered = text.lower()
    return " in total " in lowered and " mapped (" in lowered


def parse_alignment_stats(path: str) -> AlignmentStats:
    stats = AlignmentStats(path=str(path))
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("SN\t"):
            parse_samtools_stats_line(stats, stripped)
        else:
            parse_flagstat_line(stats, stripped)
    if stats.recognized_metrics == 0:
        raise ValueError(f"{path!r} does not look like samtools stats or flagstat output")
    return stats


def parse_samtools_stats_line(stats: AlignmentStats, line: str) -> None:
    fields = line.split("\t")
    if len(fields) < 3:
        return
    key = fields[1].rstrip(":").lower()
    value = parse_number(fields[2])
    if value is None:
        return
    mapping = {
        "raw total sequences": "total_reads",
        "reads mapped": "mapped_reads",
        "reads duplicated": "duplicate_reads",
        "insert size average": "insert_size_average",
        "average length": "average_length",
        "error rate": "error_rate",
    }
    attribute = mapping.get(key)
    if attribute is not None:
        setattr(stats, attribute, value)
        stats.recognized_metrics += 1


def parse_flagstat_line(stats: AlignmentStats, line: str) -> None:
    match = FLAGSTAT_COUNT_RE.match(line)
    if match is None:
        return
    count = float(match.group(1))
    label = match.group(2).lower()
    if " in total " in f" {label} ":
        stats.total_reads = count
    elif label.startswith("mapped "):
        stats.mapped_reads = count
    elif label.startswith("duplicates"):
        stats.duplicate_reads = count
    elif label.startswith("paired in sequencing"):
        stats.paired_reads = count
    elif label.startswith("properly paired"):
        stats.properly_paired_reads = count
    else:
        return
    stats.recognized_metrics += 1


def parse_number(value: str) -> float | None:
    try:
        return float(value.strip())
    except ValueError:
        return None


def fraction(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator is None or denominator == 0:
        return None
    return numerator / denominator


def add_ratio_metric(
    metrics: dict[str, float],
    scores: list[float],
    name: str,
    left: float | None,
    right: float | None,
) -> None:
    if left is None and right is None:
        return
    score = magnitude_ratio(left or 0.0, right or 0.0)
    metrics[name] = score
    scores.append(score)


def add_rate_similarity(
    metrics: dict[str, float],
    scores: list[float],
    name: str,
    left: float | None,
    right: float | None,
) -> None:
    if left is None and right is None:
        return
    score = 0.0 if left is None or right is None else clamp01(1.0 - abs(left - right))
    metrics[name] = score
    scores.append(score)


def magnitude_ratio(left: float, right: float) -> float:
    if left == 0 and right == 0:
        return 1.0
    return min(abs(left), abs(right)) / max(abs(left), abs(right))


def mean(values: list[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def stats_to_dict(stats: AlignmentStats) -> dict[str, float | str | None]:
    return {
        "path": stats.path,
        "total_reads": stats.total_reads,
        "mapped_reads": stats.mapped_reads,
        "duplicate_reads": stats.duplicate_reads,
        "paired_reads": stats.paired_reads,
        "properly_paired_reads": stats.properly_paired_reads,
        "insert_size_average": stats.insert_size_average,
        "average_length": stats.average_length,
        "error_rate": stats.error_rate,
        "alignment_rate": stats.alignment_rate,
        "duplicate_rate": stats.duplicate_rate,
        "proper_pair_rate": stats.proper_pair_rate,
    }


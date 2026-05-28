from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from biocompare.comparators.base import Comparator
from biocompare.core.report import ConcordanceReport
from biocompare.core.utils import clamp01
from biocompare.detection.filetype import detect_file_type


@dataclass(frozen=True, slots=True)
class Interval:
    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(slots=True)
class BedData:
    path: str
    intervals: list[Interval]
    merged_by_chrom: dict[str, list[Interval]]


class BEDComparator(Comparator):
    """Comparator for BED-like genomic interval files."""

    name = "bed"
    supported_types = ("bed",)

    def can_handle(self, file_a: str, file_b: str, **kwargs: object) -> bool:
        requested_type = kwargs.get("file_type")
        if requested_type == "bed":
            return True
        if requested_type is not None:
            return False
        return detect_file_type(file_a).kind == "bed" and detect_file_type(file_b).kind == "bed"

    def compare(self, file_a: str, file_b: str, **kwargs: object) -> ConcordanceReport:
        min_reciprocal_overlap = optional_float(kwargs.get("min_reciprocal_overlap"), 0.0, "min_reciprocal_overlap")
        if min_reciprocal_overlap > 1.0:
            raise ValueError("min_reciprocal_overlap must be between 0.0 and 1.0")

        bed_a = read_bed(file_a)
        bed_b = read_bed(file_b)
        warnings: list[str] = []
        if not bed_a.intervals:
            warnings.append("file_a contains no intervals.")
        if not bed_b.intervals:
            warnings.append("file_b contains no intervals.")

        chromosomes_a = set(bed_a.merged_by_chrom)
        chromosomes_b = set(bed_b.merged_by_chrom)
        shared_chromosomes = sorted(chromosomes_a & chromosomes_b)
        total_bp_a = total_bp(bed_a.merged_by_chrom)
        total_bp_b = total_bp(bed_b.merged_by_chrom)
        intersection_bp = intersect_bp(bed_a.merged_by_chrom, bed_b.merged_by_chrom)
        union_bp = total_bp_a + total_bp_b - intersection_bp
        bp_jaccard = intersection_bp / union_bp if union_bp else 1.0
        file_a_coverage = intersection_bp / total_bp_a if total_bp_a else 1.0 if total_bp_b == 0 else 0.0
        file_b_coverage = intersection_bp / total_bp_b if total_bp_b else 1.0 if total_bp_a == 0 else 0.0

        matched_a = matched_interval_count(bed_a.intervals, bed_b.intervals, min_reciprocal_overlap)
        matched_b = matched_interval_count(bed_b.intervals, bed_a.intervals, min_reciprocal_overlap)
        interval_recall = matched_a / len(bed_a.intervals) if bed_a.intervals else 1.0 if not bed_b.intervals else 0.0
        interval_precision = matched_b / len(bed_b.intervals) if bed_b.intervals else 1.0 if not bed_a.intervals else 0.0
        interval_f1 = f1(interval_precision, interval_recall)

        lengths_a = [interval.length for interval in bed_a.intervals]
        lengths_b = [interval.length for interval in bed_b.intervals]
        mean_length_ratio = magnitude_ratio(mean(lengths_a), mean(lengths_b))
        median_length_ratio = magnitude_ratio(median(lengths_a), median(lengths_b))
        interval_count_ratio = magnitude_ratio(len(bed_a.intervals), len(bed_b.intervals))
        merged_interval_count_ratio = magnitude_ratio(merged_interval_count(bed_a.merged_by_chrom), merged_interval_count(bed_b.merged_by_chrom))

        metrics = {
            "bp_jaccard": bp_jaccard,
            "file_a_bp_coverage": file_a_coverage,
            "file_b_bp_coverage": file_b_coverage,
            "interval_precision_file_b_vs_a": interval_precision,
            "interval_recall_file_b_vs_a": interval_recall,
            "interval_f1_file_b_vs_a": interval_f1,
            "mean_length_ratio": mean_length_ratio,
            "median_length_ratio": median_length_ratio,
            "interval_count_ratio": interval_count_ratio,
            "merged_interval_count_ratio": merged_interval_count_ratio,
        }
        overall = clamp01(
            0.45 * bp_jaccard
            + 0.35 * interval_f1
            + 0.10 * mean_length_ratio
            + 0.10 * interval_count_ratio
        )
        details = {
            "coordinate_system": "BED 0-based half-open",
            "min_reciprocal_overlap": min_reciprocal_overlap,
            "file_a_intervals": len(bed_a.intervals),
            "file_b_intervals": len(bed_b.intervals),
            "file_a_merged_intervals": merged_interval_count(bed_a.merged_by_chrom),
            "file_b_merged_intervals": merged_interval_count(bed_b.merged_by_chrom),
            "file_a_total_bp": total_bp_a,
            "file_b_total_bp": total_bp_b,
            "intersection_bp": intersection_bp,
            "union_bp": union_bp,
            "file_a_chromosomes": sorted(chromosomes_a),
            "file_b_chromosomes": sorted(chromosomes_b),
            "shared_chromosomes": shared_chromosomes,
            "file_a_only_chromosomes": sorted(chromosomes_a - chromosomes_b),
            "file_b_only_chromosomes": sorted(chromosomes_b - chromosomes_a),
        }
        return ConcordanceReport(
            comparator=self.__class__.__name__,
            file_a=str(file_a),
            file_b=str(file_b),
            overall_concordance=overall,
            metrics=metrics,
            details=details,
            warnings=warnings,
        )


def read_bed(path: str) -> BedData:
    intervals: list[Interval] = []
    for line_number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("track") or stripped.startswith("browser"):
            continue
        fields = stripped.split()
        if len(fields) < 3:
            raise ValueError(f"{path!r} line {line_number} has fewer than 3 BED columns")
        chrom = fields[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError as exc:
            raise ValueError(f"{path!r} line {line_number} has non-integer coordinates") from exc
        if start < 0:
            raise ValueError(f"{path!r} line {line_number} has a negative start coordinate")
        if end <= start:
            raise ValueError(f"{path!r} line {line_number} has end <= start")
        intervals.append(Interval(chrom=chrom, start=start, end=end))
    return BedData(path=str(path), intervals=intervals, merged_by_chrom=merge_by_chrom(intervals))


def merge_by_chrom(intervals: list[Interval]) -> dict[str, list[Interval]]:
    by_chrom: dict[str, list[Interval]] = {}
    for interval in intervals:
        by_chrom.setdefault(interval.chrom, []).append(interval)
    merged: dict[str, list[Interval]] = {}
    for chrom, chrom_intervals in by_chrom.items():
        sorted_intervals = sorted(chrom_intervals, key=lambda interval: (interval.start, interval.end))
        merged_chrom: list[Interval] = []
        for interval in sorted_intervals:
            if not merged_chrom or interval.start > merged_chrom[-1].end:
                merged_chrom.append(interval)
            else:
                previous = merged_chrom[-1]
                merged_chrom[-1] = Interval(chrom=chrom, start=previous.start, end=max(previous.end, interval.end))
        merged[chrom] = merged_chrom
    return merged


def total_bp(intervals_by_chrom: dict[str, list[Interval]]) -> int:
    return sum(interval.length for intervals in intervals_by_chrom.values() for interval in intervals)


def intersect_bp(left: dict[str, list[Interval]], right: dict[str, list[Interval]]) -> int:
    total = 0
    for chrom in set(left) & set(right):
        total += intersect_sorted_intervals(left[chrom], right[chrom])
    return total


def intersect_sorted_intervals(left: list[Interval], right: list[Interval]) -> int:
    total = 0
    left_index = 0
    right_index = 0
    while left_index < len(left) and right_index < len(right):
        left_interval = left[left_index]
        right_interval = right[right_index]
        overlap_start = max(left_interval.start, right_interval.start)
        overlap_end = min(left_interval.end, right_interval.end)
        if overlap_end > overlap_start:
            total += overlap_end - overlap_start
        if left_interval.end < right_interval.end:
            left_index += 1
        else:
            right_index += 1
    return total


def matched_interval_count(query: list[Interval], target: list[Interval], min_reciprocal_overlap: float) -> int:
    target_by_chrom: dict[str, list[Interval]] = {}
    for interval in target:
        target_by_chrom.setdefault(interval.chrom, []).append(interval)
    for chrom in target_by_chrom:
        target_by_chrom[chrom].sort(key=lambda interval: interval.start)

    matched = 0
    for interval in query:
        if interval_matches(interval, target_by_chrom.get(interval.chrom, []), min_reciprocal_overlap):
            matched += 1
    return matched


def interval_matches(query: Interval, targets: list[Interval], min_reciprocal_overlap: float) -> bool:
    for target in targets:
        if target.end <= query.start:
            continue
        if target.start >= query.end:
            break
        overlap = min(query.end, target.end) - max(query.start, target.start)
        if overlap <= 0:
            continue
        if min_reciprocal_overlap <= 0:
            return True
        if overlap / query.length >= min_reciprocal_overlap and overlap / target.length >= min_reciprocal_overlap:
            return True
    return False


def merged_interval_count(intervals_by_chrom: dict[str, list[Interval]]) -> int:
    return sum(len(intervals) for intervals in intervals_by_chrom.values())


def magnitude_ratio(left: float, right: float) -> float:
    if left == 0 and right == 0:
        return 1.0
    return min(abs(left), abs(right)) / max(abs(left), abs(right))


def mean(values: list[int]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)


def median(values: list[int]) -> float:
    if not values:
        return 0.0
    sorted_values = sorted(values)
    middle = len(sorted_values) // 2
    if len(sorted_values) % 2:
        return float(sorted_values[middle])
    return (sorted_values[middle - 1] + sorted_values[middle]) / 2.0


def f1(precision: float, recall: float) -> float:
    if precision + recall == 0:
        return 0.0
    return 2 * precision * recall / (precision + recall)


def optional_float(value: object, default: float, name: str) -> float:
    if value is None:
        return default
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be numeric") from exc
    if parsed < 0:
        raise ValueError(f"{name} must be non-negative")
    return parsed


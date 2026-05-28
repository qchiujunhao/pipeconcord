from __future__ import annotations

from collections import Counter
from math import sqrt
from typing import Iterable, Sequence, TypeVar

T = TypeVar("T")


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, value))


def safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    if denominator == 0:
        return default
    return numerator / denominator


def jaccard(left: Iterable[T], right: Iterable[T]) -> float:
    left_set = set(left)
    right_set = set(right)
    union = left_set | right_set
    if not union:
        return 1.0
    return len(left_set & right_set) / len(union)


def agreement_rate(left: Sequence[object], right: Sequence[object]) -> float:
    if len(left) != len(right):
        raise ValueError("agreement_rate requires sequences of equal length")
    if not left:
        return 0.0
    return sum(a == b for a, b in zip(left, right)) / len(left)


def mean_absolute_error(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right):
        raise ValueError("mean_absolute_error requires sequences of equal length")
    if not left:
        return 0.0
    return sum(abs(a - b) for a, b in zip(left, right)) / len(left)


def pearson(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right):
        raise ValueError("pearson requires sequences of equal length")
    if len(left) < 2:
        return 1.0 if left == right else 0.0

    mean_left = sum(left) / len(left)
    mean_right = sum(right) / len(right)
    centered_left = [value - mean_left for value in left]
    centered_right = [value - mean_right for value in right]
    numerator = sum(a * b for a, b in zip(centered_left, centered_right))
    denominator = sqrt(sum(a * a for a in centered_left) * sum(b * b for b in centered_right))
    if denominator == 0:
        return 1.0 if list(left) == list(right) else 0.0
    return numerator / denominator


def rankdata(values: Sequence[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda pair: pair[1])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(indexed):
        end = index
        while end + 1 < len(indexed) and indexed[end + 1][1] == indexed[index][1]:
            end += 1
        average_rank = (index + end + 2) / 2.0
        for ranked_index in range(index, end + 1):
            original_index = indexed[ranked_index][0]
            ranks[original_index] = average_rank
        index = end + 1
    return ranks


def spearman(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right):
        raise ValueError("spearman requires sequences of equal length")
    if len(left) < 2:
        return 1.0 if left == right else 0.0
    return pearson(rankdata(left), rankdata(right))


def cohen_kappa(left: Sequence[object], right: Sequence[object]) -> float:
    if len(left) != len(right):
        raise ValueError("cohen_kappa requires sequences of equal length")
    if not left:
        return 0.0

    observed = agreement_rate(left, right)
    left_counts = Counter(left)
    right_counts = Counter(right)
    total = len(left)
    labels = set(left_counts) | set(right_counts)
    expected = sum((left_counts[label] / total) * (right_counts[label] / total) for label in labels)
    if expected == 1.0:
        return 1.0 if observed == 1.0 else 0.0
    return (observed - expected) / (1.0 - expected)


def numeric_similarity(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right):
        raise ValueError("numeric_similarity requires sequences of equal length")
    if not left:
        return 0.0
    combined = list(left) + list(right)
    observed_range = max(combined) - min(combined)
    magnitude = max(abs(value) for value in combined)
    scale = max(observed_range, magnitude, 1.0)
    return clamp01(1.0 - mean_absolute_error(left, right) / scale)


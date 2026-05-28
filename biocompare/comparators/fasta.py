from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from biocompare.comparators.base import Comparator
from biocompare.core.report import ConcordanceReport
from biocompare.core.utils import clamp01, jaccard, mean_absolute_error
from biocompare.detection.filetype import detect_file_type


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    identifier: str
    sequence: str
    quality: str | None = None


@dataclass(slots=True)
class SequenceData:
    path: str
    file_kind: str
    records: dict[str, SequenceRecord]
    duplicate_ids: int


class FASTAComparator(Comparator):
    """Comparator for FASTA and FASTQ sequence files."""

    name = "fasta"
    supported_types = ("fasta", "fastq", "sequence")

    def can_handle(self, file_a: str, file_b: str, **kwargs: object) -> bool:
        requested_type = kwargs.get("file_type")
        if requested_type in self.supported_types:
            return True
        if requested_type is not None:
            return False
        return detect_file_type(file_a).kind in {"fasta", "fastq"} and detect_file_type(file_b).kind in {"fasta", "fastq"}

    def compare(self, file_a: str, file_b: str, **kwargs: object) -> ConcordanceReport:
        sequences_a = read_sequences(file_a)
        sequences_b = read_sequences(file_b)
        warnings: list[str] = []
        if sequences_a.duplicate_ids:
            warnings.append(f"file_a contains {sequences_a.duplicate_ids} duplicate sequence identifiers; kept the first occurrence.")
        if sequences_b.duplicate_ids:
            warnings.append(f"file_b contains {sequences_b.duplicate_ids} duplicate sequence identifiers; kept the first occurrence.")
        if sequences_a.file_kind == "fastq" or sequences_b.file_kind == "fastq":
            warnings.append("FASTQ quality scores are not compared yet; sequence-level metrics only.")

        ids_a = set(sequences_a.records)
        ids_b = set(sequences_b.records)
        shared_ids = sorted(ids_a & ids_b)
        record_overlap = jaccard(ids_a, ids_b)
        exact_matches = sum(
            sequences_a.records[identifier].sequence == sequences_b.records[identifier].sequence
            for identifier in shared_ids
        )
        exact_sequence_agreement = exact_matches / len(shared_ids) if shared_ids else 1.0 if not ids_a and not ids_b else 0.0

        lengths_a = [len(record.sequence) for record in sequences_a.records.values()]
        lengths_b = [len(record.sequence) for record in sequences_b.records.values()]
        shared_lengths_a = [len(sequences_a.records[identifier].sequence) for identifier in shared_ids]
        shared_lengths_b = [len(sequences_b.records[identifier].sequence) for identifier in shared_ids]
        length_mae = mean_absolute_error(shared_lengths_a, shared_lengths_b) if shared_lengths_a else 0.0
        shared_length_similarity = length_similarity(shared_lengths_a, shared_lengths_b)
        mean_length_ratio = magnitude_ratio(mean(lengths_a), mean(lengths_b))
        median_length_ratio = magnitude_ratio(median(lengths_a), median(lengths_b))
        total_bases_ratio = magnitude_ratio(sum(lengths_a), sum(lengths_b))
        gc_a = gc_fraction(record.sequence for record in sequences_a.records.values())
        gc_b = gc_fraction(record.sequence for record in sequences_b.records.values())
        gc_similarity = clamp01(1.0 - abs(gc_a - gc_b))
        record_count_ratio = magnitude_ratio(len(sequences_a.records), len(sequences_b.records))

        metrics = {
            "record_overlap": record_overlap,
            "exact_sequence_agreement": exact_sequence_agreement,
            "shared_length_similarity": shared_length_similarity,
            "shared_length_mae": length_mae,
            "mean_length_ratio": mean_length_ratio,
            "median_length_ratio": median_length_ratio,
            "total_bases_ratio": total_bases_ratio,
            "gc_content_file_a": gc_a,
            "gc_content_file_b": gc_b,
            "gc_content_similarity": gc_similarity,
            "record_count_ratio": record_count_ratio,
        }
        overall = clamp01(
            0.30 * record_overlap
            + 0.30 * exact_sequence_agreement
            + 0.15 * shared_length_similarity
            + 0.15 * gc_similarity
            + 0.10 * total_bases_ratio
        )
        details = {
            "file_a_kind": sequences_a.file_kind,
            "file_b_kind": sequences_b.file_kind,
            "file_a_records": len(sequences_a.records),
            "file_b_records": len(sequences_b.records),
            "shared_records": len(shared_ids),
            "file_a_only_records": len(ids_a - ids_b),
            "file_b_only_records": len(ids_b - ids_a),
            "file_a_total_bases": sum(lengths_a),
            "file_b_total_bases": sum(lengths_b),
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


def read_sequences(path: str) -> SequenceData:
    text = Path(path).read_text(encoding="utf-8")
    stripped = text.lstrip()
    if not stripped:
        return SequenceData(path=str(path), file_kind="fasta", records={}, duplicate_ids=0)
    if stripped.startswith(">"):
        records, duplicates = parse_fasta(text, path)
        return SequenceData(path=str(path), file_kind="fasta", records=records, duplicate_ids=duplicates)
    if stripped.startswith("@"):
        records, duplicates = parse_fastq(text, path)
        return SequenceData(path=str(path), file_kind="fastq", records=records, duplicate_ids=duplicates)
    raise ValueError(f"{path!r} is not recognizable as FASTA or FASTQ")


def parse_fasta(text: str, path: str) -> tuple[dict[str, SequenceRecord], int]:
    records: dict[str, SequenceRecord] = {}
    duplicates = 0
    current_identifier: str | None = None
    current_sequence: list[str] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            duplicates += store_fasta_record(records, current_identifier, current_sequence)
            current_identifier = stripped[1:].split()[0] if stripped[1:].split() else ""
            if not current_identifier:
                raise ValueError(f"{path!r} line {line_number} has an empty FASTA identifier")
            current_sequence = []
        else:
            if current_identifier is None:
                raise ValueError(f"{path!r} line {line_number} has sequence before a FASTA header")
            current_sequence.append(stripped.upper())
    duplicates += store_fasta_record(records, current_identifier, current_sequence)
    return records, duplicates


def store_fasta_record(records: dict[str, SequenceRecord], identifier: str | None, sequence_parts: list[str]) -> int:
    if identifier is None:
        return 0
    sequence = "".join(sequence_parts)
    if identifier in records:
        return 1
    records[identifier] = SequenceRecord(identifier=identifier, sequence=sequence)
    return 0


def parse_fastq(text: str, path: str) -> tuple[dict[str, SequenceRecord], int]:
    lines = [line.rstrip("\n") for line in text.splitlines() if line.strip()]
    if len(lines) % 4 != 0:
        raise ValueError(f"{path!r} does not contain complete 4-line FASTQ records")
    records: dict[str, SequenceRecord] = {}
    duplicates = 0
    for index in range(0, len(lines), 4):
        header, sequence, separator, quality = lines[index : index + 4]
        if not header.startswith("@"):
            raise ValueError(f"{path!r} FASTQ record starting line {index + 1} is missing '@'")
        if not separator.startswith("+"):
            raise ValueError(f"{path!r} FASTQ record starting line {index + 1} is missing '+'")
        identifier = header[1:].split()[0] if header[1:].split() else ""
        if not identifier:
            raise ValueError(f"{path!r} FASTQ record starting line {index + 1} has an empty identifier")
        if len(sequence.strip()) != len(quality.strip()):
            raise ValueError(f"{path!r} FASTQ record {identifier!r} has sequence/quality length mismatch")
        if identifier in records:
            duplicates += 1
            continue
        records[identifier] = SequenceRecord(identifier=identifier, sequence=sequence.strip().upper(), quality=quality.strip())
    return records, duplicates


def length_similarity(left: list[int], right: list[int]) -> float:
    if not left and not right:
        return 1.0
    if not left or not right:
        return 0.0
    scale = max(max(left), max(right), 1)
    return clamp01(1.0 - mean_absolute_error(left, right) / scale)


def gc_fraction(sequences) -> float:
    gc_count = 0
    base_count = 0
    for sequence in sequences:
        for base in sequence.upper():
            if base in {"A", "C", "G", "T", "N"}:
                base_count += 1
                if base in {"G", "C"}:
                    gc_count += 1
    if base_count == 0:
        return 0.0
    return gc_count / base_count


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


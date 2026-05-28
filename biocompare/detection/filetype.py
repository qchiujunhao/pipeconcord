from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class FileType:
    kind: str
    confidence: float
    delimiter: str | None = None
    reason: str = ""


def detect_file_type(path: str) -> FileType:
    file_path = Path(path)
    name = file_path.name.lower()
    text_sample = read_text_sample(file_path)
    first_line = first_nonempty_line(text_sample)
    normalized_first_line = first_line.lower()

    if normalized_first_line.startswith("##fileformat=vcf") or first_line.startswith("#CHROM"):
        return FileType("vcf", 0.99, reason="VCF header")
    if name.endswith((".vcf", ".vcf.gz")):
        return FileType("vcf", 0.9, reason="VCF extension")
    if name.endswith((".bed", ".bed.gz")):
        return FileType("bed", 0.9, delimiter="\t", reason="BED extension")
    if name.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
        return FileType("fasta", 0.9, reason="FASTA extension")
    if name.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")):
        return FileType("fastq", 0.9, reason="FASTQ extension")
    if name.endswith((".csv", ".csv.gz")):
        return FileType("csv", 0.95, delimiter=",", reason="CSV extension")
    if name.endswith((".tsv", ".tab", ".tsv.gz", ".tab.gz")):
        return FileType("tsv", 0.95, delimiter="\t", reason="TSV extension")

    delimiter = sniff_delimiter(text_sample, filename=str(path))
    if delimiter is not None:
        return FileType("table", 0.7, delimiter=delimiter, reason="consistent delimited text")
    return FileType("unknown", 0.0, reason="no known extension or recognizable header")


def read_text_sample(path: Path, size: int = 65536) -> str:
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    try:
        with opener(path, "rb") as handle:
            sample = handle.read(size)
    except FileNotFoundError:
        raise
    except OSError:
        return ""
    return sample.decode("utf-8", errors="replace")


def first_nonempty_line(text: str) -> str:
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


def sniff_delimiter(text: str, filename: str | None = None) -> str | None:
    suffix = Path(filename).name.lower() if filename else ""
    if suffix.endswith((".csv", ".csv.gz")):
        return ","
    if suffix.endswith((".tsv", ".tab", ".tsv.gz", ".tab.gz", ".bed", ".bed.gz")):
        return "\t"

    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 1:
        return None
    candidates = ["\t", ",", ";"]
    best_delimiter: str | None = None
    best_score = 0
    for delimiter in candidates:
        counts = [line.count(delimiter) for line in lines[:10]]
        positive_counts = [count for count in counts if count > 0]
        if not positive_counts:
            continue
        score = min(positive_counts) * len(positive_counts)
        if len(set(positive_counts)) == 1:
            score += 10
        if score > best_score:
            best_score = score
            best_delimiter = delimiter
    return best_delimiter

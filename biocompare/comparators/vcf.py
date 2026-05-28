from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

from biocompare.comparators.base import Comparator
from biocompare.core.report import ConcordanceReport
from biocompare.core.utils import clamp01, jaccard, pearson, spearman
from biocompare.detection.filetype import detect_file_type


GENOTYPE_SPLIT_RE = re.compile(r"[\/|]")
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


@dataclass(frozen=True, slots=True)
class VariantKey:
    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass(slots=True)
class VariantRecord:
    key: VariantKey
    identifier: str
    qual: str
    filt: str
    info: str
    format_keys: list[str]
    samples: dict[str, dict[str, str]]


@dataclass(slots=True)
class VariantObservation:
    key: VariantKey
    record: VariantRecord
    alt_index: int


@dataclass(slots=True)
class VCFData:
    path: str
    samples: list[str]
    records: dict[VariantKey, VariantRecord]
    observations: dict[VariantKey, VariantObservation]
    duplicate_variants: int
    duplicate_observations: int


class VCFComparator(Comparator):
    """Lightweight comparator for VCF variant calls."""

    name = "vcf"
    supported_types = ("vcf",)

    def can_handle(self, file_a: str, file_b: str, **kwargs: object) -> bool:
        requested_type = kwargs.get("file_type")
        if requested_type == "vcf":
            return True
        if requested_type is not None:
            return False
        return detect_file_type(file_a).kind == "vcf" and detect_file_type(file_b).kind == "vcf"

    def compare(self, file_a: str, file_b: str, **kwargs: object) -> ConcordanceReport:
        vcf_a = read_vcf(file_a)
        vcf_b = read_vcf(file_b)
        warnings = [
            "VCF variants are split by ALT allele and minimally normalized by trimming shared prefix/suffix bases; reference-based left alignment is not performed."
        ]
        if vcf_a.duplicate_variants:
            warnings.append(f"file_a contains {vcf_a.duplicate_variants} duplicate variant keys; kept the first occurrence.")
        if vcf_b.duplicate_variants:
            warnings.append(f"file_b contains {vcf_b.duplicate_variants} duplicate variant keys; kept the first occurrence.")
        if vcf_a.duplicate_observations:
            warnings.append(f"file_a contains {vcf_a.duplicate_observations} duplicate normalized allele observations; kept the first occurrence.")
        if vcf_b.duplicate_observations:
            warnings.append(f"file_b contains {vcf_b.duplicate_observations} duplicate normalized allele observations; kept the first occurrence.")

        variants_a = set(vcf_a.observations)
        variants_b = set(vcf_b.observations)
        shared_variants = sorted(variants_a & variants_b, key=lambda key: (key.chrom, key.pos, key.ref, key.alt))
        variant_jaccard = jaccard(variants_a, variants_b)
        shared_samples = [sample for sample in vcf_a.samples if sample in set(vcf_b.samples)]
        sample_overlap = jaccard(vcf_a.samples, vcf_b.samples)
        if vcf_a.samples or vcf_b.samples:
            warnings.append("Genotype comparison ignores phasing and compares allele-specific alternate dosage.")
        if (vcf_a.samples or vcf_b.samples) and not shared_samples:
            warnings.append("No shared VCF samples were found; genotype metrics are unavailable.")

        genotype_matches, genotype_compared = genotype_concordance_counts(vcf_a, vcf_b, shared_variants, shared_samples)
        genotype_concordance = genotype_matches / genotype_compared if genotype_compared else variant_jaccard

        af_a: list[float] = []
        af_b: list[float] = []
        af_samples_a = shared_samples or vcf_a.samples
        af_samples_b = shared_samples or vcf_b.samples
        for key in shared_variants:
            freq_a = alternate_allele_frequency(vcf_a.observations[key], af_samples_a)
            freq_b = alternate_allele_frequency(vcf_b.observations[key], af_samples_b)
            if freq_a is None or freq_b is None:
                continue
            af_a.append(freq_a)
            af_b.append(freq_b)
        af_pearson = pearson(af_a, af_b) if af_a else 0.0
        af_spearman = spearman(af_a, af_b) if af_a else 0.0
        af_score = clamp01((af_spearman + 1.0) / 2.0) if len(af_a) >= 2 else genotype_concordance

        titv_a = titv_ratio(vcf_a.observations.values())
        titv_b = titv_ratio(vcf_b.observations.values())
        titv_similarity = magnitude_ratio(titv_a, titv_b)

        metrics = {
            "variant_jaccard": variant_jaccard,
            "sample_overlap": sample_overlap,
            "genotype_concordance": genotype_concordance,
            "allele_frequency_pearson": af_pearson,
            "allele_frequency_spearman": af_spearman,
            "titv_ratio_file_a": titv_a,
            "titv_ratio_file_b": titv_b,
            "titv_ratio_similarity": titv_similarity,
        }
        gt_score = genotype_concordance
        overall = clamp01(0.40 * variant_jaccard + 0.35 * gt_score + 0.15 * af_score + 0.10 * titv_similarity)
        details = {
            "normalization": "split_alt_alleles_and_trim_shared_bases",
            "file_a_raw_records": len(vcf_a.records),
            "file_b_raw_records": len(vcf_b.records),
            "file_a_variants": len(vcf_a.observations),
            "file_b_variants": len(vcf_b.observations),
            "shared_variants": len(shared_variants),
            "file_a_only_variants": len(variants_a - variants_b),
            "file_b_only_variants": len(variants_b - variants_a),
            "file_a_samples": vcf_a.samples,
            "file_b_samples": vcf_b.samples,
            "shared_samples": shared_samples,
            "genotype_calls_compared": genotype_compared,
            "allele_frequency_variants_compared": len(af_a),
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


def read_vcf(path: str) -> VCFData:
    samples: list[str] = []
    records: dict[VariantKey, VariantRecord] = {}
    duplicates = 0
    for line_number, line in enumerate(read_text_lines(path), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("##"):
            continue
        if stripped.startswith("#CHROM"):
            fields = stripped.split("\t")
            samples = fields[9:] if len(fields) > 9 else []
            continue
        fields = stripped.split("\t")
        if len(fields) < 8:
            raise ValueError(f"{path!r} line {line_number} has fewer than 8 VCF columns")
        chrom, pos_raw, identifier, ref, alt, qual, filt, info = fields[:8]
        try:
            pos = int(pos_raw)
        except ValueError as exc:
            raise ValueError(f"{path!r} line {line_number} has non-integer POS") from exc
        format_keys = fields[8].split(":") if len(fields) > 8 else []
        sample_values = fields[9:]
        sample_map: dict[str, dict[str, str]] = {}
        for sample, raw_value in zip(samples, sample_values):
            values = raw_value.split(":")
            sample_map[sample] = {key: values[index] if index < len(values) else "" for index, key in enumerate(format_keys)}
        key = VariantKey(chrom=chrom, pos=pos, ref=ref, alt=alt)
        if key in records:
            duplicates += 1
            continue
        records[key] = VariantRecord(
            key=key,
            identifier=identifier,
            qual=qual,
            filt=filt,
            info=info,
            format_keys=format_keys,
            samples=sample_map,
        )
    observations, duplicate_observations = build_observations(records)
    return VCFData(
        path=str(path),
        samples=samples,
        records=records,
        observations=observations,
        duplicate_variants=duplicates,
        duplicate_observations=duplicate_observations,
    )


def read_text_lines(path: str) -> list[str]:
    file_path = Path(path)
    opener = gzip.open if file_path.name.lower().endswith(".gz") else open
    with opener(file_path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.readlines()


def genotype_concordance_counts(
    vcf_a: VCFData,
    vcf_b: VCFData,
    shared_variants: list[VariantKey],
    shared_samples: list[str],
) -> tuple[int, int]:
    matches = 0
    compared = 0
    for key in shared_variants:
        observation_a = vcf_a.observations[key]
        observation_b = vcf_b.observations[key]
        for sample in shared_samples:
            dosage_a = alternate_dosage(observation_a.record.samples.get(sample, {}).get("GT"), observation_a.alt_index)
            dosage_b = alternate_dosage(observation_b.record.samples.get(sample, {}).get("GT"), observation_b.alt_index)
            if dosage_a is None or dosage_b is None:
                continue
            compared += 1
            if dosage_a == dosage_b:
                matches += 1
    return matches, compared


def build_observations(records: dict[VariantKey, VariantRecord]) -> tuple[dict[VariantKey, VariantObservation], int]:
    observations: dict[VariantKey, VariantObservation] = {}
    duplicates = 0
    for record in records.values():
        for alt_index, alt in enumerate(record.key.alt.split(","), start=1):
            if alt in {"", ".", "*"}:
                continue
            normalized_key = normalize_variant_key(record.key.chrom, record.key.pos, record.key.ref, alt)
            if normalized_key in observations:
                duplicates += 1
                continue
            observations[normalized_key] = VariantObservation(key=normalized_key, record=record, alt_index=alt_index)
    return observations, duplicates


def normalize_variant_key(chrom: str, pos: int, ref: str, alt: str) -> VariantKey:
    normalized_ref = ref.upper()
    normalized_alt = alt.upper()
    normalized_pos = pos
    if is_symbolic_alt(normalized_alt):
        return VariantKey(chrom=chrom, pos=normalized_pos, ref=normalized_ref, alt=normalized_alt)

    while len(normalized_ref) > 1 and len(normalized_alt) > 1 and normalized_ref[-1] == normalized_alt[-1]:
        normalized_ref = normalized_ref[:-1]
        normalized_alt = normalized_alt[:-1]
    while len(normalized_ref) > 1 and len(normalized_alt) > 1 and normalized_ref[0] == normalized_alt[0]:
        normalized_ref = normalized_ref[1:]
        normalized_alt = normalized_alt[1:]
        normalized_pos += 1
    return VariantKey(chrom=chrom, pos=normalized_pos, ref=normalized_ref, alt=normalized_alt)


def is_symbolic_alt(alt: str) -> bool:
    return alt.startswith("<") or "[" in alt or "]" in alt


def alternate_dosage(genotype: str | None, alt_index: int) -> int | None:
    if genotype is None or genotype in {"", ".", "./.", ".|."}:
        return None
    allele_tokens = GENOTYPE_SPLIT_RE.split(genotype)
    if not allele_tokens or any(token == "." for token in allele_tokens):
        return None
    return sum(token == str(alt_index) for token in allele_tokens)


def alternate_allele_frequency(observation: VariantObservation, samples: list[str]) -> float | None:
    alt_alleles = 0
    called_alleles = 0
    for sample in samples:
        genotype = observation.record.samples.get(sample, {}).get("GT")
        if genotype is None:
            continue
        allele_tokens = GENOTYPE_SPLIT_RE.split(genotype)
        for token in allele_tokens:
            if token == "." or token == "":
                continue
            called_alleles += 1
            if token == str(observation.alt_index):
                alt_alleles += 1
    if called_alleles == 0:
        return None
    return alt_alleles / called_alleles


def titv_ratio(observations) -> float:
    transitions = 0
    transversions = 0
    for observation in observations:
        if len(observation.key.ref) != 1 or len(observation.key.alt) != 1:
            continue
        ref = observation.key.ref.upper()
        alt_base = observation.key.alt.upper()
        if ref == alt_base:
            continue
        if (ref, alt_base) in TRANSITIONS:
            transitions += 1
        else:
            transversions += 1
    if transitions == 0 and transversions == 0:
        return 0.0
    if transversions == 0:
        return float(transitions)
    return transitions / transversions


def magnitude_ratio(left: float, right: float) -> float:
    if left == 0 and right == 0:
        return 1.0
    return min(abs(left), abs(right)) / max(abs(left), abs(right))

from __future__ import annotations

from dataclasses import dataclass

from biocompare.comparators.base import Comparator
from biocompare.comparators.table import ParsedTable, parse_float, read_table
from biocompare.core.report import ConcordanceReport
from biocompare.core.utils import clamp01, jaccard, mean_absolute_error, numeric_similarity, pearson, spearman
from biocompare.detection.filetype import detect_file_type


GENE_ALIASES = {
    "gene",
    "genes",
    "geneid",
    "gene_id",
    "feature",
    "featureid",
    "feature_id",
    "target_id",
    "transcript",
    "transcript_id",
    "id",
}


@dataclass(slots=True)
class CountMatrix:
    table: ParsedTable
    gene_column: str
    sample_columns: list[str]
    rows_by_gene: dict[str, dict[str, str]]
    duplicate_genes: int


class CountsComparator(Comparator):
    """Comparator for gene-by-sample count or expression matrices."""

    name = "counts"
    supported_types = ("counts", "expression")

    def can_handle(self, file_a: str, file_b: str, **kwargs: object) -> bool:
        requested_type = kwargs.get("file_type")
        if requested_type in self.supported_types:
            return True
        if requested_type is not None:
            return False

        type_a = detect_file_type(file_a)
        type_b = detect_file_type(file_b)
        if type_a.kind not in {"table", "csv", "tsv"} or type_b.kind not in {"table", "csv", "tsv"}:
            return False
        try:
            matrix_a = parse_count_matrix(read_table(file_a))
            matrix_b = parse_count_matrix(read_table(file_b))
        except ValueError:
            return False
        return len(set(matrix_a.sample_columns) & set(matrix_b.sample_columns)) >= 2

    def compare(self, file_a: str, file_b: str, **kwargs: object) -> ConcordanceReport:
        delimiter = kwargs.get("delimiter")
        if delimiter is not None and not isinstance(delimiter, str):
            raise TypeError("delimiter must be a string")
        gene_column = string_kwarg(kwargs, "gene_column")
        sample_columns = sample_columns_kwarg(kwargs.get("sample_columns"))

        matrix_a = parse_count_matrix(read_table(file_a, delimiter=delimiter), gene_column=gene_column, sample_columns=sample_columns)
        matrix_b = parse_count_matrix(read_table(file_b, delimiter=delimiter), gene_column=gene_column, sample_columns=sample_columns)

        warnings: list[str] = []
        if matrix_a.duplicate_genes:
            warnings.append(f"file_a contains {matrix_a.duplicate_genes} duplicate gene identifiers; kept the first occurrence.")
        if matrix_b.duplicate_genes:
            warnings.append(f"file_b contains {matrix_b.duplicate_genes} duplicate gene identifiers; kept the first occurrence.")

        genes_a = set(matrix_a.rows_by_gene)
        genes_b = set(matrix_b.rows_by_gene)
        shared_genes = sorted(genes_a & genes_b)
        samples_a = set(matrix_a.sample_columns)
        samples_b = set(matrix_b.sample_columns)
        shared_samples = [sample for sample in matrix_a.sample_columns if sample in samples_b]
        if not shared_genes:
            warnings.append("No shared gene identifiers were found.")
        if not shared_samples:
            warnings.append("No shared sample columns were found.")

        gene_overlap = jaccard(genes_a, genes_b)
        sample_overlap = jaccard(samples_a, samples_b)

        sample_pearsons: list[float] = []
        sample_spearmans: list[float] = []
        sample_similarities: list[float] = []
        sample_maes: list[float] = []
        library_size_ratios: list[float] = []
        metrics: dict[str, float] = {
            "gene_overlap": gene_overlap,
            "sample_overlap": sample_overlap,
        }

        for sample in shared_samples:
            values_a, values_b = sample_vectors(matrix_a, matrix_b, shared_genes, sample)
            if not values_a:
                continue
            sample_pearson = pearson(values_a, values_b)
            sample_spearman = spearman(values_a, values_b)
            sample_similarity = numeric_similarity(values_a, values_b)
            sample_mae = mean_absolute_error(values_a, values_b)
            library_ratio = magnitude_ratio(sum(values_a), sum(values_b))
            metrics[f"sample.{sample}.pearson"] = sample_pearson
            metrics[f"sample.{sample}.spearman"] = sample_spearman
            metrics[f"sample.{sample}.similarity"] = sample_similarity
            metrics[f"sample.{sample}.mae"] = sample_mae
            metrics[f"sample.{sample}.library_size_ratio"] = library_ratio
            sample_pearsons.append(sample_pearson)
            sample_spearmans.append(sample_spearman)
            sample_similarities.append(sample_similarity)
            sample_maes.append(sample_mae)
            library_size_ratios.append(library_ratio)

        gene_profile_pearsons: list[float] = []
        gene_profile_spearmans: list[float] = []
        if len(shared_samples) >= 2:
            for gene in shared_genes:
                values_a, values_b = gene_vectors(matrix_a, matrix_b, gene, shared_samples)
                if len(values_a) >= 2:
                    gene_profile_pearsons.append(pearson(values_a, values_b))
                    gene_profile_spearmans.append(spearman(values_a, values_b))

        zero_jaccard = zero_pattern_jaccard(matrix_a, matrix_b, shared_genes, shared_samples)
        metrics["mean_sample_pearson"] = mean(sample_pearsons)
        metrics["mean_sample_spearman"] = mean(sample_spearmans)
        metrics["mean_sample_similarity"] = mean(sample_similarities)
        metrics["mean_sample_mae"] = mean(sample_maes)
        metrics["mean_library_size_ratio"] = mean(library_size_ratios)
        metrics["mean_gene_profile_pearson"] = mean(gene_profile_pearsons)
        metrics["mean_gene_profile_spearman"] = mean(gene_profile_spearmans)
        metrics["zero_pattern_jaccard"] = zero_jaccard

        sample_corr_score = clamp01((metrics["mean_sample_spearman"] + 1.0) / 2.0)
        gene_profile_score = clamp01((metrics["mean_gene_profile_spearman"] + 1.0) / 2.0)
        overall = clamp01(
            0.20 * gene_overlap
            + 0.15 * sample_overlap
            + 0.25 * sample_corr_score
            + 0.20 * metrics["mean_sample_similarity"]
            + 0.10 * metrics["mean_library_size_ratio"]
            + 0.10 * gene_profile_score
        )

        details = {
            "file_a_rows": len(matrix_a.table.rows),
            "file_b_rows": len(matrix_b.table.rows),
            "file_a_gene_column": matrix_a.gene_column,
            "file_b_gene_column": matrix_b.gene_column,
            "file_a_sample_columns": matrix_a.sample_columns,
            "file_b_sample_columns": matrix_b.sample_columns,
            "shared_genes": len(shared_genes),
            "file_a_only_genes": len(genes_a - genes_b),
            "file_b_only_genes": len(genes_b - genes_a),
            "shared_samples": shared_samples,
            "file_a_only_samples": sorted(samples_a - samples_b),
            "file_b_only_samples": sorted(samples_b - samples_a),
            "compared_cells": len(shared_genes) * len(shared_samples),
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


def parse_count_matrix(
    table: ParsedTable,
    *,
    gene_column: str | None = None,
    sample_columns: list[str] | None = None,
) -> CountMatrix:
    selected_gene_column = gene_column or detect_gene_column(table)
    if selected_gene_column not in table.columns:
        raise ValueError(f"gene column {selected_gene_column!r} is not present in {table.path!r}")
    selected_sample_columns = sample_columns or detect_sample_columns(table, selected_gene_column)
    missing_samples = [sample for sample in selected_sample_columns if sample not in table.columns]
    if missing_samples:
        raise ValueError(f"sample columns are not present in {table.path!r}: {', '.join(missing_samples)}")
    if len(selected_sample_columns) < 2:
        raise ValueError(f"Count matrix {table.path!r} needs at least two numeric sample columns")

    rows, duplicates = rows_by_gene(table, selected_gene_column)
    validate_numeric_matrix(rows, selected_sample_columns, table.path)
    return CountMatrix(
        table=table,
        gene_column=selected_gene_column,
        sample_columns=selected_sample_columns,
        rows_by_gene=rows,
        duplicate_genes=duplicates,
    )


def detect_gene_column(table: ParsedTable) -> str:
    normalized_lookup = {normalize_column(column): column for column in table.columns}
    for alias in GENE_ALIASES:
        normalized_alias = normalize_column(alias)
        if normalized_alias in normalized_lookup:
            return normalized_lookup[normalized_alias]

    first_column = table.columns[0]
    first_values = [row.get(first_column, "") for row in table.rows if row.get(first_column, "") != ""]
    if first_values and len(first_values) == len(set(first_values)):
        numeric_values = [parse_float(value) for value in first_values]
        if not all(value is not None for value in numeric_values):
            return first_column
    raise ValueError(f"Could not detect gene identifier column in {table.path!r}")


def detect_sample_columns(table: ParsedTable, gene_column: str) -> list[str]:
    sample_columns: list[str] = []
    for column in table.columns:
        if column == gene_column:
            continue
        values = [row.get(column, "") for row in table.rows]
        if values and all(parse_float(value) is not None for value in values if value != "") and all(value != "" for value in values):
            sample_columns.append(column)
    return sample_columns


def rows_by_gene(table: ParsedTable, gene_column: str) -> tuple[dict[str, dict[str, str]], int]:
    rows: dict[str, dict[str, str]] = {}
    duplicates = 0
    for row in table.rows:
        gene = row.get(gene_column, "").strip()
        if not gene:
            continue
        if gene in rows:
            duplicates += 1
            continue
        rows[gene] = row
    return rows, duplicates


def validate_numeric_matrix(rows: dict[str, dict[str, str]], sample_columns: list[str], path: str) -> None:
    for gene, row in rows.items():
        for sample in sample_columns:
            if parse_float(row.get(sample, "")) is None:
                raise ValueError(f"Non-numeric value for gene {gene!r}, sample {sample!r} in {path!r}")


def sample_vectors(
    matrix_a: CountMatrix,
    matrix_b: CountMatrix,
    genes: list[str],
    sample: str,
) -> tuple[list[float], list[float]]:
    values_a: list[float] = []
    values_b: list[float] = []
    for gene in genes:
        value_a = parse_float(matrix_a.rows_by_gene[gene][sample])
        value_b = parse_float(matrix_b.rows_by_gene[gene][sample])
        if value_a is None or value_b is None:
            continue
        values_a.append(value_a)
        values_b.append(value_b)
    return values_a, values_b


def gene_vectors(
    matrix_a: CountMatrix,
    matrix_b: CountMatrix,
    gene: str,
    samples: list[str],
) -> tuple[list[float], list[float]]:
    values_a: list[float] = []
    values_b: list[float] = []
    for sample in samples:
        value_a = parse_float(matrix_a.rows_by_gene[gene][sample])
        value_b = parse_float(matrix_b.rows_by_gene[gene][sample])
        if value_a is None or value_b is None:
            continue
        values_a.append(value_a)
        values_b.append(value_b)
    return values_a, values_b


def zero_pattern_jaccard(
    matrix_a: CountMatrix,
    matrix_b: CountMatrix,
    genes: list[str],
    samples: list[str],
) -> float:
    zero_a: set[tuple[str, str]] = set()
    zero_b: set[tuple[str, str]] = set()
    for gene in genes:
        for sample in samples:
            value_a = parse_float(matrix_a.rows_by_gene[gene][sample])
            value_b = parse_float(matrix_b.rows_by_gene[gene][sample])
            if value_a == 0:
                zero_a.add((gene, sample))
            if value_b == 0:
                zero_b.add((gene, sample))
    return jaccard(zero_a, zero_b)


def magnitude_ratio(left: float, right: float) -> float:
    if left == 0 and right == 0:
        return 1.0
    return min(abs(left), abs(right)) / max(abs(left), abs(right))


def sample_columns_kwarg(value: object) -> list[str] | None:
    if value is None:
        return None
    if isinstance(value, str):
        columns = [column.strip() for column in value.split(",") if column.strip()]
    elif isinstance(value, list):
        columns = value
    else:
        raise TypeError("sample_columns must be a comma-separated string or list of strings")
    if not columns or not all(isinstance(column, str) for column in columns):
        raise ValueError("sample_columns must contain at least one string column name")
    return columns


def string_kwarg(kwargs: dict[str, object], name: str) -> str | None:
    value = kwargs.get(name)
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a string")
    return value


def normalize_column(column: str) -> str:
    return "".join(character for character in column.lower() if character.isalnum())


def mean(values: list[float]) -> float:
    if not values:
        return 0.0
    return sum(values) / len(values)

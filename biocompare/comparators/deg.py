from __future__ import annotations

from dataclasses import dataclass

from biocompare.comparators.base import Comparator
from biocompare.comparators.table import ParsedTable, parse_float, read_table
from biocompare.core.report import ConcordanceReport
from biocompare.core.utils import clamp01, jaccard, mean_absolute_error, pearson, spearman
from biocompare.detection.filetype import detect_file_type


GENE_ALIASES = {
    "gene",
    "genes",
    "geneid",
    "gene_id",
    "geneidentifier",
    "gene_name",
    "genename",
    "symbol",
    "genesymbol",
    "feature",
    "featureid",
    "feature_id",
    "id",
}
LOGFC_ALIASES = {
    "logfc",
    "log2fc",
    "log2foldchange",
    "logfoldchange",
    "log_fold_change",
    "lfc",
    "estimate",
}
PADJ_ALIASES = {
    "padj",
    "adjp",
    "adjpvalue",
    "adjustedpvalue",
    "adjusted_p_value",
    "p_adjust",
    "p_adj",
    "fdr",
    "qvalue",
    "qval",
}
PVALUE_ALIASES = {
    "pvalue",
    "pval",
    "p_value",
    "p.value",
    "p",
}


@dataclass(slots=True)
class DegColumns:
    gene: str
    logfc: str
    pvalue: str
    uses_adjusted_pvalue: bool


class DEGComparator(Comparator):
    """Comparator for differential expression result tables."""

    name = "deg"
    supported_types = ("deg",)

    def can_handle(self, file_a: str, file_b: str, **kwargs: object) -> bool:
        requested_type = kwargs.get("file_type")
        if requested_type == "deg":
            return True
        if requested_type is not None:
            return False

        type_a = detect_file_type(file_a)
        type_b = detect_file_type(file_b)
        if type_a.kind not in {"table", "csv", "tsv"} or type_b.kind not in {"table", "csv", "tsv"}:
            return False
        try:
            table_a = read_table(file_a)
            table_b = read_table(file_b)
            detect_deg_columns(table_a)
            detect_deg_columns(table_b)
        except ValueError:
            return False
        return True

    def compare(self, file_a: str, file_b: str, **kwargs: object) -> ConcordanceReport:
        alpha = optional_float(kwargs.get("alpha"), 0.05, "alpha")
        lfc_threshold = optional_float(kwargs.get("lfc_threshold"), 0.0, "lfc_threshold")
        top_n = optional_int(kwargs.get("top_n"), 50, "top_n")
        delimiter = kwargs.get("delimiter")
        if delimiter is not None and not isinstance(delimiter, str):
            raise TypeError("delimiter must be a string")

        table_a = read_table(file_a, delimiter=delimiter)
        table_b = read_table(file_b, delimiter=delimiter)
        columns_a = detect_deg_columns(
            table_a,
            gene_column=string_kwarg(kwargs, "gene_column"),
            logfc_column=string_kwarg(kwargs, "logfc_column"),
            pvalue_column=string_kwarg(kwargs, "pvalue_column"),
            padj_column=string_kwarg(kwargs, "padj_column"),
        )
        columns_b = detect_deg_columns(
            table_b,
            gene_column=string_kwarg(kwargs, "gene_column"),
            logfc_column=string_kwarg(kwargs, "logfc_column"),
            pvalue_column=string_kwarg(kwargs, "pvalue_column"),
            padj_column=string_kwarg(kwargs, "padj_column"),
        )

        warnings: list[str] = []
        if not columns_a.uses_adjusted_pvalue:
            warnings.append("file_a has no adjusted p-value column; using raw p-value for significance.")
        if not columns_b.uses_adjusted_pvalue:
            warnings.append("file_b has no adjusted p-value column; using raw p-value for significance.")

        rows_a, duplicate_a = rows_by_gene(table_a, columns_a.gene)
        rows_b, duplicate_b = rows_by_gene(table_b, columns_b.gene)
        if duplicate_a:
            warnings.append(f"file_a contains {duplicate_a} duplicate gene identifiers; kept the first occurrence.")
        if duplicate_b:
            warnings.append(f"file_b contains {duplicate_b} duplicate gene identifiers; kept the first occurrence.")

        genes_a = set(rows_a)
        genes_b = set(rows_b)
        shared_genes = sorted(genes_a & genes_b)
        row_overlap = jaccard(genes_a, genes_b)
        if not shared_genes:
            warnings.append("No shared gene identifiers were found.")

        shared_lfc_a: list[float] = []
        shared_lfc_b: list[float] = []
        direction_pairs: list[tuple[int, int]] = []
        for gene in shared_genes:
            lfc_a = parse_float(rows_a[gene].get(columns_a.logfc, ""))
            lfc_b = parse_float(rows_b[gene].get(columns_b.logfc, ""))
            if lfc_a is None or lfc_b is None:
                continue
            shared_lfc_a.append(lfc_a)
            shared_lfc_b.append(lfc_b)
            sign_a = direction(lfc_a, lfc_threshold)
            sign_b = direction(lfc_b, lfc_threshold)
            if sign_a != 0 or sign_b != 0:
                direction_pairs.append((sign_a, sign_b))

        sig_a = significant_genes(rows_a, columns_a, alpha=alpha, lfc_threshold=lfc_threshold)
        sig_b = significant_genes(rows_b, columns_b, alpha=alpha, lfc_threshold=lfc_threshold)
        shared_sig = sorted(sig_a & sig_b)
        up_a = directional_gene_set(rows_a, columns_a, sig_a, 1, lfc_threshold)
        up_b = directional_gene_set(rows_b, columns_b, sig_b, 1, lfc_threshold)
        down_a = directional_gene_set(rows_a, columns_a, sig_a, -1, lfc_threshold)
        down_b = directional_gene_set(rows_b, columns_b, sig_b, -1, lfc_threshold)

        top_a = top_genes(rows_a, columns_a, top_n)
        top_b = top_genes(rows_b, columns_b, top_n)
        logfc_pearson = pearson(shared_lfc_a, shared_lfc_b) if shared_lfc_a else 0.0
        logfc_spearman = spearman(shared_lfc_a, shared_lfc_b) if shared_lfc_a else 0.0
        logfc_mae = mean_absolute_error(shared_lfc_a, shared_lfc_b) if shared_lfc_a else 0.0
        shared_direction_agreement = direction_agreement(direction_pairs)
        significant_direction_agreement = significant_direction(rows_a, rows_b, columns_a, columns_b, shared_sig, lfc_threshold, sig_a, sig_b)
        significant_jaccard = jaccard(sig_a, sig_b)
        lfc_score = clamp01((logfc_spearman + 1.0) / 2.0)
        overall = clamp01(
            0.15 * row_overlap
            + 0.35 * significant_jaccard
            + 0.25 * significant_direction_agreement
            + 0.25 * lfc_score
        )

        metrics = {
            "gene_overlap": row_overlap,
            "significant_jaccard": significant_jaccard,
            "significant_precision_file_b_vs_a": safe_fraction(len(sig_a & sig_b), len(sig_b), default=1.0 if not sig_a else 0.0),
            "significant_recall_file_b_vs_a": safe_fraction(len(sig_a & sig_b), len(sig_a), default=1.0 if not sig_b else 0.0),
            "significant_f1_file_b_vs_a": f1_score(len(sig_a & sig_b), len(sig_b), len(sig_a)),
            "upregulated_jaccard": jaccard(up_a, up_b),
            "downregulated_jaccard": jaccard(down_a, down_b),
            "top_genes_jaccard": jaccard(top_a, top_b),
            "logfc_pearson": logfc_pearson,
            "logfc_spearman": logfc_spearman,
            "logfc_mae": logfc_mae,
            "shared_direction_agreement": shared_direction_agreement,
            "significant_direction_agreement": significant_direction_agreement,
        }
        details = {
            "alpha": alpha,
            "lfc_threshold": lfc_threshold,
            "top_n": top_n,
            "file_a_rows": len(table_a.rows),
            "file_b_rows": len(table_b.rows),
            "file_a_columns": {
                "gene": columns_a.gene,
                "logfc": columns_a.logfc,
                "pvalue": columns_a.pvalue,
                "uses_adjusted_pvalue": columns_a.uses_adjusted_pvalue,
            },
            "file_b_columns": {
                "gene": columns_b.gene,
                "logfc": columns_b.logfc,
                "pvalue": columns_b.pvalue,
                "uses_adjusted_pvalue": columns_b.uses_adjusted_pvalue,
            },
            "shared_genes": len(shared_genes),
            "file_a_only_genes": len(genes_a - genes_b),
            "file_b_only_genes": len(genes_b - genes_a),
            "file_a_significant": len(sig_a),
            "file_b_significant": len(sig_b),
            "shared_significant": len(shared_sig),
            "file_a_upregulated": len(up_a),
            "file_b_upregulated": len(up_b),
            "file_a_downregulated": len(down_a),
            "file_b_downregulated": len(down_b),
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


def detect_deg_columns(
    table: ParsedTable,
    *,
    gene_column: str | None = None,
    logfc_column: str | None = None,
    pvalue_column: str | None = None,
    padj_column: str | None = None,
) -> DegColumns:
    gene = require_column(table, gene_column, GENE_ALIASES, "gene identifier")
    logfc = require_column(table, logfc_column, LOGFC_ALIASES, "log fold-change")
    padj = find_column(table.columns, padj_column, PADJ_ALIASES)
    if padj is not None:
        return DegColumns(gene=gene, logfc=logfc, pvalue=padj, uses_adjusted_pvalue=True)
    pvalue = require_column(table, pvalue_column, PVALUE_ALIASES, "p-value")
    return DegColumns(gene=gene, logfc=logfc, pvalue=pvalue, uses_adjusted_pvalue=False)


def require_column(table: ParsedTable, requested: str | None, aliases: set[str], label: str) -> str:
    column = find_column(table.columns, requested, aliases)
    if column is None:
        raise ValueError(f"Could not detect {label} column in {table.path!r}")
    return column


def find_column(columns: list[str], requested: str | None, aliases: set[str]) -> str | None:
    if requested is not None:
        if requested not in columns:
            raise ValueError(f"Requested column {requested!r} is not present")
        return requested
    normalized_lookup = {normalize_column(column): column for column in columns}
    for alias in aliases:
        normalized_alias = normalize_column(alias)
        if normalized_alias in normalized_lookup:
            return normalized_lookup[normalized_alias]
    return None


def normalize_column(column: str) -> str:
    return "".join(character for character in column.lower() if character.isalnum())


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


def significant_genes(
    rows: dict[str, dict[str, str]],
    columns: DegColumns,
    *,
    alpha: float,
    lfc_threshold: float,
) -> set[str]:
    significant: set[str] = set()
    for gene, row in rows.items():
        pvalue = parse_float(row.get(columns.pvalue, ""))
        logfc = parse_float(row.get(columns.logfc, ""))
        if pvalue is None or logfc is None:
            continue
        if pvalue <= alpha and abs(logfc) >= lfc_threshold:
            significant.add(gene)
    return significant


def directional_gene_set(
    rows: dict[str, dict[str, str]],
    columns: DegColumns,
    genes: set[str],
    target_direction: int,
    lfc_threshold: float,
) -> set[str]:
    selected: set[str] = set()
    for gene in genes:
        logfc = parse_float(rows[gene].get(columns.logfc, ""))
        if logfc is None:
            continue
        if direction(logfc, lfc_threshold) == target_direction:
            selected.add(gene)
    return selected


def top_genes(rows: dict[str, dict[str, str]], columns: DegColumns, top_n: int) -> set[str]:
    ranked: list[tuple[float, str]] = []
    for gene, row in rows.items():
        pvalue = parse_float(row.get(columns.pvalue, ""))
        if pvalue is not None:
            ranked.append((pvalue, gene))
    ranked.sort(key=lambda item: (item[0], item[1]))
    return {gene for _, gene in ranked[:top_n]}


def direction(value: float, threshold: float) -> int:
    if value > threshold:
        return 1
    if value < -threshold:
        return -1
    return 0


def direction_agreement(pairs: list[tuple[int, int]]) -> float:
    if not pairs:
        return 0.0
    return sum(left == right for left, right in pairs) / len(pairs)


def significant_direction(
    rows_a: dict[str, dict[str, str]],
    rows_b: dict[str, dict[str, str]],
    columns_a: DegColumns,
    columns_b: DegColumns,
    shared_sig: list[str],
    lfc_threshold: float,
    sig_a: set[str],
    sig_b: set[str],
) -> float:
    if not shared_sig:
        return 1.0 if not sig_a and not sig_b else 0.0
    pairs: list[tuple[int, int]] = []
    for gene in shared_sig:
        lfc_a = parse_float(rows_a[gene].get(columns_a.logfc, ""))
        lfc_b = parse_float(rows_b[gene].get(columns_b.logfc, ""))
        if lfc_a is None or lfc_b is None:
            continue
        pairs.append((direction(lfc_a, lfc_threshold), direction(lfc_b, lfc_threshold)))
    return direction_agreement(pairs)


def safe_fraction(numerator: int, denominator: int, *, default: float) -> float:
    if denominator == 0:
        return default
    return numerator / denominator


def f1_score(true_positive: int, predicted: int, expected: int) -> float:
    precision = safe_fraction(true_positive, predicted, default=1.0 if expected == 0 else 0.0)
    recall = safe_fraction(true_positive, expected, default=1.0 if predicted == 0 else 0.0)
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


def optional_int(value: object, default: int, name: str) -> int:
    if value is None:
        return default
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if parsed <= 0:
        raise ValueError(f"{name} must be positive")
    return parsed


def string_kwarg(kwargs: dict[str, object], name: str) -> str | None:
    value = kwargs.get(name)
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a string")
    return value


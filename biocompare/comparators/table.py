from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from pathlib import Path

from biocompare.comparators.base import Comparator
from biocompare.core.report import ConcordanceReport
from biocompare.core.utils import (
    agreement_rate,
    clamp01,
    cohen_kappa,
    jaccard,
    mean_absolute_error,
    numeric_similarity,
    pearson,
    spearman,
)
from biocompare.detection.filetype import detect_file_type, sniff_delimiter


@dataclass(slots=True)
class ParsedTable:
    path: str
    columns: list[str]
    rows: list[dict[str, str]]
    delimiter: str


class TableComparator(Comparator):
    """Generic comparator for delimited tabular outputs."""

    name = "table"
    supported_types = ("table", "csv", "tsv")

    def can_handle(self, file_a: str, file_b: str, **kwargs: object) -> bool:
        requested_type = kwargs.get("file_type")
        if requested_type in self.supported_types:
            return True
        type_a = detect_file_type(file_a)
        type_b = detect_file_type(file_b)
        table_types = set(self.supported_types)
        return type_a.kind in table_types and type_b.kind in table_types

    def compare(self, file_a: str, file_b: str, **kwargs: object) -> ConcordanceReport:
        key_column = kwargs.get("key_column")
        if key_column is not None and not isinstance(key_column, str):
            raise TypeError("key_column must be a string")
        delimiter = kwargs.get("delimiter")
        if delimiter is not None and not isinstance(delimiter, str):
            raise TypeError("delimiter must be a string")

        table_a = read_table(file_a, delimiter=delimiter)
        table_b = read_table(file_b, delimiter=delimiter)
        warnings: list[str] = []

        shared_columns = [column for column in table_a.columns if column in table_b.columns]
        if not shared_columns:
            raise ValueError("tables do not share any column names")

        selected_key = key_column or auto_key_column(table_a, table_b)
        if selected_key and selected_key not in shared_columns:
            raise ValueError(f"key column {selected_key!r} is not present in both tables")

        pairs, row_overlap, row_details = align_rows(table_a, table_b, selected_key)
        if selected_key is None and len(table_a.rows) != len(table_b.rows):
            warnings.append("No key column was detected; compared rows by position and ignored extra rows.")
        if selected_key is not None and row_details["compared_rows"] == 0:
            warnings.append(f"No shared row identifiers found in key column {selected_key!r}.")

        compared_columns = [column for column in shared_columns if column != selected_key]
        only_a = [column for column in table_a.columns if column not in table_b.columns]
        only_b = [column for column in table_b.columns if column not in table_a.columns]
        if only_a:
            warnings.append(f"Columns only in file_a were not compared: {', '.join(only_a)}")
        if only_b:
            warnings.append(f"Columns only in file_b were not compared: {', '.join(only_b)}")

        metrics: dict[str, float] = {"row_overlap": row_overlap}
        numeric_columns: list[str] = []
        categorical_columns: list[str] = []
        column_scores: list[float] = []

        for column in compared_columns:
            left_values, right_values = paired_values(pairs, column)
            if not left_values:
                warnings.append(f"Column {column!r} had no paired non-empty values.")
                continue

            numeric_pairs = [(parse_float(left), parse_float(right)) for left, right in zip(left_values, right_values)]
            if all(left is not None and right is not None for left, right in numeric_pairs):
                left_numbers = [left for left, _ in numeric_pairs if left is not None]
                right_numbers = [right for _, right in numeric_pairs if right is not None]
                score = numeric_similarity(left_numbers, right_numbers)
                metrics[f"numeric.{column}.pearson"] = pearson(left_numbers, right_numbers)
                metrics[f"numeric.{column}.spearman"] = spearman(left_numbers, right_numbers)
                metrics[f"numeric.{column}.mae"] = mean_absolute_error(left_numbers, right_numbers)
                metrics[f"numeric.{column}.similarity"] = score
                numeric_columns.append(column)
                column_scores.append(score)
            else:
                agreement = agreement_rate(left_values, right_values)
                kappa = cohen_kappa(left_values, right_values)
                metrics[f"categorical.{column}.agreement"] = agreement
                metrics[f"categorical.{column}.cohen_kappa"] = kappa
                categorical_columns.append(column)
                column_scores.append(clamp01((agreement + max(kappa, 0.0)) / 2.0))

        column_score_mean = sum(column_scores) / len(column_scores) if column_scores else 0.0
        metrics["column_score_mean"] = column_score_mean
        overall = row_overlap if not column_scores else clamp01(0.25 * row_overlap + 0.75 * column_score_mean)

        details = {
            "file_a_rows": len(table_a.rows),
            "file_b_rows": len(table_b.rows),
            "file_a_columns": table_a.columns,
            "file_b_columns": table_b.columns,
            "shared_columns": shared_columns,
            "compared_columns": compared_columns,
            "numeric_columns": numeric_columns,
            "categorical_columns": categorical_columns,
            "key_column": selected_key,
            "file_a_delimiter": printable_delimiter(table_a.delimiter),
            "file_b_delimiter": printable_delimiter(table_b.delimiter),
            **row_details,
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


def read_table(path: str, delimiter: str | None = None) -> ParsedTable:
    table_path = Path(path)
    text = table_path.read_text(encoding="utf-8-sig")
    if not text.strip():
        raise ValueError(f"{path!r} is empty")
    detected_delimiter = delimiter or sniff_delimiter(text, filename=str(path))
    if detected_delimiter is None:
        raise ValueError(f"Could not detect a CSV/TSV delimiter for {path!r}")

    reader = csv.reader(io.StringIO(text), delimiter=detected_delimiter)
    try:
        header = next(reader)
    except StopIteration as exc:
        raise ValueError(f"{path!r} does not contain a header row") from exc

    columns = [column.strip() for column in header]
    if not all(columns):
        raise ValueError(f"{path!r} contains an empty column name")
    if len(set(columns)) != len(columns):
        raise ValueError(f"{path!r} contains duplicate column names")

    rows: list[dict[str, str]] = []
    for raw_row in reader:
        if not raw_row or not any(value.strip() for value in raw_row):
            continue
        padded_row = list(raw_row[: len(columns)])
        if len(padded_row) < len(columns):
            padded_row.extend([""] * (len(columns) - len(padded_row)))
        rows.append({column: padded_row[index].strip() for index, column in enumerate(columns)})

    return ParsedTable(path=str(path), columns=columns, rows=rows, delimiter=detected_delimiter)


def auto_key_column(table_a: ParsedTable, table_b: ParsedTable) -> str | None:
    shared_columns = [column for column in table_a.columns if column in table_b.columns]
    preferred = {
        "id",
        "gene",
        "gene_id",
        "geneid",
        "feature",
        "feature_id",
        "sample",
        "sample_id",
        "variant",
        "name",
    }
    ordered_columns = sorted(shared_columns, key=lambda column: (column.lower() not in preferred, table_a.columns.index(column)))
    for column in ordered_columns:
        values_a = [row[column] for row in table_a.rows if row[column] != ""]
        values_b = [row[column] for row in table_b.rows if row[column] != ""]
        if not values_a or not values_b:
            continue
        if len(values_a) == len(set(values_a)) and len(values_b) == len(set(values_b)):
            if set(values_a) & set(values_b):
                return column
    return None


def align_rows(
    table_a: ParsedTable,
    table_b: ParsedTable,
    key_column: str | None,
) -> tuple[list[tuple[dict[str, str], dict[str, str]]], float, dict[str, object]]:
    if key_column is None:
        compared_rows = min(len(table_a.rows), len(table_b.rows))
        max_rows = max(len(table_a.rows), len(table_b.rows))
        row_overlap = compared_rows / max_rows if max_rows else 1.0
        return (
            list(zip(table_a.rows[:compared_rows], table_b.rows[:compared_rows])),
            row_overlap,
            {
                "alignment": "position",
                "compared_rows": compared_rows,
                "file_a_only_rows": max(len(table_a.rows) - compared_rows, 0),
                "file_b_only_rows": max(len(table_b.rows) - compared_rows, 0),
            },
        )

    rows_a = {row[key_column]: row for row in table_a.rows if row[key_column] != ""}
    rows_b = {row[key_column]: row for row in table_b.rows if row[key_column] != ""}
    keys_a = set(rows_a)
    keys_b = set(rows_b)
    common_keys = sorted(keys_a & keys_b)
    row_overlap = jaccard(keys_a, keys_b)
    return (
        [(rows_a[key], rows_b[key]) for key in common_keys],
        row_overlap,
        {
            "alignment": "key",
            "compared_rows": len(common_keys),
            "file_a_only_rows": len(keys_a - keys_b),
            "file_b_only_rows": len(keys_b - keys_a),
        },
    )


def paired_values(
    pairs: list[tuple[dict[str, str], dict[str, str]]],
    column: str,
) -> tuple[list[str], list[str]]:
    left_values: list[str] = []
    right_values: list[str] = []
    for left_row, right_row in pairs:
        left_value = left_row.get(column, "")
        right_value = right_row.get(column, "")
        if left_value == "" or right_value == "":
            continue
        left_values.append(left_value)
        right_values.append(right_value)
    return left_values, right_values


def parse_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except ValueError:
        return None
    if parsed != parsed or parsed in (float("inf"), float("-inf")):
        return None
    return parsed


def printable_delimiter(delimiter: str) -> str:
    if delimiter == "\t":
        return "\\t"
    return delimiter


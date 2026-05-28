from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from biocompare.core.engine import ComparisonEngine
from biocompare.core.report import ConcordanceReport


@dataclass(slots=True)
class BatchItem:
    label: str
    file_a: str
    file_b: str
    file_type: str | None = None


@dataclass(slots=True)
class BatchResult:
    item: BatchItem
    report: ConcordanceReport | None = None
    error: str | None = None

    @property
    def status(self) -> str:
        return "error" if self.error else "ok"


def read_manifest(path: str) -> list[BatchItem]:
    manifest_path = Path(path)
    delimiter = "\t" if manifest_path.suffix.lower() in {".tsv", ".tab"} else ","
    with manifest_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"{path!r} is empty")
        required = {"file_a", "file_b"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path!r} is missing required columns: {', '.join(sorted(missing))}")
        items: list[BatchItem] = []
        for index, row in enumerate(reader, start=1):
            file_a = (row.get("file_a") or "").strip()
            file_b = (row.get("file_b") or "").strip()
            if not file_a or not file_b:
                raise ValueError(f"{path!r} manifest row {index} must include file_a and file_b")
            label = (row.get("label") or row.get("name") or f"comparison_{index}").strip()
            file_type = (row.get("type") or row.get("file_type") or "").strip() or None
            items.append(
                BatchItem(
                    label=label,
                    file_a=resolve_manifest_path(manifest_path, file_a),
                    file_b=resolve_manifest_path(manifest_path, file_b),
                    file_type=file_type,
                )
            )
    return items


def run_batch(
    manifest_path: str,
    *,
    engine: ComparisonEngine | None = None,
    stop_on_error: bool = False,
    default_file_type: str | None = None,
    **kwargs: object,
) -> list[BatchResult]:
    selected_engine = engine or ComparisonEngine()
    results: list[BatchResult] = []
    for item in read_manifest(manifest_path):
        try:
            report = selected_engine.compare(
                item.file_a,
                item.file_b,
                file_type=item.file_type or default_file_type,
                **kwargs,
            )
            results.append(BatchResult(item=item, report=report))
        except Exception as exc:
            if stop_on_error:
                raise
            results.append(BatchResult(item=item, error=str(exc)))
    return results


def batch_summary(results: list[BatchResult]) -> dict[str, Any]:
    successful = [result for result in results if result.report is not None]
    failed = [result for result in results if result.error is not None]
    concordances = [result.report.overall_concordance for result in successful if result.report is not None]
    return {
        "total": len(results),
        "ok": len(successful),
        "error": len(failed),
        "min_concordance": min(concordances) if concordances else None,
        "mean_concordance": sum(concordances) / len(concordances) if concordances else None,
        "max_concordance": max(concordances) if concordances else None,
    }


def resolve_manifest_path(manifest_path: Path, value: str) -> str:
    path = Path(value)
    if path.is_absolute():
        return str(path)
    candidate = manifest_path.parent / path
    if candidate.exists():
        return str(candidate)
    return value


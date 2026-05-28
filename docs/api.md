# API Reference

## `ConcordanceReport`

All comparators return `biocompare.core.report.ConcordanceReport`.

Fields:

- `comparator`: comparator class name.
- `file_a`, `file_b`: input paths.
- `overall_concordance`: 0-1 summary score.
- `metrics`: flat dictionary of numeric metrics.
- `details`: comparator-specific structured metadata.
- `warnings`: assumptions, skipped content, or caveats.

## `ComparisonEngine`

```python
from biocompare import ComparisonEngine

report = ComparisonEngine().compare(
    "tests/fixtures/degs_tool_a.tsv",
    "tests/fixtures/degs_tool_b.tsv",
    file_type="deg",
    alpha=0.05,
)
```

The engine registers built-in comparators and then asks the registry to find a
comparator that can handle the file pair. Pass `file_type` when auto-detection
would be ambiguous.

## Batch API

```python
from biocompare.core.batch import run_batch, batch_summary

results = run_batch("tests/fixtures/batch_manifest.tsv")
summary = batch_summary(results)
```

Manifest columns:

- `file_a`: required.
- `file_b`: required.
- `label`: optional.
- `type` or `file_type`: optional comparator override.

## Report Writers

```python
from biocompare.io.report_writers import report_to_json, report_to_html

json_text = report_to_json(report)
html_text = report_to_html(report)
```

Single reports support JSON, text, and HTML. Batch reports support JSON, TSV,
text, and HTML.


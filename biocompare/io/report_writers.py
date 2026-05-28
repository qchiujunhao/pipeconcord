from __future__ import annotations

import json
from pathlib import Path
from html import escape

from biocompare.core.batch import BatchResult, batch_summary
from biocompare.core.report import ConcordanceReport


def report_to_json(report: ConcordanceReport, *, indent: int = 2) -> str:
    return json.dumps(report.to_dict(), indent=indent, sort_keys=True)


def report_to_text(report: ConcordanceReport) -> str:
    lines = [
        f"Comparator: {report.comparator}",
        f"File A: {report.file_a}",
        f"File B: {report.file_b}",
        f"Overall concordance: {report.overall_concordance:.4f}",
        "Metrics:",
    ]
    for name, value in sorted(report.metrics.items()):
        lines.append(f"  {name}: {value:.4f}")
    if report.warnings:
        lines.append("Warnings:")
        lines.extend(f"  {warning}" for warning in report.warnings)
    return "\n".join(lines)


def report_to_html(report: ConcordanceReport) -> str:
    warnings = (
        "<ul>" + "".join(f"<li>{html_escape(warning)}</li>" for warning in report.warnings) + "</ul>"
        if report.warnings
        else "<p>None</p>"
    )
    metrics_rows = "".join(
        f"<tr><th>{html_escape(name)}</th><td>{value:.6f}</td></tr>"
        for name, value in sorted(report.metrics.items())
    )
    detail_rows = "".join(
        f"<tr><th>{html_escape(name)}</th><td><code>{html_escape(json.dumps(value, sort_keys=True))}</code></td></tr>"
        for name, value in sorted(report.details.items())
    )
    return html_document(
        "biocompare report",
        f"""
        <header>
          <p class="eyebrow">Single Comparison</p>
          <h1>{html_escape(report.comparator)}</h1>
          <p class="score">Overall concordance <strong>{report.overall_concordance:.4f}</strong></p>
        </header>
        <main>
          <section>
            <h2>Inputs</h2>
            <table>
              <tr><th>File A</th><td>{html_escape(report.file_a)}</td></tr>
              <tr><th>File B</th><td>{html_escape(report.file_b)}</td></tr>
            </table>
          </section>
          <section>
            <h2>Metrics</h2>
            <table>{metrics_rows}</table>
          </section>
          <section>
            <h2>Warnings</h2>
            {warnings}
          </section>
          <section>
            <h2>Details</h2>
            <table>{detail_rows}</table>
          </section>
        </main>
        """,
    )


def write_report(report: ConcordanceReport, path: str, *, fmt: str | None = None) -> None:
    output_path = Path(path)
    selected_format = fmt or output_path.suffix.lower().lstrip(".") or "json"
    if selected_format == "json":
        output_path.write_text(report_to_json(report) + "\n", encoding="utf-8")
        return
    if selected_format in {"txt", "text"}:
        output_path.write_text(report_to_text(report) + "\n", encoding="utf-8")
        return
    if selected_format in {"html", "htm"}:
        output_path.write_text(report_to_html(report) + "\n", encoding="utf-8")
        return
    raise ValueError(f"Unsupported report format {selected_format!r}")


def batch_to_json(results: list[BatchResult], *, indent: int = 2) -> str:
    payload = {
        "summary": batch_summary(results),
        "results": [
            {
                "label": result.item.label,
                "file_a": result.item.file_a,
                "file_b": result.item.file_b,
                "type": result.item.file_type,
                "status": result.status,
                "error": result.error,
                "report": result.report.to_dict() if result.report is not None else None,
            }
            for result in results
        ],
    }
    return json.dumps(payload, indent=indent, sort_keys=True)


def batch_to_tsv(results: list[BatchResult]) -> str:
    lines = [
        "\t".join(
            [
                "label",
                "status",
                "comparator",
                "overall_concordance",
                "file_a",
                "file_b",
                "error",
                "warnings",
            ]
        )
    ]
    for result in results:
        report = result.report
        fields = [
            result.item.label,
            result.status,
            report.comparator if report is not None else "",
            f"{report.overall_concordance:.6f}" if report is not None else "",
            result.item.file_a,
            result.item.file_b,
            result.error or "",
            " | ".join(report.warnings) if report is not None else "",
        ]
        lines.append("\t".join(escape_tsv_field(field) for field in fields))
    return "\n".join(lines)


def batch_to_text(results: list[BatchResult]) -> str:
    summary = batch_summary(results)
    lines = [
        f"Batch comparisons: {summary['total']}",
        f"OK: {summary['ok']}",
        f"Errors: {summary['error']}",
    ]
    if summary["mean_concordance"] is not None:
        lines.append(f"Mean concordance: {summary['mean_concordance']:.4f}")
        lines.append(f"Min concordance: {summary['min_concordance']:.4f}")
        lines.append(f"Max concordance: {summary['max_concordance']:.4f}")
    lines.append("Results:")
    for result in results:
        if result.report is None:
            lines.append(f"  {result.item.label}: error - {result.error}")
        else:
            lines.append(
                f"  {result.item.label}: {result.report.comparator} "
                f"{result.report.overall_concordance:.4f}"
            )
    return "\n".join(lines)


def batch_to_html(results: list[BatchResult]) -> str:
    summary = batch_summary(results)
    summary_rows = "".join(
        f"<tr><th>{html_escape(format_summary_label(key))}</th><td>{format_summary_value(value)}</td></tr>"
        for key, value in summary.items()
    )
    result_rows = "".join(batch_result_html_row(result) for result in results)
    return html_document(
        "biocompare batch report",
        f"""
        <header>
          <p class="eyebrow">Batch Comparison</p>
          <h1>biocompare batch report</h1>
          <p class="score">{summary['ok']} of {summary['total']} comparisons completed</p>
        </header>
        <main>
          <section>
            <h2>Summary</h2>
            <table>{summary_rows}</table>
          </section>
          <section>
            <h2>Results</h2>
            <table>
              <thead>
                <tr>
                  <th>Label</th>
                  <th>Status</th>
                  <th>Comparator</th>
                  <th>Concordance</th>
                  <th>File A</th>
                  <th>File B</th>
                  <th>Warnings / Error</th>
                </tr>
              </thead>
              <tbody>{result_rows}</tbody>
            </table>
          </section>
        </main>
        """,
    )


def write_batch(results: list[BatchResult], path: str, *, fmt: str | None = None) -> None:
    output_path = Path(path)
    selected_format = fmt or output_path.suffix.lower().lstrip(".") or "tsv"
    if selected_format == "json":
        output_path.write_text(batch_to_json(results) + "\n", encoding="utf-8")
        return
    if selected_format in {"tsv", "tab"}:
        output_path.write_text(batch_to_tsv(results) + "\n", encoding="utf-8")
        return
    if selected_format in {"txt", "text"}:
        output_path.write_text(batch_to_text(results) + "\n", encoding="utf-8")
        return
    if selected_format in {"html", "htm"}:
        output_path.write_text(batch_to_html(results) + "\n", encoding="utf-8")
        return
    raise ValueError(f"Unsupported batch report format {selected_format!r}")


def escape_tsv_field(value: object) -> str:
    return str(value).replace("\t", " ").replace("\n", " ")


def batch_result_html_row(result: BatchResult) -> str:
    report = result.report
    comparator = report.comparator if report is not None else ""
    concordance = f"{report.overall_concordance:.4f}" if report is not None else ""
    messages = result.error or (" | ".join(report.warnings) if report is not None else "")
    status_class = "status-error" if result.error else "status-ok"
    return (
        "<tr>"
        f"<th>{html_escape(result.item.label)}</th>"
        f"<td class=\"{status_class}\">{html_escape(result.status)}</td>"
        f"<td>{html_escape(comparator)}</td>"
        f"<td>{html_escape(concordance)}</td>"
        f"<td>{html_escape(result.item.file_a)}</td>"
        f"<td>{html_escape(result.item.file_b)}</td>"
        f"<td>{html_escape(messages)}</td>"
        "</tr>"
    )


def format_summary_value(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.6f}"
    if value is None:
        return ""
    return html_escape(value)


def format_summary_label(value: str) -> str:
    return value.replace("_", " ").capitalize()


def html_escape(value: object) -> str:
    return escape(str(value), quote=True)


def html_document(title: str, body: str) -> str:
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html_escape(title)}</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f7f8fa;
      --panel: #ffffff;
      --text: #1f2933;
      --muted: #607080;
      --border: #d8dee6;
      --accent: #0f766e;
      --danger: #b42318;
    }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--text);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.5;
    }}
    header, main {{
      width: min(1120px, calc(100vw - 48px));
      margin: 0 auto;
    }}
    header {{
      padding: 40px 0 24px;
    }}
    h1, h2, p {{
      margin-top: 0;
    }}
    h1 {{
      margin-bottom: 8px;
      font-size: 32px;
      line-height: 1.15;
    }}
    h2 {{
      margin-bottom: 12px;
      font-size: 18px;
    }}
    section {{
      margin-bottom: 24px;
      padding: 20px;
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 8px;
      overflow-x: auto;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 14px;
    }}
    th, td {{
      padding: 10px 12px;
      border-bottom: 1px solid var(--border);
      text-align: left;
      vertical-align: top;
    }}
    th {{
      color: var(--muted);
      font-weight: 600;
      white-space: nowrap;
    }}
    tr:last-child th, tr:last-child td {{
      border-bottom: 0;
    }}
    code {{
      white-space: pre-wrap;
      word-break: break-word;
    }}
    .eyebrow {{
      margin-bottom: 8px;
      color: var(--accent);
      font-size: 13px;
      font-weight: 700;
      letter-spacing: .04em;
      text-transform: uppercase;
    }}
    .score {{
      color: var(--muted);
      font-size: 18px;
    }}
    .score strong, .status-ok {{
      color: var(--accent);
    }}
    .status-error {{
      color: var(--danger);
      font-weight: 700;
    }}
  </style>
</head>
<body>
{body}
</body>
</html>"""

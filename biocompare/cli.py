from __future__ import annotations

import argparse
import sys

from biocompare.core.batch import run_batch
from biocompare.core.engine import ComparisonEngine
from biocompare.io.report_writers import (
    batch_to_html,
    batch_to_json,
    batch_to_text,
    batch_to_tsv,
    report_to_html,
    report_to_json,
    report_to_text,
    write_batch,
    write_report,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="biocompare",
        description="Generate semantic concordance reports for bioinformatics outputs.",
    )
    subparsers = parser.add_subparsers(dest="command")

    compare_parser = subparsers.add_parser("compare", help="Compare one pair of files.")
    add_common_options(compare_parser)
    compare_parser.add_argument("file_a", help="First output file to compare.")
    compare_parser.add_argument("file_b", help="Second output file to compare.")
    compare_parser.add_argument("-o", "--output", help="Write the report to a file instead of stdout.")
    compare_parser.add_argument("--format", choices=["html", "json", "text"], default="json", help="Output format.")

    batch_parser = subparsers.add_parser("batch", help="Compare file pairs listed in a CSV/TSV manifest.")
    add_common_options(batch_parser)
    batch_parser.add_argument("manifest", help="CSV/TSV manifest with file_a and file_b columns.")
    batch_parser.add_argument("--min-concordance", type=float, help="Fail if any successful comparison is below this threshold.")
    batch_parser.add_argument("--stop-on-error", action="store_true", help="Stop on the first failed comparison.")
    batch_parser.add_argument("-o", "--output", help="Write the batch report to a file instead of stdout.")
    batch_parser.add_argument("--format", choices=["html", "json", "tsv", "text"], default="tsv", help="Batch output format.")

    return parser


def add_common_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-t", "--type", dest="file_type", help="Force a comparator/file type such as bam_stats, bed, counts, deg, fasta, fastq, table, csv, tsv, or vcf.")
    parser.add_argument("--key", dest="key_column", help="Column to use for row alignment.")
    parser.add_argument("--delimiter", help="Force a delimiter for tabular files.")
    parser.add_argument("--alpha", type=float, help="DEG adjusted p-value threshold. Default: 0.05.")
    parser.add_argument("--lfc-threshold", type=float, help="DEG absolute log-fold-change threshold. Default: 0.0.")
    parser.add_argument("--top-n", type=int, help="Number of top-ranked DEG genes to compare. Default: 50.")
    parser.add_argument("--gene-column", help="Gene identifier column override for DEG/count matrices.")
    parser.add_argument("--sample-columns", help="Counts comparator sample columns as a comma-separated list.")
    parser.add_argument("--min-reciprocal-overlap", type=float, help="BED interval match threshold. Default: 0.0 for any overlap.")
    parser.add_argument("--logfc-column", help="DEG log-fold-change column override.")
    parser.add_argument("--padj-column", help="DEG adjusted p-value column override.")
    parser.add_argument("--pvalue-column", help="DEG raw p-value column override when adjusted p-values are absent.")


def main(argv: list[str] | None = None) -> int:
    args_list = list(sys.argv[1:] if argv is None else argv)
    if args_list and args_list[0] not in {"compare", "batch", "-h", "--help"}:
        args_list = ["compare", *args_list]
    parser = build_parser()
    args = parser.parse_args(args_list)
    if args.command is None:
        parser.print_help()
        return 2

    try:
        if args.command == "batch":
            return run_batch_command(args)
        return run_compare_command(args)
    except Exception as exc:
        parser.exit(2, f"biocompare: error: {exc}\n")
    return 0


def run_compare_command(args: argparse.Namespace) -> int:
    engine = ComparisonEngine()
    report = engine.compare(
        args.file_a,
        args.file_b,
        **comparison_kwargs(args),
    )
    if args.output:
        write_report(report, args.output, fmt=args.format)
    elif args.format == "html":
        print(report_to_html(report))
    elif args.format == "text":
        print(report_to_text(report))
    else:
        print(report_to_json(report))
    return 0


def run_batch_command(args: argparse.Namespace) -> int:
    results = run_batch(
        args.manifest,
        stop_on_error=args.stop_on_error,
        default_file_type=args.file_type,
        **comparison_kwargs(args, include_file_type=False),
    )
    if args.output:
        write_batch(results, args.output, fmt=args.format)
    elif args.format == "html":
        print(batch_to_html(results))
    elif args.format == "json":
        print(batch_to_json(results))
    elif args.format == "text":
        print(batch_to_text(results))
    else:
        print(batch_to_tsv(results))
    if args.min_concordance is not None and not 0.0 <= args.min_concordance <= 1.0:
        raise ValueError("min-concordance must be between 0.0 and 1.0")
    below_threshold = [
        result
        for result in results
        if result.report is not None and args.min_concordance is not None and result.report.overall_concordance < args.min_concordance
    ]
    return 1 if any(result.error for result in results) or below_threshold else 0


def comparison_kwargs(args: argparse.Namespace, *, include_file_type: bool = True) -> dict[str, object]:
    kwargs: dict[str, object] = {
        "key_column": args.key_column,
        "delimiter": args.delimiter,
        "alpha": args.alpha,
        "lfc_threshold": args.lfc_threshold,
        "top_n": args.top_n,
        "gene_column": args.gene_column,
        "sample_columns": args.sample_columns,
        "min_reciprocal_overlap": args.min_reciprocal_overlap,
        "logfc_column": args.logfc_column,
        "padj_column": args.padj_column,
        "pvalue_column": args.pvalue_column,
    }
    if include_file_type:
        kwargs["file_type"] = args.file_type
    return kwargs


if __name__ == "__main__":
    raise SystemExit(main())

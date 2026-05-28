# Design

`biocompare` is intentionally small and plugin-oriented. The core package does
four things:

1. Detect common file types.
2. Select a comparator.
3. Produce a `ConcordanceReport`.
4. Write reports in machine-readable and human-readable formats.

## Comparator Contract

Every comparator subclasses `Comparator` and implements:

- `can_handle(file_a, file_b, **kwargs) -> bool`
- `compare(file_a, file_b, **kwargs) -> ConcordanceReport`

The summary score is comparator-specific but must stay in the 0-1 range.
Comparator-specific formulas are reflected in metric names and details.

## Why Dependency-Free First

The built-in comparators use the Python standard library. This keeps the tool
easy to install in CI and pipeline test environments. Heavier format support can
be added through optional dependencies or third-party entry point plugins.

## Comparator Selection

Selection order favors specific scientific comparators before generic fallback:

1. DEG and expression-like tables.
2. Count matrices.
3. Domain formats such as BED, FASTA/FASTQ, VCF, and BAM stats text.
4. Generic CSV/TSV table comparison.

Use `--type` when a file can be interpreted more than one way.

## VCF Normalization Scope

The VCF comparator splits ALT alleles and trims shared prefix/suffix bases.
It does not perform reference-based left alignment. For production variant
benchmarking, normalize upstream with tools such as `bcftools norm` and use
specialized truth-set benchmarking where appropriate.


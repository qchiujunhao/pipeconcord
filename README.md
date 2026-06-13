# biocompare

`biocompare` is a Python toolkit for comparing bioinformatics pipeline outputs
with semantic, format-aware metrics. Instead of only checking whether files are
byte-for-byte identical, it measures whether two runs agree in biologically or
analytically meaningful ways.

Project website: <https://qchiujunhao.github.io/biocompare/>

Status: alpha. The core comparison model and initial comparators are usable, but
APIs and metrics may change as more bioinformatics formats and workflows are
validated.

This repository currently implements the Phase 1 vertical slice:

- a shared `ConcordanceReport` model
- comparator registry with plugin entry point support
- file type detection for common bioinformatics and tabular formats
- a differential expression result comparator
- a count/expression matrix comparator
- a normalized expression matrix comparator
- a BED interval comparator
- a FASTA/FASTQ sequence comparator
- a lightweight VCF comparator with ALT splitting and minimal allele trimming
- a `samtools flagstat`/`samtools stats` comparator
- a generic CSV/TSV table comparator
- JSON/text report writers
- a command-line interface
- automated tests with `unittest`

## Quickstart

```bash
python3 -m biocompare compare tests/fixtures/table_a.tsv tests/fixtures/table_b.tsv --key gene_id
```

Write a report to disk:

```bash
python3 -m biocompare compare tests/fixtures/table_a.tsv tests/fixtures/table_b.tsv \
  --key gene_id \
  --output report.json
```

Compare differential expression tables:

```bash
python3 -m biocompare compare tests/fixtures/degs_tool_a.tsv tests/fixtures/degs_tool_b.tsv \
  --type deg \
  --alpha 0.05
```

Compare count matrices:

```bash
python3 -m biocompare compare tests/fixtures/counts_a.tsv tests/fixtures/counts_b.tsv \
  --type counts
```

Compare normalized expression matrices:

```bash
python3 -m biocompare compare tests/fixtures/expression_tpm_a.tsv tests/fixtures/expression_tpm_b.tsv \
  --type expression
```

Compare BED intervals:

```bash
python3 -m biocompare compare tests/fixtures/peaks_a.bed tests/fixtures/peaks_b.bed \
  --type bed \
  --min-reciprocal-overlap 0.5
```

Compare FASTA sequences:

```bash
python3 -m biocompare compare tests/fixtures/sequences_a.fa tests/fixtures/sequences_b.fa \
  --type fasta
```

Compare VCF calls:

```bash
python3 -m biocompare compare tests/fixtures/calls_a.vcf tests/fixtures/calls_b.vcf \
  --type vcf
```

Optionally provide a reference FASTA for simple repeated-indel left alignment:

```bash
python3 -m biocompare compare calls_a.vcf calls_b.vcf \
  --type vcf \
  --reference-fasta reference.fa
```

Compare alignment summary statistics:

```bash
python3 -m biocompare compare tests/fixtures/flagstat_a.txt tests/fixtures/flagstat_b.txt \
  --type bam_stats
```

Run the tests:

```bash
python3 -m unittest discover -s tests
```

Install development tooling:

```bash
python3 -m pip install -e ".[dev]"
python3 -m ruff check .
python3 -m coverage run -m unittest discover -s tests
python3 -m coverage report
```

Run a batch comparison from a manifest:

```bash
python3 -m biocompare batch tests/fixtures/batch_manifest.tsv --format text
```

The manifest must contain `file_a` and `file_b` columns. Optional columns are
`label` and `type`.

Use `--min-concordance` in CI to fail when any comparison falls below a chosen
threshold:

```bash
python3 -m biocompare batch tests/fixtures/batch_manifest.tsv --min-concordance 0.95
```

Write an HTML report:

```bash
python3 -m biocompare compare tests/fixtures/peaks_a.bed tests/fixtures/peaks_b.bed \
  --type bed \
  --format html \
  --output report.html
```

Batch reports also support `--format html`.

## Plugin Model

Comparators subclass `biocompare.comparators.base.Comparator` and return a
`biocompare.core.report.ConcordanceReport`. Third-party packages can register
comparators with the `biocompare.comparators` entry point group.

## Documentation

Additional documentation is available on the project website and in `docs/`,
including API notes, design rationale, and tutorials for regression testing,
RNA-seq outputs, variant calls, and BED peak comparisons.

## Citation and Paper Draft

Citation metadata is available in `CITATION.cff`. A draft JOSS-style paper is
available under `paper/`.

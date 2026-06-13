# PipeConcord

`pipeconcord` is a Python toolkit for comparing bioinformatics pipeline outputs
with semantic, format-aware metrics. Instead of only checking whether files are
byte-for-byte identical, it measures whether two runs agree in biologically or
analytically meaningful ways.

Project website: <https://qchiujunhao.github.io/pipeconcord/>

Status: alpha. The core comparison model and initial comparators are usable, but
APIs and metrics may change as more bioinformatics formats and workflows are
validated.

Rename note: the first alpha release used the package name `biocompare`. Current
and future releases use `pipeconcord` to avoid confusion with an unrelated
life-science product directory.

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

Install from PyPI:

```bash
python3 -m pip install pipeconcord
```

The package is published at <https://pypi.org/project/pipeconcord/>.
PipeConcord requires Python 3.10 or newer.

To test unreleased changes from the default branch, install from GitHub:

```bash
python3 -m pip install git+https://github.com/qchiujunhao/pipeconcord.git
```

Create two small TSV files:

```bash
cat > old.tsv <<'EOF'
gene_id	value
A	1.0
B	2.0
EOF

cat > new.tsv <<'EOF'
gene_id	value
A	1.1
B	2.0
EOF
```

Compare them by `gene_id`:

```bash
pipeconcord compare old.tsv new.tsv --key gene_id
```

Write a report to disk:

```bash
pipeconcord compare old.tsv new.tsv \
  --key gene_id \
  --output report.json
```

Compare differential expression tables:

```bash
pipeconcord compare old_de.tsv new_de.tsv \
  --type deg \
  --alpha 0.05
```

Compare count matrices:

```bash
pipeconcord compare old_counts.tsv new_counts.tsv \
  --type counts
```

Compare normalized expression matrices:

```bash
pipeconcord compare old_tpm.tsv new_tpm.tsv \
  --type expression
```

Compare BED intervals:

```bash
pipeconcord compare old_peaks.bed new_peaks.bed \
  --type bed \
  --min-reciprocal-overlap 0.5
```

Compare FASTA sequences:

```bash
pipeconcord compare old_sequences.fa new_sequences.fa \
  --type fasta
```

Compare VCF calls:

```bash
pipeconcord compare old_calls.vcf new_calls.vcf \
  --type vcf
```

Optionally provide a reference FASTA for simple repeated-indel left alignment:

```bash
pipeconcord compare calls_a.vcf calls_b.vcf \
  --type vcf \
  --reference-fasta reference.fa
```

Compare alignment summary statistics:

```bash
pipeconcord compare old_flagstat.txt new_flagstat.txt \
  --type bam_stats
```

Run a batch comparison from a manifest:

```bash
pipeconcord batch manifest.tsv --format text
```

The manifest must contain `file_a` and `file_b` columns. Optional columns are
`label` and `type`.

Use `--min-concordance` in CI to fail when any comparison falls below a chosen
threshold:

```bash
pipeconcord batch manifest.tsv --min-concordance 0.95
```

Write an HTML report:

```bash
pipeconcord compare old_peaks.bed new_peaks.bed \
  --type bed \
  --format html \
  --output report.html
```

Batch reports also support `--format html`.

## Development

Install the repository in editable mode with development tools:

```bash
python3 -m pip install -e ".[dev]"
```

Run the tests:

```bash
python3 -m unittest discover -s tests
```

Run lint and coverage:

```bash
python3 -m ruff check .
python3 -m coverage run -m unittest discover -s tests
python3 -m coverage report
```

## Plugin Model

Comparators subclass `pipeconcord.comparators.base.Comparator` and return a
`pipeconcord.core.report.ConcordanceReport`. Third-party packages can register
comparators with the `pipeconcord.comparators` entry point group.

## Documentation

Additional documentation is available on the project website and in `docs/`,
including API notes, design rationale, and tutorials for regression testing,
RNA-seq outputs, variant calls, and BED peak comparisons.

## Citation and Paper Draft

Citation metadata is available in `CITATION.cff`. A draft JOSS-style paper is
available under `paper/`.

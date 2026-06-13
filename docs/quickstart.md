# Quickstart

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

Compare them and align rows by `gene_id`:

```bash
pipeconcord compare old.tsv new.tsv --key gene_id
```

The JSON report includes:

- `overall_concordance`: a 0-1 summary score
- `metrics.row_overlap`: Jaccard overlap of row identifiers
- `numeric.<column>.*`: Pearson, Spearman, MAE, and similarity metrics
- `categorical.<column>.*`: exact agreement and Cohen's kappa
- `warnings`: assumptions or skipped content

Compare differential expression outputs:

```bash
pipeconcord compare old_de.tsv new_de.tsv --type deg
```

The DEG comparator auto-detects common gene, log-fold-change, and adjusted
p-value column names. It reports significant gene set overlap, up/down set
overlap, top-ranked gene overlap, logFC correlations, and direction agreement.

Compare count or expression matrices:

```bash
pipeconcord compare old_counts.tsv new_counts.tsv --type counts
```

The counts comparator expects a gene-by-sample table with one gene identifier
column and at least two numeric sample columns. It reports gene/sample overlap,
per-sample correlations, per-sample magnitude similarity, library size ratios,
gene profile correlations, and zero-pattern overlap.

Compare normalized expression matrices:

```bash
pipeconcord compare old_tpm.tsv new_tpm.tsv --type expression
```

The expression comparator is intended for TPM, FPKM, CPM, and other normalized
expression matrices. It reports gene/sample overlap, sample profile
correlations, expression magnitude similarity, distribution similarity, top
expressed gene overlap, and gene profile similarity across samples.

Compare BED intervals:

```bash
pipeconcord compare old_peaks.bed new_peaks.bed --type bed
```

The BED comparator assumes standard BED 0-based half-open coordinates. It
reports base-pair Jaccard, reciprocal coverage, interval precision/recall/F1,
and interval length/count similarity. Use `--min-reciprocal-overlap 0.5` when
intervals should count as matched only if both intervals overlap by at least
50 percent.

Compare FASTA or FASTQ sequences:

```bash
pipeconcord compare old_sequences.fa new_sequences.fa --type fasta
```

The sequence comparator aligns records by identifier and reports record overlap,
exact sequence agreement, length similarity, total-base ratio, and GC-content
similarity. FASTQ quality scores are validated but not scored yet.

Compare VCF calls:

```bash
pipeconcord compare old_calls.vcf new_calls.vcf --type vcf
```

The VCF comparator splits multiallelic ALT alleles and trims shared prefix/suffix
bases before comparing normalized `CHROM/POS/REF/ALT` keys. It reports variant
Jaccard, sample overlap, allele-specific genotype concordance, allele-frequency
correlation, and Ti/Tv ratio similarity. Reference-based left alignment is not
performed by default. Use `--reference-fasta reference.fa` to left-align simple
repeated indels against a provided reference. Normalize upstream with dedicated
VCF tooling when complex representation differences matter.

Compare alignment summary statistics:

```bash
pipeconcord compare old_flagstat.txt new_flagstat.txt --type bam_stats
```

The BAM stats comparator intentionally does not parse BAM files. It compares
text outputs from `samtools flagstat` or `samtools stats`, including total/mapped
read ratios, alignment-rate similarity, duplicate-rate similarity, proper-pair
rate similarity, insert-size ratio, average-read-length ratio, and error-rate
similarity when those fields are present.

Run a batch manifest:

```bash
pipeconcord batch manifest.tsv --format tsv
```

Manifest files can be CSV or TSV and must include `file_a` and `file_b`.
Optional `label` and `type` columns name each comparison and force comparator
selection per row.

For regression testing, add `--min-concordance 0.95`; the command exits with a
non-zero status if any successful comparison falls below that threshold.

Write HTML reports:

```bash
pipeconcord compare old_peaks.bed new_peaks.bed \
  --type bed \
  --format html \
  --output report.html

pipeconcord batch manifest.tsv \
  --format html \
  --output batch_report.html
```

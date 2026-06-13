# PipeConcord

`pipeconcord` compares bioinformatics pipeline outputs with semantic metrics
instead of byte-for-byte snapshots. Each comparator returns the same
`ConcordanceReport` shape, so one-off checks and batch regression tests can
consume results consistently.

Rename note: the first alpha release used the package name `biocompare`. Current
and future releases use `pipeconcord` to avoid confusion with an unrelated
life-science product directory.

The documentation site is intended for users evaluating the project as a tool
for regression testing, pipeline validation, and reproducible scientific
comparisons.

## Current Comparators

| Type | Inputs | Main metrics |
| --- | --- | --- |
| `table` | CSV/TSV tables | row overlap, numeric similarity, categorical agreement |
| `deg` | differential expression tables | significant set overlap, logFC correlation, direction agreement |
| `counts` | gene-by-sample count matrices | sample correlations, library ratios, zero-pattern overlap |
| `expression` | TPM/FPKM/CPM matrices | expression profile correlation, distribution similarity, top gene overlap |
| `bed` | BED3+ intervals | bp Jaccard, interval precision/recall/F1 |
| `fasta`/`fastq` | sequence records | record overlap, exact sequence agreement, GC/length similarity |
| `vcf` | VCF calls | normalized allele Jaccard, genotype concordance, AF correlation, Ti/Tv |
| `bam_stats` | `samtools flagstat` / `samtools stats` text | alignment-rate, duplicate-rate, and count-ratio similarity |

## Common Commands

Install from PyPI:

```bash
python3 -m pip install pipeconcord
```

PyPI package: <https://pypi.org/project/pipeconcord/>

```bash
pipeconcord compare file_a.tsv file_b.tsv --type table --key gene_id
pipeconcord batch manifest.tsv --min-concordance 0.95
pipeconcord compare file_a.bed file_b.bed --type bed --format html --output report.html
```

## Documentation

- [Quickstart](quickstart.md)
- [API](api.md)
- [Design](design.md)
- [Writing Comparators](writing_comparators.md)
- [Regression Testing Tutorial](tutorials/regression_testing.md)
- [RNA-seq Tutorial](tutorials/compare_rnaseq_pipelines.md)
- [Variant Caller Tutorial](tutorials/compare_variant_callers.md)
- [Peak Caller Tutorial](tutorials/compare_peak_callers.md)

## Source Repository

- [GitHub repository](https://github.com/qchiujunhao/pipeconcord)

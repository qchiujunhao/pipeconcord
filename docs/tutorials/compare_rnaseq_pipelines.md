# Comparing RNA-seq Pipeline Outputs

RNA-seq pipelines commonly produce raw count matrices, normalized expression
matrices, and differential expression tables. Compare each artifact with the
comparator that matches its scientific meaning.

## Counts

```bash
python3 -m biocompare compare old/counts.tsv new/counts.tsv --type counts
```

Use count metrics to check sample-wise correlation, library-size consistency,
gene overlap, and zero-pattern agreement.

## Normalized Expression

```bash
python3 -m biocompare compare old/tpm.tsv new/tpm.tsv --type expression
```

Use expression metrics to check profile correlation, magnitude similarity,
distribution similarity, and top expressed gene overlap.

## Differential Expression

```bash
python3 -m biocompare compare old/de.tsv new/de.tsv --type deg --alpha 0.05
```

Use DEG metrics to compare significant gene sets, log-fold-change rank
correlation, and up/down direction agreement. Keep contrast direction explicit:
flipping the reference condition flips logFC signs.


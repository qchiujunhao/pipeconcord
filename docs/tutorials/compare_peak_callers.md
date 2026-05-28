# Comparing Peak Caller Outputs

Use the BED comparator for ATAC-seq, ChIP-seq, and other interval outputs:

```bash
python3 -m biocompare compare peaks_a.bed peaks_b.bed --type bed
```

The comparator assumes standard BED 0-based half-open intervals and reports:

- base-pair Jaccard
- reciprocal bp coverage
- interval precision, recall, and F1
- interval length and count similarity

Set a reciprocal overlap threshold when small edge overlaps should not count as
matching intervals:

```bash
python3 -m biocompare compare peaks_a.bed peaks_b.bed \
  --type bed \
  --min-reciprocal-overlap 0.5
```


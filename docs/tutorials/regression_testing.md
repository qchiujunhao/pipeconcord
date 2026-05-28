# Regression Testing With Batch Mode

Batch mode is designed for pipeline regression tests where each workflow run
produces several output files.

Create a manifest:

```tsv
label	file_a	file_b	type
de_table	old/de.tsv	new/de.tsv	deg
counts	old/counts.tsv	new/counts.tsv	counts
peaks	old/peaks.bed	new/peaks.bed	bed
```

Run the batch check:

```bash
python3 -m biocompare batch manifest.tsv --format tsv --output concordance.tsv
```

Fail CI when any successful comparison falls below a threshold:

```bash
python3 -m biocompare batch manifest.tsv --min-concordance 0.95
```

Use `--format html --output concordance.html` for review artifacts attached to
CI runs.


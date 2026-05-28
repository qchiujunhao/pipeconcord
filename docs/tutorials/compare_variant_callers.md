# Comparing Variant Caller Outputs

For lightweight VCF concordance:

```bash
python3 -m biocompare compare caller_a.vcf caller_b.vcf --type vcf
```

The VCF comparator reports:

- normalized allele-level variant Jaccard
- sample overlap
- allele-specific genotype concordance
- allele-frequency correlation
- Ti/Tv ratio similarity

The built-in normalization splits multiallelic ALT records and trims shared
prefix/suffix bases. It does not do reference-based left alignment. Normalize
complex indels upstream when representation differences should not count as
discordance.

For truth-set benchmarking, use specialized tools such as `hap.py` or
`vcfeval`; `biocompare` is intended for pipeline regression and broad semantic
comparison.


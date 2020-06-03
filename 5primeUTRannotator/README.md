# 5' UTR annotator

We have upgraded this plugin. To use the new plugin, please go to https://github.com/ImperialCardioGenetics/UTRannotator

## Description

A VEP Plugin to annotate 5' UTR variants specifically, variants that create upstream AUGs (uAUGs) and disrupt stop codons of upstream open reading frames (uORFs)

Currently the script only annotates SNVs. We will update soon with a version that also annotates InDels.

## Usage

`./vep -i variants.vcf --plugin five_prime_UTR_annotator,uORF_starts_ends_GRCh37_PUBLIC.txt`

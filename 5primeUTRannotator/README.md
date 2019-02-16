# 5' UTR annotator

## Description

A VEP Plugin to annotate 5' UTR variants specifically, variants that create upstream AUGs (uAUGs) and disrupt stop codons of upstream open reading frames (uORFs)

Currently the script only annotates SNVs. We will update soon with a version that also annotates InDels.

## Usage

`./vep -i variants.vcf --plugin five_prime_UTR_annotator,uORF_starts_ends_GRCh37_PUBLIC.txt`

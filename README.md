# uORFs

Data and code relating to Whiffin *et al*. 2019 Characterising the loss-of-function impact of 5' untranslated region variants in whole genome sequence data from 15,708 individuals [https://www.biorxiv.org/content/10.1101/543504v1].
  
## Contents

### 5'UTR annotator

A VEP plugin that annotates 5'UTR variants for predicted effects on translation

NB: this script is under active development and currently only works on SNVs

### data files

(1) Files of all possible uAUG-creating and stop-removing variants used in the paper analysis

- uAUG-creating_all_possible_annotated.txt
- stop-removing_all_possible_annotated.txt

(2) Classifications for all genes according to likelihood that uORF-pertubation is a important disease mechanism

Genes were classified into nine categories according to the following logic:
- Class 0 - genes with no annotated 5’UTR on the canonical transcript
- Class 1 - genes with no possible uAUG-creating or stop-removing SNVs identified
- Class 2 - remaining genes with no possible SNVs of predicted high-impact
- Class 3 - remaining genes where the UTR has a high-confidence oORF (strong/moderate Kozak or documented in sorfs.orf) indicating creating a second would be of low-impact
- Class 4 - remaining genes where one or more identified high-impact SNVs have AF > 0.1% in gnomAD (genomes AC>15)
- Class 5 - remaining genes that are intolerant to LoF variants
- Class 8 - remaining genes curated as haploinsufficient by ClinGen or the MacArthur lab, curated as acting via a loss-of-function mechanism in DDG2P or with >=10 high-confidence Pathogenic LoF variants in ClinVar (known LoF disease genes)
- Class 7 - remaining genes intolerant to LoF variants or with >=2 high-confidence Pathogenic LoF variants in ClinVar
- Class 6 - all genes not classified into any other class

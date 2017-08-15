# TCGA Glioblastoma Gene Signature Analysis 

## Summary

Completes gene signature analysis for TCGA Provisional data. 

![Survival Curve](gbm_survival.pdf)

## Reproducibility and Data

The data can be downloaded using this link: 
http://www.cbioportal.org/study?id=gbm_tcga#summary. Alternatively, to download the data 
and reproduce all analyses: 

```
bash data/data_files.sh

Rscript surv_gene_sig.R
```
## Contact 

* About the code: Lia Harrington (lia.x.harrington.gr@dartmouth.edu)

* About the project or collaboration: Dr. Mark Israel (mark.a.israel@dartmouth.edu) 
and Damian Almiron (damian.a.almiron.gr@dartmouth.edu)
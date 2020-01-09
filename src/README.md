# Focused BH

Software implementing Focused BH with the outer nodes filter. This is a multiple testing procedure intended for tree- or DAG-structured hypotheses, such as those arising from ICD codes in a phenome-wide association study.

## Reference 

E. Katsevich, C. Sabatti, and M. Bogomolov. "Controlling FDR while highlighting selected discoveries," 2019. Available on [arXiv](https://arxiv.org/abs/1809.01792)

## Files 

* FocusedBH.R:      implementation of Focused BH for the outer nodes filter
* Test_FocusedBH.R: demonstrates the use of the Focused BH code
* aux.R:            auxiliary functions (necessary only for Test_FocusedBH.R)
* coding6.tsv:      tree-structured self-reported phenotype coding from UK Biobank

## Author

* Eugene Katsevich

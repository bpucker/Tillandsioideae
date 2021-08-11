# Tillandsioideae

## 1) Prediction of peptide sequences ##
The combination of ORFpredictor, ORFfinder, and TransDecoder are integrated as previously described [(Haak et al., 2018)](https://doi.org/10.3389/fmolb.2018.00062).

## 2) Identification of the corresponding CDS ##
The script pep2cds.py was developed and applied to get the corresponding CDS of each predicted peptide sequence.

## 3) Running BUSCO ##
A previously developed wrapper script ([Pucker and Brockington, 2018](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5360-z)) was applied to run BUSCO on all peptide sequence collections. These BUSCO results were the basis for the construction of phylogenetic trees.



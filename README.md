# UKBCaseFinder
A tool for identifying patients in UK biobank given the definition of disease phenotypes (icd10, icd9, opcs or cancer histology).

Here we combine 3 data sources: hospital in-patient episode records (hesin), death records, and cancer records to identify all 
the patients and their earlist onset date.

# Installing UKBCaseFinder

To install the development version from github (the package devtools is required):
```
library(devtools)
install_github("yeyixuan/UKBCaseFinder")
```
Once installed, the package can be loaded using:
```
library(UKBCaseFinder)
```

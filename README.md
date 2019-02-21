# PrimedRPA

This is a commanline tool for the creation and filtering of primer and exo probe sets for use in Recombinase Polymerase Amplification as described here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty701/5068159

## Prerequisites & Dependencies

The PrimedRPA tool is dependent on several packages and programs which are outlined below:

- Clustal Omega (Download from: http://www.clustal.org/omega/)
- BLAST (Download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- Python 3.6  
- Pandas 0.20.3
- sys 3.6.3
- Bio 1.70
- glob2 0.5	

## Installation

Installation steps

```
git clone https://github.com/MH-LSHTM/PrimedRPA.git
cd PrimedRPA
bash PrimedRPA_Installation.sh

```
## PrimedRPA Overview


## Introduction


## Parameters File Breakdown



## Test Run

Within the cloned repositry there is a directory named Test. Within this directory there is a pre-configurated PrimedRPA_Parameters.txt file alongside a test fasta sequence. Navigate into this directory and run the command:

```
PrimedRPA PrimedRPA_Parameters.txt
```

If all dependancies are install the program shall run and the follow output should be obtained:

```

Update this output with the final version.

```
## Additional Notes 

Clustalo - Is not required if you are using a single input fasta file only.

BLAST - Is not required if a background binding check is not necessary. 

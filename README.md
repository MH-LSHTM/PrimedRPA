# PrimedRPA

This is a commanline tool for the creation and filtering of primer and exo probe sets for use in Recombinase Polymerase Amplification as described here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty701/5068159

## Prerequisites & Dependencies

The PrimedRPA tool is dependent on several packages and programs which are outlined below:

- Clustal Omega (Download from: http://www.clustal.org/omega/)
- Python 3.6  
- Pandas 0.20.3
- sys 3.6.3
- Bio 1.70

## Installation

Installation steps

```
git clone https://github.com/MH-LSHTM/PrimedRPA.git
cd PrimedRPA

#Make file executable
chmod +x ./PrimedRPA.py

#Save path to file
Path=$(pwd)

FileName="PrimedRPA.py'"

aliastext="alias PrimedRPA='"

#Create full path with name
Fullpath=${aliastext}${Path}${FileName}

# Copy bashrc to back it up
cp ~/.bashrc ~/.bashrc_prior_PrimedRPA

# Create alias for PrimedRPA
echo $Fullpath >> ~/.bashrc

# Copy bashrc to back it up
cp ~/.bash_profile ~/.bashrc_profile_PrimedRPA

# Create alias for PrimedRPA
echo $Fullpath >> ~/.bash_profile

```
## PrimedRPA Overview

## Test Run

Within the cloned repositry there is a directory named Test. Within this directory there is a pre-configurated PrimedRPA_Parameters.txt file alongside a test fasta sequence. Navigate into this directory and run the command:

```
PrimedRPA PrimedRPA_Parameters.txt
```
If all dependancies are install the program shall run and the follow output should be obtained:

```
-------------------------------------------
----------------PrimedRPA------------------
---Finding RPA Primer and Exo Probe Sets---
-------------Higgins M et al.--------------
-------------------------------------------

Received Parameters:
Input Fasta: Input_Example_1.fasta
Homology of Conserved Target DNA: 99%
Desired Primer Length: 32
Desired Probe Length: 52
Desired Amplicon Length: 200
Minimum GC Content: 30%
Maximum GC Content: 70%
Tolerated Length (self-binding): 5
Tolerated Length (secondary structure): 5
Performing Background Binding Check: no
Tolerated Bumber of Background Binding Nucleotides: 24
Split File Size: 3,007,000
Background Check Files: ['Input_Example_1.fasta']

Potential primer and probe binding sites: 1,391
Potential primer pairs after GC content filter: 1,410
Number of primer sets after removing repeated regions: 1,312
Number of sets with viable probe: 787
Number of sets remaining after self-binding filter: 117
Number of sets after filtering by secondary structure: 1
Writing Primer and Probe set list to 'RPA_Primers_and_Probes_Set_List.xlsx'
RPA Primers & Probes Set List Generated, No Background Check Performed.

```

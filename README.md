# PrimedRPA

PrimedRPA is a commanline tool for the creation and filtering of primer and exo probe sets for use in Recombinase Polymerase Amplification as described here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty701/5068159

## Prerequisites & Dependencies

The PrimedRPA tool is dependent on several packages and programs which are outlined below:

- Clustal Omega (Download from: http://www.clustal.org/omega/)
- BLAST (Download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- Python 3.6  
- Pandas 0.20.3
- sys 3.6.3
- Bio 1.70
- glob2 0.5	


Please ensure you have all of the above dependancies pre-installed before beginning the PrimedRPA installation process. 

## Installation

Installation steps

```
git clone https://github.com/MH-LSHTM/PrimedRPA.git
cd PrimedRPA
bash PrimedRPA_Installation.sh

```
## PrimedRPA Overview

Within the PrimedRPA folder you will find the following 3 files:

PrimedPRA.py - The core python script for the PrimedRPA tool
PrimedRPA_Parameters.txt - The user editable parameters file which guides the PrimedRPA tool
PrimedRPA_Installation.sh - The installation file which creates the 'PrimedRPA' alias

Note - The PrimedRPA_Parameters file can be renamed however for the remained of this document it shall be reffered to as PrimedRPA_Parameters.txt

Through modifying the PrimedRPA_Parameters.txt file the user defines the input files, filtering parameters and if necessary background files. The command outlined below commences the run:

```

PrimedRPA PrimedRPA_Parameters.txt

```

Note - The input and background fasta files outlined by user in the PrimedRPA_Parameters.txt file must be in the same working directory as the PrimedRPA_Parameter.txt file for the tool to work properly. 

## PrimedRPA_Parameters Breakdown

In this section we will breakdown each of the user defined parameters which are necessary for the PrimedRPA tool.

Run Name - For every PrimedRPA run all of the resulting & intermediate files generated will be stored in a sub-directory called ./<Run_Name>_Results 

Input File - This should be the full name of the input fasta file. The input fasta file will be interpretted to determine if a single or multiple sequences are present. If a single sequence is present the entire sequences shall be parsed to be the sequence of interest. However, if multiple sequences are present then a Clustalo alignment shall be generated and assess to identify homologous subsequences which shall then be parsed as sequences of interest.

Homology Percentage Threshold - This threshold is necessary if the input file contains multiple sequences and will act as the cut off for subsequences which are not consevered.

Desire Primed Length - TwistDX specify that the ideal RPA primers are between 30 & 35 nucleotides

Desire Probe Length - TwistDx specify that the ideal RPA probe is between 46 to 52 nucleotides

Desire Amplicon Length - TwistDx recommend an amplicon length between 200bp - 500bp

Minimum & Maximum GC Content - The recommended GC content range is between 30%-70%

Internal-Binding Threshold - This is the tolerated length of sequences which could result in binding between primers & probe within the same set.

Secondary-Structure Threshold - This is the tolerated length of sequence which could result in the formation of secondary structure features and disrupt amplication e.g. hairpins.

Background binding Check - Here you can define if you would like to perform a background binding check to see if any off-site binding to possible background sequences would occur. Please see the section below for further detail on the background binding check.

## Test Run

Within the cloned repositry there is a directory named Test. Within this directory there is a pre-configurated PrimedRPA_Parameters.txt file alongside a test fasta sequence. Navigate into this directory and run the command:

```
PrimedRPA PrimedRPA_Parameters.txt
```

If all dependancies are install the program shall run and the follow output should be obtained:

```

Update this output with the final version.

```

## Output Summary

The choice to or to not run a background binding check will affect the output obtained. If no background binding check is run then you will recieve a single output file with the nomenclature: 

<Run_Name>_PrimedRPA_Output_Sets.xlsx

However, if you run a background binding check then the following directories & files shall be generated:

./<Run_Name>_Results
	./<Run_Name>_BLAST_DB/ - directory containing the local BLAST database created. 
	./<Run_Name>_BLAST_Outputs/ - directory containing the alignment outputs generated for each of the primers & probe sets. 
	./<Run_Name>_Oligo_Fasta/ - POSSIBLY BUILD THIS IN TO DELETE AS NOT REALLY NECESSARY
	.<Run_Name>_Aligned_Input.fasta - This is the output file generated from the clustalo alignment of the input file.
	./<Run_Name>_PrimedRPA_Output_Sets.csv - This is the final output containing the Primed RPA sets of interest. 


## Background Binding Check

The background binding check is completed using BLAST and utilises the bit-score generated when comparing each primer or probe against the background sequences. For each oligo an overall percentage alignment is calculated indicating the likelihood of off-site binding to one or more of the background sequences. If the oligo has a direct match to just one or more background binding sequences the overall percentage alignment will be fixed at 100%. 

The overall set score is taken as the mean of the overall percentage alignments for each of the oligo's which make up a given PrimedRPA primer & probe set. This overall set score is used in the filtering process utilising the cut-off percentage set by the user in the PrimedRPA_Parameters.txt file.


## Additional Notes 

Clustalo - Is not required if you are using a single input fasta file only.

BLAST - Is not required if a background binding check is not necessary. 

Bit Score - The bit score, S', is derived from the raw alignment score, S, taking the statistical properties of the scoring system into account. Because bit scores are normalized with respect to the scoring system, they can be used to compare alignment scores from different searches. - Definition taken from NCBI - https://www.ncbi.nlm.nih.gov/books/NBK62051/ 
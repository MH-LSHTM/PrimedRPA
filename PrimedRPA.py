#! /usr/bin/env python3

##########################################################################################################################
#  PrimedRPA: RPA Primer and Probe Set Finder                      														 #
#  Higgins M et al. Submitted. 2018                                 													 #
#                                                                   													 #
#  Dependencies:                                                    												     #
#     Clustal Omega (Download from: http://www.clustal.org/omega/)  													 #
#	  BLASTn (Download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)	 #
#     Python 3.6                                                    													 #
#     Pandas 0.20.3                                                 													 #
#     sys 3.6.3                                                     													 #
#     Bio 1.70                                                      													 #
#	  glob2 0.5																											 #
##########################################################################################################################

import os
import sys
import glob

print('-------------------------------------------')
print('----------------PrimedRPA------------------')
print('---Finding RPA Primer and Exo Probe Sets---')
print('-------------Higgins M et al.--------------')
print('-------------------------------------------\n')

if '-h' in sys.argv or '--help' in sys.argv or len(sys.argv) == 1:
	print('Usage: ./PrimedRPA.py <parameters_file>\n')
	print('Notes: PrimedRPA expects to find certain files in the current working directory, these include:')
	print('  - parameters_file: Runtime parameters and specific names should be stored here.')
	print('  - <input.fasta>: This is your input file, its specific name should be declared in your parameters file.')
	print('  - <background.fasta>: This is your background check file, its specific name should be specified in your parameters file.\n')
	print('Additional information can be found at pathogenseq.lshtm.ac.uk/PrimedRPA.html')
	print('If used, please cite:\n\tHiggins M et al. Submitted. 2018')
	sys.exit()

import math
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import Bio.SeqIO as SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline

parametersFile = sys.argv[1]

if os.path.isfile(parametersFile):
	# Parse parameters file, assign variables accordingly.
	paraFile = open(parametersFile,"r")
	iu = []
	for line in paraFile.readlines():
		if ">" in line:
			n = line.strip('\n')
			h = n.strip('>')
			iu.append(h)
	u = iu[1:]
	Targetsequencefile = str(u[1])
	if os.path.isfile(Targetsequencefile):
		fileInfo = os.stat(Targetsequencefile)
		fileSize = fileInfo.st_size
	else:
		print("Error: Input fasta not found. Please ensure that '{}' is in the current working directory.".format(Targetsequencefile))
		sys.exit()
	NumberofConCatseq = open(Targetsequencefile,'r').read().count('>')
	singleOrMultipleFiles = 'y' if NumberofConCatseq == 1 else 'n'
	RunName = str(u[0])

	if len(glob.glob('{}*'.format(RunName))) != 0:
		print('Error: Please choose a different Run Name in the parameters file')
		sys.exit()

	UserConservedDecimal = int(u[2])
	DesiredPrimerLength = int(u[3])
	DesiredProbeLength = int(u[4])
	AmpliconSize = int(u[5])
	MinimumGCcontentFilter = int(u[6])
	MaximumGCcontentFilter = int(u[7])
	PrimerProbeSelfComplementaryBinding = int(u[8]) #This will be done in base pairs e.g 4 bp
	Lengthofregionswhichcantoleratesecondarystructure = int(u[9])#This will be done in base pairs e.g 4 bp
	BackgroundCheckRequired = str(u[10])
	PriorBackgroundPercentageCutOff = int(u[11])
	BackgroundGenomeFiles = u[12:]
	NumberOfSeqinINputfile = NumberofConCatseq + 1

	BackgroundPercentageCutOff = int(u[11])/100

	print(('Received Parameters:\n'+
		'Input Fasta: {Targetsequencefile}\n'+ 
		'Homology of Conserved Target DNA: {UserConservedDecimal}%\n'+
		'Desired Primer Length: {DesiredPrimerLength:,}\n'+
		'Desired Probe Length: {DesiredProbeLength:,}\n'+
		'Desired Amplicon Length: {AmpliconSize:,}\n'+ 
		'Minimum GC Content: {MinimumGCcontentFilter}%\n'+
		'Maximum GC Content: {MaximumGCcontentFilter}%\n'+
		'Tolerated Length (self-binding): {PrimerProbeSelfComplementaryBinding:,}\n'+
		'Tolerated Length (secondary structure): {Lengthofregionswhichcantoleratesecondarystructure:,}\n'+
		'Performing Background Binding Check: {BackgroundCheckRequired}\n'+
		'Background Binding Percentage Cut-Off: {PriorBackgroundPercentageCutOff}\n'
		'Background Check Files: {BackgroundGenomeFiles}\n').format(**globals()))
else:
	print("Error: Parameters file not found. Please ensure that '{}'' is in the current working directory.".format(parametersFile))
	sys.exit()

if fileSize > 4000000:
	print("Input file to large, reduce its size to below 4MB")
	sys.exit()

TargetSequenceListFormat = []
TargetSequenceListFormat.append(Targetsequencefile)

defaultOutputFile = '{}_Aligned_Input.fasta'.format(RunName)

# Wrapper for clustal omega
def RunningClustalo1(ListOfFileNames, overwriteOutput=True):
	for FileName in ListOfFileNames:
		if overwriteOutput:
			OutputName = defaultOutputFile
		else:
			StrippedFileName = FileName.strip(".fasta")
			OutputName = StrippedFileName + "Output"
		command = "clustalo -i '{0}' -o {1} --outfmt=clustal".format(FileName, OutputName)
		result = subprocess.call([command], stdout=subprocess.PIPE, shell=True,)

def getComplement(seq,reverse=False,rule='N2N'):
	seqComp = ""
	for base in seq:
		base = base.upper()
		if base == "A":
			seqComp += "T"
		elif base == "T":
			seqComp += "A"
		elif base == "C":
			seqComp += "G"
		elif base == "G":
			seqComp += "C"
		elif base == "N":
			if rule == 'N2-':
				seqComp += '-'
			elif rule == 'N2N':
				seqComp += "N"
		elif base == "-":
			seqComp += "N"
	if reverse:
		return(seqComp[::-1])
	else:
		return(seqComp)

def CheckingAlignedOutputFile(Outputfile,desiredAmpliconLength,singleOrMultipleFiles,NumberOfSeqinINputfile,UserConservedPercentage):
	listofconservedblocknumbers = []
	UserConservedDecimal = (UserConservedPercentage/100) * (desiredAmpliconLength+20)
	lineofinterest = NumberOfSeqinINputfile + 1
	Openfile = open(Outputfile,"r")
	lineslist = Openfile.readlines()
	precleanedlineslist = lineslist[2:]

	newlist = []
	psclist = []

	positionofspace = precleanedlineslist[1].find(' ')
	yuouy = precleanedlineslist[1][positionofspace:].count(' ')

	for seq in precleanedlineslist:
		bugfix1 = seq[(positionofspace+yuouy):].replace('\n','')
		newlist.append(bugfix1)

		psclist.append(''.join([base.upper() for base in bugfix1 if base.upper() in ['A','C','G','T']])) # Does tolerate nucleotide wildcards, need to build this in. 

	Numberofblocks = int((len(newlist)/lineofinterest))
	ListofNumbersoflinesofinterest = [x * lineofinterest for x in range(1, Numberofblocks)]

	Subsetinterestlist = []
	for hhg in ListofNumbersoflinesofinterest:
		Subsetinterestlist.append(newlist[hhg-NumberOfSeqinINputfile])

	targetSequenceofinterest=''.join(Subsetinterestlist)
	FragmentsoftargetSequenceofintereste = [targetSequenceofinterest[i:i+(desiredAmpliconLength+20)] for i in range(len(targetSequenceofinterest)-(desiredAmpliconLength+19))]


	astriceinterestlist = []
	for whg in ListofNumbersoflinesofinterest:
		astriceinterestlist.append(newlist[whg-1])
	astriceinterestsequence=''.join(astriceinterestlist)

	listofconservedregions = []

	Fragmentsofastriceinterestsequence = [astriceinterestsequence[i:i+(desiredAmpliconLength+20)] for i in range(len(astriceinterestsequence)-(desiredAmpliconLength+19))]
	faiscounter = 0
	for fais in Fragmentsofastriceinterestsequence:
		if fais.count('*') >= UserConservedDecimal:
			if '-' not in FragmentsoftargetSequenceofintereste[faiscounter]:
				listofconservedregions.append(FragmentsoftargetSequenceofintereste[faiscounter])
		faiscounter += 1

	return listofconservedregions

def CheckingNormalOutputFile(Outputfile,desiredAmpliconLength,singleOrMultipleFiles):
	with open(Outputfile,"r") as Openfile:
		lineslist = Openfile.readlines()
		sequencestartcounter = 0
		passedregionslist, sequencestartlist = [], []
		startoffirstseq = 1
		endoffirstseq = len(lineslist)-1
		if singleOrMultipleFiles == "n":
			for line in lineslist:
				sequencestartcounter = sequencestartcounter + 1
				if ">" in line:
					sequencestartlist.append(sequencestartcounter)
		if singleOrMultipleFiles == "n":
			startoffirstseq, endoffirstseq = sequencestartlist
		Newlineslist = lineslist[startoffirstseq:(endoffirstseq-1)]

	CompleteSequenceString = ""
	for seqline in Newlineslist:
		cleanedseqline = seqline.strip("\n")
		supercleanedseqline = cleanedseqline.strip()
		CompleteSequenceString += supercleanedseqline
	Fragmentsofseq = [CompleteSequenceString[i:i+(desiredAmpliconLength+20)] for i in range(len(CompleteSequenceString)-(desiredAmpliconLength+19))]
	for segment in Fragmentsofseq:
		if len(segment) > (desiredAmpliconLength+18):
			matchingbasecounter = 0
			cleanedsegment = segment.strip("\n")
			for base in cleanedsegment:
				if base == "A" or base == "G" or base == "C" or base == "T":
					matchingbasecounter += 1 
			if matchingbasecounter > (desiredAmpliconLength+16):
				passedregionslist.append(cleanedsegment)
	return passedregionslist

if singleOrMultipleFiles == "n":
	RunningClustalo1(TargetSequenceListFormat)
	Outputfilename = defaultOutputFile
	PotentialRegionsOfInterest = CheckingAlignedOutputFile(Outputfilename,AmpliconSize,singleOrMultipleFiles,NumberOfSeqinINputfile,UserConservedDecimal)
elif singleOrMultipleFiles == "y":
	Outputfilename = Targetsequencefile
	PotentialRegionsOfInterest = CheckingNormalOutputFile(Outputfilename,AmpliconSize,singleOrMultipleFiles)

CleanedFilteredSegmentsOfInterest = []
for segments in PotentialRegionsOfInterest:
	filteredSeg = segments.replace('\n','')
	CleanedFilteredSegmentsOfInterest.append(filteredSeg)

# Generating Primer Probe Binding Site
PrimerProbeBindingSiteSets = []
for string in CleanedFilteredSegmentsOfInterest:
	maximumvalue = len(string) - AmpliconSize
	Combinedprimersprobelength = DesiredProbeLength + DesiredPrimerLength + DesiredPrimerLength
	PSSD = int((AmpliconSize-Combinedprimersprobelength)/2)
	FowardPrimerBindingsiteset = list([string[i:i+DesiredPrimerLength] for i in range(maximumvalue)])
	ProbeBindingSitesets = list([string[(i+PSSD+DesiredPrimerLength):(i+DesiredPrimerLength+DesiredProbeLength+PSSD)] for i in range(maximumvalue)])
	ReversePrimerBindingSiteSet = list([string[(i+(AmpliconSize-DesiredPrimerLength)):(i+AmpliconSize)] for i in range(maximumvalue)])
	PrimerProbesetsZip = zip(FowardPrimerBindingsiteset,ProbeBindingSitesets,ReversePrimerBindingSiteSet)
	PrimerProbesets = list(PrimerProbesetsZip)
	PrimerProbeBindingSiteSets.append(PrimerProbesets)
print("Potential primer and probe binding sites: {:,}".format(len(PrimerProbeBindingSiteSets)))

# GC content Filter
GCFilteredSets = []
for PPSet in PrimerProbeBindingSiteSets:
	for seqs in PPSet:
		SN = 0
		FPBS = seqs[0]
		PBS = seqs[1]
		RPBS = seqs[2]
		CombinedSeq = "%s%s%s" % (FPBS,PBS,RPBS)
		for base in CombinedSeq:
			if base == "G" or base == "C":
				SN = SN + 1
		if seqs not in GCFilteredSets:
			if ((MinimumGCcontentFilter/100)*((DesiredPrimerLength*2)+DesiredProbeLength)) <= SN <= ((MaximumGCcontentFilter/100)*((DesiredPrimerLength*2)+DesiredProbeLength)):
				GCFilteredSets.append(seqs)
print("Potential primer pairs after GC content filter: {:,}".format(len(GCFilteredSets)))

# Filtering out sets with repetitive regions
PriorRRFilteredPPBindingSiteList = []
repSet = ["NNN", "AAAAA", "CCCCC", "TTTTT", "GGGGG"]

for forwardprimer, probe, reverseprimer in GCFilteredSets:
	# Check for a repetitive region
	for i, repType in enumerate(repSet):
		if repType in forwardprimer:
			break
		elif repType in reverseprimer:
			break
		elif repType in probe:
			break
		
		# If we survived all potential breaks, ie. no repetitive regions were found, add the set
		if i == (len(repSet)-1):
			PriorRRFilteredPPBindingSiteList.append([forwardprimer, probe, reverseprimer])

# Remove sets with too many Ns in probe, forward primer, or reverse primer
RRFilteredPPBindingSiteList = []
Nlimit = 4
for fPrimer, probe, rPrimer in PriorRRFilteredPPBindingSiteList:
	if probe.count('N') < Nlimit and fPrimer.count('N') < Nlimit and rPrimer.count('N') < Nlimit:
		RRFilteredPPBindingSiteList.append([fPrimer, probe, rPrimer])
print("Number of primer sets after removing repeated regions: {:,}".format(len(RRFilteredPPBindingSiteList)))

# Finding Fluorophone And Quencher Sites
GoodProbeFilteredPPSets = []
for PPSet in RRFilteredPPBindingSiteList:
	ProbeBindingSeq = PPSet[1]
	proberegionlength = len(ProbeBindingSeq)
	minIndexPosition = int(DesiredProbeLength*0.45)
	maxIndexPosition = int((DesiredProbeLength*0.75))
	basenumber = 0
	for base in ProbeBindingSeq:
		if (basenumber + 4) < proberegionlength:
			basenumber += 1
			if minIndexPosition < basenumber < maxIndexPosition and base == "T": 
				if (basenumber + 2) < proberegionlength:
					if ProbeBindingSeq[(basenumber+2)] == "T" or ProbeBindingSeq[(basenumber+3)] == "T":
						if PPSet not in GoodProbeFilteredPPSets:
							GoodProbeFilteredPPSets.append(PPSet)
print("Number of sets with viable probe: {:,}".format(len(GoodProbeFilteredPPSets)))

# Correct primer and probe sets
ActualPPSeqSets = []
for FPBS, PBS, RPBS in GoodProbeFilteredPPSets:
	FinalRPBS = getComplement(RPBS,reverse=True)
	ActualPPSeqSets.append([FPBS,PBS,FinalRPBS])

# Identifying Primer Probe Self Binding Sets
SelfbindingFilteredSets = []
for ForwardPrimer, Probe, ReversePrimer in ActualPPSeqSets:
	skip = False
	ppset = [ForwardPrimer, Probe, ReversePrimer]
	FPStringSubsets = [ForwardPrimer[i:i+PrimerProbeSelfComplementaryBinding] for i in range(len(ForwardPrimer)-(PrimerProbeSelfComplementaryBinding-1))]

	ComplementaryBindingSequences = []
	for substring in FPStringSubsets:
		ComplementaryBindingSequence = getComplement(substring)
		if ComplementaryBindingSequence in ReversePrimer or ComplementaryBindingSequence in Probe:
			skip = True # We want to break an iteration of the upper loop, but can only break the current loop. A bool gate passes on that message.
			break			
	if not skip and ppset not in SelfbindingFilteredSets:
			SelfbindingFilteredSets.append(ppset)
print("Number of sets remaining after self-binding filter: {:,}".format(len(SelfbindingFilteredSets)))


# Filter based on secondary structure
SecondaryStructureFilteredPPsets = []
for ppset in SelfbindingFilteredSets:
	skip = False
	counter = 0
	completeset = ppset
	ForwardPrimer, Probe, ReversePrimer = ppset

	FPStringSubsets = set([ForwardPrimer[i:i+Lengthofregionswhichcantoleratesecondarystructure] for i in range(len(ForwardPrimer)-(Lengthofregionswhichcantoleratesecondarystructure-1))])
	RPStringSubsets = set([ReversePrimer[i:i+Lengthofregionswhichcantoleratesecondarystructure] for i in range(len(ReversePrimer)-(Lengthofregionswhichcantoleratesecondarystructure-1))])
	PStringSubsets = set([Probe[i:i+Lengthofregionswhichcantoleratesecondarystructure] for i in range(len(Probe)-(Lengthofregionswhichcantoleratesecondarystructure-1))])

	FScomplementaryseq = [getComplement(s,reverse=True) for s in FPStringSubsets]
	RScomplementaryseq = [getComplement(s,reverse=True) for s in RPStringSubsets]
	Pcomplementaryseq = [getComplement(s,reverse=True) for s in PStringSubsets]

	for seq in FScomplementaryseq:
		if seq in ForwardPrimer:
			skip = True
			break

	for seq in RScomplementaryseq:
		if seq in ReversePrimer:
			skip = True
			break

	for seq in Pcomplementaryseq:
		if seq in Probe:
			skip = True
			break

	if not skip:
		SecondaryStructureFilteredPPsets.append(ppset)

print("Number of sets after filtering by secondary structure: {:,}".format(len(SecondaryStructureFilteredPPsets)))

if BackgroundCheckRequired.lower() == "no":
	def PreppingUncheckedPPSetsForList(IndexsForEachPPSet,listposition):
		List = []
		for PPsetIndexes in IndexsForEachPPSet:
			Site = PPsetIndexes[listposition]
			List.append(Site)
		return List

	FPStartUncheckedList = PreppingUncheckedPPSetsForList(SecondaryStructureFilteredPPsets,0)
	PrStartUncheckedList = PreppingUncheckedPPSetsForList(SecondaryStructureFilteredPPsets,1)
	RPStartUncheckedList = PreppingUncheckedPPSetsForList(SecondaryStructureFilteredPPsets,2)

	Table3 = {"Forward_Primer":FPStartUncheckedList,"Probe":PrStartUncheckedList, "Reverse_Primer":RPStartUncheckedList}
	DF3 = pd.DataFrame(data=Table3)
	print("Writing Primer and Probe set list to '{}'".format("RPA_Primers_and_Probes_Set_List.xlsx"))
	writer3 = pd.ExcelWriter('{}_PrimedRPA_Output_Sets.xlsx'.format(RunName)) 
	DF3.to_excel(writer3,"Sheet1")
	writer3.save()
	print("RPA Primers & Probes Set List Generated, No Background Check Performed.")
	if singleOrMultipleFiles == "n":
		os.remove(defaultOutputFile)
	quit()


else: # Given that we have background specific functions, let's only define them if they're needed.


	def BLASTCScoreGeneration(Dataframe, referenceprimerorprobe, outputCSVSetNumber, typeoligo):

		# The first section of this function is to generate the length adjusted reference bit score
		FullReferenceID = 'Oligo_Header_'+typeoligo+str(outputCSVSetNumber) 
	
		ReferenceOligoDF = Dataframe[Dataframe[0]==FullReferenceID]
		ListConverted_ReferencePrimerDF = ReferenceOligoDF.values.tolist()

		IteminsideList = ListConverted_ReferencePrimerDF[0]
		ReferenceBitScore = IteminsideList[1]
		ReferencePrimerLength = len(referenceprimerorprobe)
		LengthAdjustedBitScore = ReferenceBitScore/ReferencePrimerLength
		#print(LengthAdjustedBitScore)

		# This generates the  DF which only contains the score associated with background sequences and no other oligo sequences
		BackgroundSeqOnlyDF = Dataframe[~Dataframe[0].str.startswith('Oligo_Header_')]
		BackgroundSeq_BitScore_List = BackgroundSeqOnlyDF[1].tolist()

		if len(BackgroundSeq_BitScore_List) != 0: # This is to account for if Blastn finds no sort of alignment between the oligo and background sequences. 

			SumOfAdjustedBitScoreDifferences = 0 

			DirectMatchFailSafe = 0 #This is to account for if a complete direct match is found.

			for BackgroundSeqBitScore in BackgroundSeq_BitScore_List:
				BackgroundSeqAdjustedBitScore = BackgroundSeqBitScore/ReferencePrimerLength
				#print(BackgroundSeqAdjustedBitScore)

				if BackgroundSeqAdjustedBitScore == LengthAdjustedBitScore:
					DirectMatchFailSafe += 1

				else:
					SumOfAdjustedBitScoreDifferences = SumOfAdjustedBitScoreDifferences + BackgroundSeqAdjustedBitScore


			if DirectMatchFailSafe == 0:

				PercentageObtained = SumOfAdjustedBitScoreDifferences/(LengthAdjustedBitScore*len(BackgroundSeq_BitScore_List))

			else:
				PercentageObtained = 1

		else:
			PercentageObtained = 0

		PercentageObtained=PercentageObtained*100
		# This returns the average percentage (as a decimal) of the oligo which aligns to the background sequence.
		return PercentageObtained


	# This checks to make sure all background fasta files indicated are present.
	FastafilesInWD = glob.glob('*.fasta')
	FafilesInWD = glob.glob('*.fa')

	AllPossibleBackgroundFilesInWD = FastafilesInWD + FafilesInWD

	MissingBackgroundFiles = []
	for sbgf in BackgroundGenomeFiles:
		if sbgf not in AllPossibleBackgroundFilesInWD:
			MissingBackgroundFiles.append(sbgf)

	if len(MissingBackgroundFiles) != 0:
		print('Missing Background Files Identified')
		print(MissingBackgroundFiles)
		sys.exit()


	OligoFastaFiles = []


	SetNumberCounter = 0
	SetNumberCounterList = []
	for PPSetBC in SecondaryStructureFilteredPPsets:
		
		SetNumberCounterList.append(SetNumberCounter)

		ForwardPrimerFastaFile = open(str(SetNumberCounter)+'_FP.fasta','w')
		ForwardPrimerFastaFile.write('>Oligo_Header_FP'+str(SetNumberCounter)+'\n')
		ForwardPrimerFastaFile.write(PPSetBC[0]+'\n')
		OligoFastaFiles.append(str(SetNumberCounter)+'_FP')
		ForwardPrimerFastaFile.close()

		ProbePrimerFastaFile = open(str(SetNumberCounter)+'_PR.fasta','w')
		ProbePrimerFastaFile.write('>Oligo_Header_PR'+str(SetNumberCounter)+'\n')
		ProbePrimerFastaFile.write(PPSetBC[1]+'\n')
		OligoFastaFiles.append(str(SetNumberCounter)+'_PR')
		ProbePrimerFastaFile.close()

		ReversePrimerFastaFile = open(str(SetNumberCounter)+'_RP.fasta','w')
		ReversePrimerFastaFile.write('>Oligo_Header_RP'+str(SetNumberCounter)+'\n')
		ReversePrimerFastaFile.write(PPSetBC[2]+'\n')
		OligoFastaFiles.append(str(SetNumberCounter)+'_RP')
		ReversePrimerFastaFile.close()

		SetNumberCounter += 1


	Appendedfastaformat = []
	for debug1 in OligoFastaFiles:
		Appendedfastaformat.append(debug1+'.fasta')


	#Append the oligo fasta files to the background file list ready to create the database 
	BackgroundandOligoGenomeFiles = BackgroundGenomeFiles + Appendedfastaformat

	MergedBackgroundFastaName = '{}_Merged_Background_Fasta.fasta'.format(RunName)

	# Create Merged Background File ready for Database creation. 
	MergedBackgroundFasta = open( MergedBackgroundFastaName,'w')
	# Loop to create concatenate all possible background files
	for bf in BackgroundandOligoGenomeFiles:

		backgroundfileraw = open(bf,'r')
		for bfrline in backgroundfileraw.readlines():
			if bfrline != '\n':
				MergedBackgroundFasta.write(bfrline)

	MergedBackgroundFasta.close()


	# Create BLAST Database
	CreateDatabaseCommnad = 'makeblastdb -in {} -parse_seqids -dbtype nucl -out PrimedRPA_BLAST_Background_Seq_Database'.format(MergedBackgroundFastaName)
	CreateDatabaseCommnadresult = subprocess.call([CreateDatabaseCommnad], stdout=subprocess.PIPE, shell=True,)


	BLAST_DB_Name = '{}_BLAST_DB'.format(RunName)
	# Move BLAST Database to dedicated directory
	os.makedirs(BLAST_DB_Name)
	mvdbcommand = 'mv PrimedRPA_BLAST_Background_Seq_Database* ./{}'.format(BLAST_DB_Name)
	mvdbcommandresult = subprocess.call([mvdbcommand], stdout=subprocess.PIPE, shell=True,)


	BLAST_Outputs_Name = '{}_BLAST_Outputs'.format(RunName)
	# Blasting all primers & probes against the database created. 
	os.makedirs(BLAST_Outputs_Name)
	for olgiofile in OligoFastaFiles:
		fastaoligofile = olgiofile + '.fasta'
		BLASTSearchCommand = 'blastn -task "blastn-short" -query ./{0} -db ./{1}/PrimedRPA_BLAST_Background_Seq_Database -out ./{2}/{3}"_output.csv" -outfmt "10 sseqid bitscore evalue"'.format(fastaoligofile, BLAST_DB_Name, BLAST_Outputs_Name, olgiofile)
		BLASTSearchCommandResult = subprocess.call([BLASTSearchCommand], stdout=subprocess.PIPE, shell=True,)


	BLASTbackgroundcheckresultsFiles = glob.glob('./{}/*_output.csv'.format(BLAST_Outputs_Name))


	ForwardPrimerList = []
	ForwardPrimerScoreList = []
	ProbeList = []
	ProbeScoreList = []
	ReversePrimerList = []
	ReversePrimerScoreList = []
	OverallScore = []



	# Begining Scoring system for primer & probe sets.
	for outputCSVSetNumber in SetNumberCounterList:

		AssociatedPrimerSet = SecondaryStructureFilteredPPsets[outputCSVSetNumber]

		FowardPrimerDataFrame = pd.read_csv('./{}/'.format(BLAST_Outputs_Name)+str(outputCSVSetNumber)+'_FP_output.csv', header=None)
		ProbeDataFrame = pd.read_csv('./{}/'.format(BLAST_Outputs_Name)+str(outputCSVSetNumber)+'_PR_output.csv', header=None)
		ReversePrimerDataFrame = pd.read_csv('./{}/'.format(BLAST_Outputs_Name)+str(outputCSVSetNumber)+'_RP_output.csv', header=None)


		# Need to now update the function to account for the change in file name system
		ForwardPrimerScore = BLASTCScoreGeneration(FowardPrimerDataFrame, AssociatedPrimerSet[0], outputCSVSetNumber, 'FP')
		ProbeScore = BLASTCScoreGeneration(ProbeDataFrame, AssociatedPrimerSet[1], outputCSVSetNumber, 'PR')
		ReversePrimerScore = BLASTCScoreGeneration(ReversePrimerDataFrame, AssociatedPrimerSet[2], outputCSVSetNumber, 'RP')

		OverallScoreSetScore = (ForwardPrimerScore + ProbeScore + ReversePrimerScore)/3 #This totals the percentages and divides by 3 to get an overall score percentage. 

		ForwardPrimerList.append(AssociatedPrimerSet[0])
		ForwardPrimerScoreList.append(ForwardPrimerScore)
		ProbeList.append(AssociatedPrimerSet[1])
		ProbeScoreList.append(ProbeScore)
		ReversePrimerList.append(AssociatedPrimerSet[2])
		ReversePrimerScoreList.append(ReversePrimerScore)
		OverallScore.append(OverallScoreSetScore)


	Olig_Fasta_Name = '{}_Oligo_Fasta'.format(RunName)
	os.makedirs(Olig_Fasta_Name)
	for tidyup in Appendedfastaformat:
		tucommand = 'mv {0} ./{1}'.format(tidyup, Olig_Fasta_Name)
		tucommandresults = subprocess.call([tucommand], stdout=subprocess.PIPE, shell=True,)

	# Final pipeline output 
	PreFinalOutputDF = pd.DataFrame(data={'Forward Primer':ForwardPrimerList, 'Reverse Primer':ReversePrimerList, 'Probe':ProbeList, 'Forward Primer Score':ForwardPrimerScoreList, 'Probe Score':ProbeScoreList, 'Reverse Primer Score':ReversePrimerScoreList ,'Overall Set Score':OverallScore}).sort_values(by=['Overall Set Score'])

	CutOffAdjustedOutput = PreFinalOutputDF[PreFinalOutputDF['Overall Set Score']<PriorBackgroundPercentageCutOff]
	
	DFShape = CutOffAdjustedOutput.shape

	if DFShape == (0,7):
		print("Number of sets after filtering by background binding: {:,} \n\n".format(DFShape[0]))


		cutofffailoption = input("Do you wish to export the PrimedRPA sets which did not pass the background binding cut-off (yes/no): ")


		if cutofffailoption == 'no':


			#Clean-up 
			AllNewFiles = glob.glob('{}*'.format(RunName))
			FinalCleanUpFolder = '{}_Results'.format(RunName)
			os.makedirs(FinalCleanUpFolder)
			for thankyou in AllNewFiles:
				thankyoucommand = 'mv {0} ./{1}'.format(thankyou, FinalCleanUpFolder)
				thankyoucommandresult = subprocess.call([thankyoucommand], stdout=subprocess.PIPE, shell=True,) 

			print('\n\nAnalysis Complete\n\n')

			sys.exit()

			
		else:

			print('\n\nCreating table containing sets above the background binding cut-off') # Potentially change this into a user operated option where they either click enter or dont.

			PreFinalOutputDF.to_csv('{}_PrimedRPA_Output_Sets.csv'.format(RunName), columns=['Overall Set Score',"Forward Primer","Reverse Primer","Probe","Forward Primer Score", "Reverse Primer Score", "Probe Score"], index=False)

			#Clean-up 
			AllNewFiles = glob.glob('{}*'.format(RunName))
			FinalCleanUpFolder = '{}_Results'.format(RunName)
			os.makedirs(FinalCleanUpFolder)
			for thankyou in AllNewFiles:
				thankyoucommand = 'mv {0} ./{1}'.format(thankyou, FinalCleanUpFolder)
				thankyoucommandresult = subprocess.call([thankyoucommand], stdout=subprocess.PIPE, shell=True,) 

			print('\n\nAnalysis Complete\n\n')
			sys.exit()


	else:


		print("Number of sets after filtering by background binding: {:,}".format(DFShape[0]))
		CutOffAdjustedOutput.to_csv('{}_PrimedRPA_Output_Sets.csv'.format(RunName), columns=['Overall Set Score',"Forward Primer","Reverse Primer","Probe","Forward Primer Score", "Reverse Primer Score", "Probe Score"], index=False)

		#Clean-up 
		AllNewFiles = glob.glob('{}*'.format(RunName))
		FinalCleanUpFolder = '{}_Results'.format(RunName)
		os.makedirs(FinalCleanUpFolder)
		for thankyou in AllNewFiles:
			thankyoucommand = 'mv {0} ./{1}'.format(thankyou, FinalCleanUpFolder)
			thankyoucommandresult = subprocess.call([thankyoucommand], stdout=subprocess.PIPE, shell=True,) 

		print('\n\nAnalysis Complete\n\n')


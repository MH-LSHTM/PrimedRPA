#! /usr/bin/env python3

#####################################################################
#  PrimedRPA: RPA Primer and Probe Set Finder                       #
#  Higgins M et al. Submitted. 2018                                 #
#                                                                   #
#  Dependencies:                                                    #
#     Clustal Omega (Download from: http://www.clustal.org/omega/)  #
#     Python 3.6                                                    #
#     Pandas 0.20.3                                                 #
#     sys 3.6.3                                                     #
#     Bio 1.70                                                      #
#####################################################################

import os
import sys

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
defaultOutputFile = 'Output.fasta'

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
	Targetsequencefile = str(u[0])
	if os.path.isfile(Targetsequencefile):
		fileInfo = os.stat(Targetsequencefile)
		fileSize = fileInfo.st_size
	else:
		print("Error: Input fasta not found. Please ensure that '{}' is in the current working directory.".format(Targetsequencefile))
		sys.exit()
	NumberofConCatseq = open(Targetsequencefile,'r').read().count('>')
	singleOrMultipleFiles = 'y' if NumberofConCatseq == 1 else 'n'
	UserConservedDecimal = int(u[1])
	DesiredPrimerLength = int(u[2])
	DesiredProbeLength = int(u[3])
	AmpliconSize = int(u[4])
	MinimumGCcontentFilter = int(u[5])
	MaximumGCcontentFilter = int(u[6])
	PrimerProbeSelfComplementaryBinding = int(u[7]) #This will be done in base pairs e.g 4 bp
	Lengthofregionswhichcantoleratesecondarystructure = int(u[8])#This will be done in base pairs e.g 4 bp
	BackgroundCheckRequired = str(u[9])
	NumberOfBackgroundMatchesWhichCanBeTolerated = int(u[10])
	SplitFileSize = int(u[11])
	BackgroundGenomeFiles = u[12:]
	NumberOfSeqinINputfile = NumberofConCatseq + 1

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
		'Tolerated Bumber of Background Binding Nucleotides: {NumberOfBackgroundMatchesWhichCanBeTolerated}\n'+
		'Split File Size: {SplitFileSize:,}\n'+
		'Background Check Files: {BackgroundGenomeFiles}\n').format(**globals()))
else:
	print("Error: Parameters file not found. Please ensure that '{}'' is in the current working directory.".format(parametersFile))
	sys.exit()

if fileSize > 4000000:
	print("Input file to large, reduce its size to below 4MB")
	sys.exit()

#print("Desired amplicon length: {}".format(AmpliconSize))

TargetSequenceListFormat = []
TargetSequenceListFormat.append(Targetsequencefile)

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
	UserConservedDecimal = UserConservedPercentage/100
	lineofinterest = NumberOfSeqinINputfile + 1
	print('Line of interest: {}'.format(lineofinterest))
	Openfile = open(Outputfile,"r")
	lineslist = Openfile.readlines()
	precleanedlineslist = lineslist[2:]

	psclist = []
	for seq in precleanedlineslist:
		positionofspace = seq.find(' ')
		bugfix1 = seq[positionofspace:]
		psclist.append(''.join([base.upper() for base in bugfix1 if base.upper() in ['A','C','G','T']]))

	Numberofblocks = int((len(precleanedlineslist)/lineofinterest))
	ListofNumbersoflinesofinterest = [x * lineofinterest for x in range(1, Numberofblocks)]

	Alignedlinesofinterest = []
	for linepositionnumberfloat in ListofNumbersoflinesofinterest:
		Alignedlinesofinterest.append(precleanedlineslist[(int(linepositionnumberfloat)-1)])

	Conservedblocknumber = []
	cpcounter = 0 
	for astricesline in Alignedlinesofinterest:
		
		alcounter = 0 
		pcleanedastricesline = astricesline.replace("\n","")
		cleanedastricesline = pcleanedastricesline.strip()
		lengthofcleanedastricesline = len(cleanedastricesline)
		numberofastricesrequired = int(lengthofcleanedastricesline*UserConservedDecimal)
		
		for position in cleanedastricesline:
			if position == "*":
				alcounter += 1

		cpcounter += 1

		if alcounter >= numberofastricesrequired:
			Conservedblocknumber.append(cpcounter)

	lengthofconservedblocknumberlist = len(Conservedblocknumber)

	uno = list(range(1,lengthofconservedblocknumberlist))

	listofconservedregions = []
	nextblockcounter = -1

	for blockcounter in Conservedblocknumber:
		conservedsequencestring = "" + psclist[(blockcounter * lineofinterest) - 1]
		nextblockcounter += 1
		for nextnumber in uno:
			checknextnumber = blockcounter + nextnumber
			nextNnextB = nextnumber + nextblockcounter
			if nextNnextB < lengthofconservedblocknumberlist:
				if checknextnumber == Conservedblocknumber[nextNnextB]:
					conservedsequencestring += psclist[(checknextnumber * lineofinterest) - 2]
			else:
				break
		listofconservedregions.append(conservedsequencestring)
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
		#print(CombinedSeq)
		for base in CombinedSeq:
			if base == "G" or base == "C":
				SN = SN + 1
		if seqs not in GCFilteredSets:
			if ((MinimumGCcontentFilter/100)*118) <= SN <= ((MaximumGCcontentFilter/100)*118):
				GCFilteredSets.append(seqs)
print("Potential primer pairs after GC content filter: {:,}".format(len(GCFilteredSets)))

# Filtering out sets with repetitive regions
PriorRRFilteredPPBindingSiteList = []
repSet = ["NNN", "AAAAAAA", "CCCCCCC", "TTTTTTT", "GGGGGGG"]

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
#print(SelfbindingFilteredSets)

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
	writer3 = pd.ExcelWriter("_RPA_Primers_and_Probes__Set_List.xlsx")
	DF3.to_excel(writer3,"Sheet1")
	writer3.save()
	print("RPA Primers & Probes Set List Generated, No Background Check Performed.")
	if singleOrMultipleFiles == "n":
		os.remove(defaultOutputFile)
	quit()

else: # Given that we have background specific functions, let's only define them if they're needed.
	def CreatingSplitFiles(CleanedBackgroundFile,RoundedUpFileSplitRequired,NumberOfLinesPerDividedFile,targetSequence,TargetSequenceName):
		with open("BSSplitFile1" + TargetSequenceName + ".fasta","w") as BSSplitFile:
			BSSplitFile.write(">PrimerOrProbe\n"+targetSequence+'\n>BackgroundSequence\n')
			seqlines = CleanedBackgroundFile[1:NumberOfLinesPerDividedFile]
			for line in seqlines:
				BSSplitFile.write(line + '\n')

		for i in range(2,RoundedUpFileSplitRequired):
			seqlines = CleanedBackgroundFile[(NumberOfLinesPerDividedFile*(i-1)):(NumberOfLinesPerDividedFile*i)]
			with open(("BSSplitFile" + str(i) + TargetSequenceName + ".fasta"),  "w") as BSSFile:
				BSSFile.write(">PrimerOrProbe\n"+targetSequence+"\n>BackgroundSequence\n")
				for seq in seqlines:
					BSSFile.write(seq + '\n')
		
	def CheckingOutPutFiles(ListOfOutputFileNames,NumberOfBackgroundMatchesWhichCanBeTolerated,NumberOfSeqinINputfile):		
		listofmatchesprimers =[]
		Lineofinterest = 4

		Numberoflinesinblock = 4
		for outputFile in ListOfOutputFileNames:
			Openfile = open(outputFile,"r")
			lineslist = Openfile.readlines()
			precleanedlineslist = lineslist[2:] ##chsne

			Numberofblocks = int(((len(precleanedlineslist))/Numberoflinesinblock))
			ListofNumbersoflinesofinterest = [x * Lineofinterest for x in range(1, Numberofblocks)]

			linesofinterest = []
			for linepositionnumberfloat in ListofNumbersoflinesofinterest:
				linesofinterest.append(precleanedlineslist[(int(linepositionnumberfloat)-1)])

			linenumber = 0
			MatchingLinesList = []
			for line in linesofinterest:
				cleanedline = line.replace("\n","")
				if linenumber < (Numberofblocks - 2):
					linenumber += 1
					nextline = linesofinterest[linenumber]
					cleannextline = nextline.replace("\n","")
					doubleLineString = line + cleannextline

					if doubleLineString.count("*") > NumberOfBackgroundMatchesWhichCanBeTolerated:
						listofmatchesprimers.append("BAD")
		return listofmatchesprimers


	def ScoringSets(ListOfOutputFileNames,NumberOfBackgroundMatchesWhichCanBeTolerated,NumberOfSeqinINputfile):		
		listofmatchesprimers =[]
		Lineofinterest = 4
		Numberoflinesinblock = 4
		SetScore = 0
		for outputFile in ListOfOutputFileNames:
			Openfile = open(outputFile,"r")
			lineslist = Openfile.readlines()
			precleanedlineslist = lineslist[2:] ##chsne

			Numberofblocks = int(((len(precleanedlineslist))/Numberoflinesinblock))
			ListofNumbersoflinesofinterest = [x * Lineofinterest for x in range(1, Numberofblocks)]

			linesofinterest = []
			for linepositionnumberfloat in ListofNumbersoflinesofinterest:
				linesofinterest.append(precleanedlineslist[(int(linepositionnumberfloat)-1)])
			linenumber = 0
			MatchingLinesList = []
			for line in linesofinterest:
				cleanedline = line.replace("\n","")
				SetScore += cleanedline.count("*")
		return SetScore

	def BackgroundBindingCheck(backgroundfile,filteredPPset,NumberOfBackgroundMatchesWhichCanBeTolerated,NumberOfSeqinINputfile,SplitFileSize):
		FinalPassedPrimerAndProbeSets = []

		for PPSet in filteredPPset:
			#print(PPSet)
			ForwardPrimer = PPSet[0]
			Probe = PPSet[1]
			ReversePrimer = PPSet[2]
			#ForwardPrimer, Probe, ReversePrimer = PPSet
			FowardPrimerName, ProbeName, ReversePrimerName = "ForwardPrimer", "Probe", "ReversePrimerName"
			fileInfo = os.stat(backgroundfile)
			fileSize = fileInfo.st_size
			FileSplitRequired = (fileSize/SplitFileSize)
			#print(FileSplitRequired)
			RoundedUpFileSplitRequired = math.ceil(FileSplitRequired)
			#print(RoundedUpFileSplitRequired)

			with open(backgroundfile,"r") as backgroundseqfile:		
				CleanedBackgroundFile = []
				for line in backgroundseqfile.readlines():
					if ">" not in line:
						CleanedBackgroundFile.append(line.strip('\n'))

			NumberOfLinesPerDividedFile = int( len(CleanedBackgroundFile) / RoundedUpFileSplitRequired )

			CreatingFilesToCompareForwardPrimers = CreatingSplitFiles(CleanedBackgroundFile,RoundedUpFileSplitRequired,NumberOfLinesPerDividedFile,ForwardPrimer,FowardPrimerName)
			CreatingFilesToCompareProbe = CreatingSplitFiles(CleanedBackgroundFile,RoundedUpFileSplitRequired,NumberOfLinesPerDividedFile,Probe,ProbeName)
			CreatingFilesToCompareReversePrimers = CreatingSplitFiles(CleanedBackgroundFile,RoundedUpFileSplitRequired,NumberOfLinesPerDividedFile,ReversePrimer,ReversePrimerName)

			ListOfFileNames = []
			for file in os.listdir():
				if "BSS" in file and file.endswith(".fasta"):
					ListOfFileNames.append(str(file))

			RunningClustalo1(ListOfFileNames, overwriteOutput=False)
			ListOfOutputFileNames = [f for f in os.listdir() if 'Output' in f and f != defaultOutputFile]
			ListOfForwardPrimerOutputFileNames = [v for v in os.listdir() if 'ForwardPrimerOutput' in v and v != defaultOutputFile]
			ListOfProbeOutputFileNames = [l for l in os.listdir() if 'ProbeOutput' in l and l != defaultOutputFile]
			ListOfReversePrimerOutputFileNames = [p for p in os.listdir() if 'ReversePrimerNameOutput' in p and p != defaultOutputFile]
			BackgroundBindingCheck = CheckingOutPutFiles(ListOfOutputFileNames,NumberOfBackgroundMatchesWhichCanBeTolerated,NumberOfSeqinINputfile)

			if len(BackgroundBindingCheck) == 0:
				FPScore = ScoringSets(ListOfForwardPrimerOutputFileNames,NumberOfBackgroundMatchesWhichCanBeTolerated,NumberOfSeqinINputfile)
				PrScore = ScoringSets(ListOfProbeOutputFileNames,NumberOfBackgroundMatchesWhichCanBeTolerated,NumberOfSeqinINputfile)
				RPScore = ScoringSets(ListOfReversePrimerOutputFileNames,NumberOfBackgroundMatchesWhichCanBeTolerated,NumberOfSeqinINputfile)
				OverallScore = ScoringSets(ListOfOutputFileNames,NumberOfBackgroundMatchesWhichCanBeTolerated,NumberOfSeqinINputfile)

				PPSet.append(FPScore)
				PPSet.append(PrScore)
				PPSet.append(RPScore)
				PPSet.append(OverallScore)
				FinalPassedPrimerAndProbeSets.append(PPSet)
				#print(FinalPassedPrimerAndProbeSets)

			for x in ListOfFileNames:
				os.remove(x)
			for y in ListOfOutputFileNames:
				os.remove(y)

		print("Background check with '{}' complete.".format(backgroundfile))
		print("Number of primer & probe sets passed: {}".format(len(FinalPassedPrimerAndProbeSets)))
		return FinalPassedPrimerAndProbeSets

	# Running background checks against multiple fasta files
	for background in BackgroundGenomeFiles:
		SecondaryStructureFilteredPPsets = BackgroundBindingCheck(background,SecondaryStructureFilteredPPsets,NumberOfBackgroundMatchesWhichCanBeTolerated,NumberOfSeqinINputfile,SplitFileSize)

	if len(SecondaryStructureFilteredPPsets) == 0:
		print("No Suitable Sets Found")
		sys.exit()

	ForwardPrimerColumn = []
	for PPsets in SecondaryStructureFilteredPPsets:
		ForwardPrimer = PPsets[0]
		ForwardPrimerColumn.append(ForwardPrimer)

	ProbeColumn = []
	for PPsets in SecondaryStructureFilteredPPsets:
		P = PPsets[1]
		ProbeColumn.append(P)

	ReversePrimerColumn = []
	for PPsets in SecondaryStructureFilteredPPsets:
		RP = PPsets[2]
		ReversePrimerColumn.append(RP)

	OverallScoreColumn = [p[6] for p in SecondaryStructureFilteredPPsets]
	FPScoreColumn = [u[3] for u in SecondaryStructureFilteredPPsets]
	PrScoreColumn = [q[4] for q in SecondaryStructureFilteredPPsets]
	RPScoreColumn = [g[5] for g in SecondaryStructureFilteredPPsets]

	if len(OverallScoreColumn) != 0:
		kjkjk = min(OverallScoreColumn)-1

	adjustedScoreColumn = []
	for mnm in OverallScoreColumn:
		gjg = mnm - kjkjk
		adjustedScoreColumn.append(gjg)

	Nameofoutput = Targetsequencefile.split('.')[0]

	DF = pd.DataFrame(data={'Forward Primer':ForwardPrimerColumn, 'Reverse Primer':ReversePrimerColumn, 'Probe':ProbeColumn, 'Forward Primer Score':FPScoreColumn, 'Probe Score':PrScoreColumn, 'Reverse Primer Score':RPScoreColumn ,'Final Score':adjustedScoreColumn}).sort_values(by=['Final Score'])

	# Write to .xlsx (.csv might be more compatible)
	print("Writing output to '{}'".format(Nameofoutput+"_Primers&Probes.xlsx"))
	writer = pd.ExcelWriter(Nameofoutput+"_Primers&Probes.xlsx")
	DF.to_excel(writer,"Sheet1")
	writer.save()

	if BackgroundCheckRequired == "yes":
		if singleOrMultipleFiles == "n":
			os.remove(defaultOutputFile)

	print("RPA Primer and Probe Table Complete!")

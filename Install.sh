
#Make file executable
chmod +x ./PrimedRPA.py

#Save path to file
Path=$(pwd)

FileName="/PrimedRPA.py'"

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


echo "Installation complete, please enter the PrimedRPA_Test directory to run the first test by typing in 'PrimedRPA PrimedRPA_Parameters.txt'"

printf "\nPlease ensure you have the following dependancies installed:\nPython 3.6\nPandas 0.20.3\nsys 3.6.3\nBio 1.70\nglob2 0.5\nClustal Omega (Download from: http://www.clustal.org/omega/)\nBLASTn (Download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)\n"


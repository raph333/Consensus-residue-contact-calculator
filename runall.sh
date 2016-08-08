# bash runall.sh data/raw_pdb_files data/ras_reference_alignment.fa 1g16 PF00071

if [ $1 == '-h'  ] || [ $1 == '--help' ] || [ $# -ne 5 ]; then
  printf 'Usage: bash %s PDB-files-directory reference-alignment reference-structure-PDB-ID Pfam-domain-ID SIFTS-file\n' $0
  printf '\npositional arguments:
  PDB-files-directory  Directory with PDB-structures for calculation of
                       residue contact network.
  reference-alignment  Alignment of the sequences of the structures for which
                       residue contact networks were created. See the
                       documentation for information about the requirements
                       of such an alignment.
  reference_structure  For all residues, the equivalent residues (PDB-numbering)
                       of the reference structure will be provided.
                       Just provide the PDB-ID of your favourite structure
                       of the dataset.
  Pfam-domain-ID       Pfam domain of interest. Important for the program to
                       know which chains to analyse in case of complex
                       structures.
  sifts_chain_pfam     SIFTS-file "pdb_chain_pfam.csv" for finding chains
                       withing the PDB-structures which contain the Pfam-
                       domain of interest.\n'
  exit 0
fi

RAW_PDB_FILES_DIR=$1  # e.g. path/to/my_pdb_files
REFERENCE_ALIGNMENT=$2  # e.g. path/to/my_alignment.fa
REFERENCE_STRUCTURE=$3  # e.g. 1g16
PFAM_DOMAIN_OF_INTEREST=$4  # e.g. PF00071
SIFTS_PDB_CHAIN_PFAM=$5  # file from https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html, update if necessary
ATOMIC_DISTANCE_CUTOFF=5  # two residues are considered to form a contact if any two atoms are witing 5 Angstrom of each other
# Set the value according to your preferences

printf 'INPUT DATA:\n'
printf 'PDB-files in directory %s\n' $RAW_PDB_FILES_DIR
printf 'reference alignment: %s\n' $REFERENCE_ALIGNMENT
printf 'referenct structure: %s\n' $REFERENCE_STRUCTURE
printf 'Pfam domain of interest: %s\n' $PFAM_DOMAIN_OF_INTEREST
printf 'File for identification of PDB-chains with Pfam domain of interest: %s\n\n' $SIFTS_PDB_CHAIN_PFAM

#printf '\nCheck input data (PDB-files, reference alignment, SIFTS-mapping-file):\n'
ipython scripts/check_data.py $RAW_PDB_FILES_DIR/ $REFERENCE_ALIGNMENT $SIFTS_PDB_CHAIN_PFAM 2> /dev/null
# 2> /dev/null: redirect output to null device (output not printed)
if [ $? -eq 0 ]; then 
	printf 'Data-check: successful.\nStarting analysis...\n'
else
	printf '\nScript check_data.py has non-zero exit status (maybe an incorrect/missing argument?).\n'
	printf 'Runall script abortet.\nPlease check the input data and restart runall.sh.\n'
	exit 1
fi

# Check if output directory already exists
if [ -d "results" ]; then
	printf '\nDirectory "%s" already exists. Do you want to remove its content and continue?\n' 'results'
	printf '(press "y" to continue or "n" to exit) '
	read CONTINUE
	if [ $CONTINUE == 'y' ]; then
		rm -rf results
	else
		printf 'Program exits.\nPlease move directory "%s" to a different location.\n' 'results'
		exit 1
	fi
fi
mkdir results


printf '\nPrepare PDB-files for residue contact calculation:\n'
printf '(extract only one chain, which contains the Pfam-domain of interest, from each input PDB-file and write it to a new PDB-file)\n'
ipython scripts/process_pdb.py $RAW_PDB_FILES_DIR $PFAM_DOMAIN_OF_INTEREST $SIFTS_PDB_CHAIN_PFAM results/processed_pdb_files 2> /dev/null
if [ $? -eq 0 ]; then 
	printf 'PDB-files prepared\n'
else
	printf '\nScript process_pdb.py has non-zero exit status (maybe an incorrect/missing argument?).\n'
	printf 'Runall script abortet.\nAfter the issue with process_pdb.py is fixed, please simply restart runall.sh.\n'
	exit 1
fi

printf '\ncalculate residue contact networks:\n'
ipython scripts/calculate_networks.py results/processed_pdb_files $ATOMIC_DISTANCE_CUTOFF
if [ $? -eq 0 ]; then 
	printf 'residue contact networks calculated and written to file.\n'
else
	printf '\nScript calculate_networks.py has non-zero exit status (maybe an incorrect/missing argument?).\n'
	printf 'Runall script abortet.\nAfter the issue with calculate_networks.py is fixed, please simply restart runall.sh.\n'
	exit 1
fi

printf '\nmap PDB-residue numbers to alignment positions:\n'
ipython scripts/map_networks.py results/processed_pdb_files $REFERENCE_ALIGNMENT $REFERENCE_STRUCTURE 2> /dev/null
if [ $? -eq 0 ]; then 
	printf 'residues mapped and mapping file written.\n'
else
	printf '\nScript map_networks.py has non-zero exit status (maybe an incorrect/missing argument?).\n'
	printf 'Runall script abortet.\nAfter the issue with map_networks.py is fixed, please simply restart runall.sh.\n'
	exit 1
fi

printf '\ncalculate consensus network using the single networks and the residue mapping file:\n'
Rscript scripts/calculate_consensus_network.R results/raw_networks.csv results/mapping.csv 2> /dev/null
if [ $? -eq 0 ]; then 
	printf 'consensus network written to file.\n'
else
	printf '\nScript calculate_consensus_network.R has non-zero exit status (maybe an incorrect/missing argument?).\n'
	printf 'Runall script abortet.\nAfter the issue with calculate_consensus_network.R is fixed, please simply restart runall.sh.\n'
	exit 1
fi

printf '\nanalysing consensus residue contact network:\n'
cd results
Rscript -e "library(knitr); knit('../scripts/analysis.Rmd')"  > /dev/null #2> /dev/null
if [ $? -eq 0 ]; then 
	printf 'Markdown file created from script analysis.Rmd.\n'
else
	printf '\nScript analysis.Rmd has non-zero exit status (maybe knitr is missing in the library?).\n'
	printf 'Runall script abortet.\nAfter the issue with analysis.Rmd is fixed, please simply restart runall.sh.\n'
	exit 1
fi
Rscript -e "library(markdown); markdownToHTML('analysis.md', 'analysis.html', options=c('use_xhml'))"
if [ $? -eq 0 ]; then 
	printf 'HTML-report created.\n'
else
	printf '\nCreation of HTML-report failed (maybe markdownToHTML is missing in library?).\n'
	printf 'Runall script abortet.\nAfter the issue is fixed, please simply restart runall.sh.\n'
	exit 1
fi
pandoc -s analysis.html -o analysis.pdf
if [ $? -eq 0 ]; then 
	printf 'PDF-report created.\n'
else
	printf '\nCreation of PDF-report failed (maybe pandoc not installed?).\n'
	printf 'Runall script abortet.\nAfter the issue is fixed, please simply restart runall.sh.\n'
	exit 1
fi

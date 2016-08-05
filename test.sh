RAW_PDB_FILES_DIR=$1
REFERENCE_ALIGNMENT=$2
REFERENCE_STRUCTURE=$3
PFAM_DOMAIN=$4
SIFTS_PDB_CHAIN_PFAM='data/pdb_chain_pfam.csv'
ATOMIC_DISTANCE_CUTOFF=5
#RESULTS_DIR='results'
#PROCESSED_PDB_FILES_DIR='results/processed_pdb_files'

printf 'INPUT DATA:\n'
printf 'PDB-files in directory %s\n' $1
printf 'reference alignment: %s\n' $2
printf 'referenct structure: %s\n' $3
printf 'Pfam domain of interest: %s\n' $4


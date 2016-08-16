# ------------------------------------------------------------------------------------------
# AUTHOR: Raphael Peer, raphael1peer@gmail.com
#
# PURPOSE:
# Individual residue contact networks are combined into a consensus residue contact network
# via a a common residue numbering system defined by a reference alignment. In short, for
# every contact present in any structure, the script determines how many other structures
# have an equivalent contact.
#
# INPUT:
# 1) 'raw' residue contact networks csv-file: For every contact, the file provided PDB-ID,
# PDB-number of the first residue and PDB-number of the second residue.
# 2) mapping csv-file: For every residue, the file provides PDB-number as well as the
# position in the reference alignment. Furthermore, the equivalent residue in the reference
# structure (i.e. same column in the reference alignment) is given.
# The two input files are written by the python scripts 'calculate_networks.py' and 
# 'map_networks.py', respectively. Further information is provided in the docstring and
# comments of those scripts.
# 3) Reference alignment (fasta file): See README-file.
#
# OUTPUT:
# 1) mapped networks csv-file: Contains the same information as the first input file with
# additional information for each contact: alignment positions, equivalent residues in the
# reference structure, amino acid, etc. All this information is provided for both residues
# forming the contact (termed A and B).
# 2) MAIN RESULT: consensus network scv-file: The file provides two alignment positions
# if at least one structure has a contact between the two residues. The number of
# structures which have an equivalent residue contact is given in the column 'contact_num'.
# The column 'conservation' gives the fraction of the number of structure which have a
# contact over the total number of structures in the dataset. The value is very small the
# contact is present in only a single structure. A contact present in all structures of
# the dataset has the value 1. The last two columns give the PDB-number of the equivalent
# residues in the reference PDB-structure. (For instance, this information can be used
# to visualize the highly conserved residue contacts ('conservation' 1 or close to 1) on
# the reference structure.)
# ------------------------------------------------------------------------------------------

# READ INPUT:
# nw <-  read.csv('../results/raw_networks.csv', stringsAsFactors=F)
# mapping <- read.csv('../results/mapping.csv', stringsAsFactors=F)
args <- commandArgs(trailingOnly = TRUE)
nw <- read.csv(args[[1]], stringsAsFactors=F)
mapping <- read.csv(args[[2]], stringsAsFactors=F)
require(Biostrings)
alignment <- readAAStringSet(args[[3]], 'fasta')


require(plyr)
nw <- ddply(nw, 'pdb_id', function(x) {  # map common numbering to network
  pdbid = unique(x$pdb_id)
  map <- subset(mapping, pdb_id == pdbid)
  map$pdb_id <- NULL
  x <- merge(x, map, by.x='res_A', by.y='resnum', all.x=T)
  x <- merge(x, map, by.x='res_B', by.y='resnum', all.x=T,  suffixes=c('_A', '_B'))
  return(x)
})


unmapped_networks <- unique(subset(nw, is.na(alignment_pos_A) | is.na(alignment_pos_B))$pdb_id)
if (length(unmapped_networks) > 0 ) {
  cat('Note: Some contacts are not mapped to an alignment position in these networks:\n',
    unique(subset(nw, is.na(alignment_pos_A) | is.na(alignment_pos_B))$pdb_id), '\n',
    'Unmapped contacts will simply be excluded form consensus network calculation.\n',
    'They should not cause a problem if the remaining dataset is large enough.')
}

nw$seq_prox <- nw$res_B - nw$res_A  # proximity in sequence between the two residues
#nw <- subset(nw, seq_prox > 4)  # exclude short range interactions (half of the contacts)
nw <- nw[,c('pdb_id', 'res_A', 'pdb_A', 'alignment_pos_A', 'aa_A', 'ref_pdb_A',
            'res_B', 'pdb_B', 'alignment_pos_B', 'aa_B', 'ref_pdb_B', 'seq_prox')]
write.csv(nw, 'results/mapped_networks.csv', row.names=F, quote=F)


nw_compact <- nw[,c('pdb_id', 'alignment_pos_A', 'alignment_pos_B')]
nw_compact <- na.omit(nw_compact)

# CHECK IF ALL STRUCTURES ARE IN ALIGNMENT:
structures <- paste0(unique(nw_compact$pdb_id), '.pdb')
not_in_alignment <- structures[! structures %in% names(alignment)]
if (length(not_in_alignment) > 0) {
  cat('\nNOTE: The following structures have no corresponding sequences in the reference
      alignment. Therefore, they are excluded from further analysis. Unless you require
      those structures in the dataset, this is not a problem. If you want those structures
      inculded in the further analysis, please add their sequences to the reference
      alignment and restart the analysis.')
}


consensus <- ddply(nw_compact, c('alignment_pos_A', 'alignment_pos_B'), function(x) {
  cons <- unique(x[,c('alignment_pos_A', 'alignment_pos_B')])
  cons$contact_num <- nrow(x)  # number of structure which have an equivalent contact
  # Fraction of structures which have an equivalent contact:
  cons$conservation <- cons$contact_num / length(unique(nw_compact$pdb_id))
  return(cons)
})

# map consensus network to PDB-positions of reference structure
alignment_to_reference <- unique(mapping[,c('alignment_pos', 'ref_pdb')])
consensus <- merge(consensus, alignment_to_reference, by.x='alignment_pos_A',
                   by.y='alignment_pos')
consensus <- merge(consensus, alignment_to_reference, by.x='alignment_pos_B',
                   by.y='alignment_pos', suffixes=c('_A', '_B'))


write.csv(consensus, 'results/consensus_network.csv', row.names=F, quote=F)

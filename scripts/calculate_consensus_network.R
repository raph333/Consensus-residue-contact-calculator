# cat(sprintf('argument1: %s\n', args[[1]]))
# cat(sprintf('argument2: %s\n', args[[2]]))
#setwd('/home/raph/resint/ras/pipeline')  # THIS LINE MUST GO, find out why working directory is not set to directory of script
# nw <-  read.csv('results/raw_networks.csv', stringsAsFactors=F)
# mapping <- read.csv('results/mapping.csv', stringsAsFactors=F)
args <- commandArgs(trailingOnly = TRUE)
nw <- read.csv(args[[1]], stringsAsFactors=F)
mapping <- read.csv(args[[2]], stringsAsFactors=F)


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
  cat('WARNING: There are some contacts not mapped to an alignment position in these networks:\n',
    unique(subset(nw, is.na(alignment_pos_A) | is.na(alignment_pos_B))$pdb_id), '\n',
    'Unmapped contacts will simply be excluded form consensus network calculation.\n',
    'They should not cause a problem if the remaining dataset is large enough.')
}

nw$seq_prox <- nw$res_B - nw$res_A  # proximity in sequence between the two contacting residues
#nw <- subset(nw, seq_prox > 4)  # exclude short range interactions (excluces about half of the contacts)
nw <- nw[,c('pdb_id', 'res_A', 'pdb_A', 'alignment_pos_A', 'aa_A', 'ref_pdb_A', 'res_B', 'pdb_B', 'alignment_pos_B', 'aa_B', 'ref_pdb_B', 'seq_prox')]
#nw <- nw[,c('pdb_id', 'res_A', 'pdb_A', 'alignment_pos_A', 'res_B', 'pdb_B', 'alignment_pos_B', 'seq_prox')]
write.csv(nw, 'results/mapped_networks.csv', row.names=F, quote=F)


require(Biostrings)
alignment <- readAAStringSet('data/ras_reference_alignment.fa', 'fasta')


nw_compact <- nw[,c('pdb_id', 'alignment_pos_A', 'alignment_pos_B')]
nw_compact <- na.omit(nw_compact)

# if (nrow(nw_compact) != nrow(na.omit(nw_compact))) {
#   cat('\nNOTE: unmapped contacts have been excluded from the analysis.\n')
#   nw_compact <- na.omit(nw_compact)
# }

# CHECK IF ALL STRUCTURES ARE IN ALIGNMENT:
check <- paste0(unique(nw_compact$pdb_id), '.pdb') %in% names(alignment)

consensus <- ddply(nw_compact, c('alignment_pos_A', 'alignment_pos_B'), function(x) {
  cons <- unique(x[,c('alignment_pos_A', 'alignment_pos_B')])
  cons$contact_num <- nrow(x)
  alignment_columns <- as.matrix(alignment)[,c(cons$alignment_pos_A, cons$alignment_pos_B)]  # columns of the two contacting residues
  colnames(alignment_columns) <- c('res_a', 'res_b')
  alignment_columns <- as.data.frame(alignment_columns)
  #cons$possible_contacts <- nrow(subset(alignment_columns, res_a != '-' & res_b != '-'))
  cons$conservation <- cons$contact_num / length(unique(nw_compact$pdb_id)) # ATTENTION: new, stricter way of calculating contact conservation
  return(cons)
})
#consensus$contact_num <- consensus$possible_contacts <- NULL

# map consensus network to PDB-positions of reference structure
alignment_to_reference <- unique(mapping[,c('alignment_pos', 'ref_pdb')])
consensus <- merge(consensus, alignment_to_reference, by.x='alignment_pos_A', by.y='alignment_pos')
consensus <- merge(consensus, alignment_to_reference, by.x='alignment_pos_B', by.y='alignment_pos', suffixes=c('_A', '_B'))


write.csv(consensus, 'results/consensus_network.csv', row.names=F, quote=F)

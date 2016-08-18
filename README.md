# Consensus residue contact calculator


## Why analyse large numbers of protein structures?

More than 100,000 protein structures have been solved experimentally and deposited in the public available protein data bank (PDB). The growing number of publicly available protein structures provides a wealth of information for structural biologists. However, analysing large numbers of structures is computationally challenging. While an impressive selection computational tools is freely available for analysis of single protein structures, applying these tools in large scale analysis of protein structures (referred to as systems structural biology) requires extensive experience in bioinformatics. The availability of easy-to-use structural biology tools would enable structural biologists to better use the ever-growing wealth of structural data in the protein data bank (PDB) - the result of millions of working hours of highly skilled scientists. With this in mind, the software at hand was created.

An alternative to the usual representation of a protein as coordinates of atoms in three dimensional space is to represent proteins as networks of contacting residues (amino acids). Residue contact networks (RCNs) define residues as nodes and contacts between residues as edges in the network. In comparison to the 3D representation, RCNs are more suited to large scale computational analysis. The software at hand uses this advantage to compare RCNs across proteins of a given protein family. Thus, it allows to extract common features from any number of related protein structures (from just a few to several hundreds).

For a good overview about the method of residue contact networks, please refer this review and the references therein:  
Lesley H. Greene. “Protein structure networks”. Briefings in Functional Genomics 11.6 (2012), pp. 469–478. DOI : 10.1093/bfgp/els039

An interesting example of how to use consensus residue contact networks is presented in this paper:  
Tilman Flock et al. “Universal allosteric mechanism for G-alpha activation by GPCRs.” Nature
524.7564 (2015), pp. 173–179. DOI : 10.1038/nature14663.


## What does this software do?

The program calculates residue contact networks for a number of related protein structures from the protein data base (PDB) provided by the user. Furthermore, for every residue contact, the software determines which fraction of structures have an equivalent contact. In order to identify common residue contacts, the protein structures provided by the user have to be related (The software will finish either way but won't find any common residue contacts if unrelated protein structures are provided). In order to identify structurally equivalent residues, the user has to provide an alignment of the sequences of the protein structures.

The user can chose a particular structure of interest. The results show, for every residue contact in the structure, the fraction of related protein structures which have an equivalent contact. This information can be used to gauge the structural importance of residue contacts thus aiding interpretation of protein structures. For instance, a residue contact present in all members of a given protein family can be expected to be important for the common fold of this protein family. Moreover, the results allow to distinguish between features unique to a given protein structure of interest and features shared across the whole protein family.

Finally, the software automatically generates a report which visualizes the results and can serve as a base for further analysis.


## Software requirements

The program was developed and has been tested on Linux. It should work on any Linux-distribution.

Before running the program, please make sure you have the following software installed on your computer:

* Python version 2.7 or newer
* IPython
* Python modules: os, sys, argparse, numpy, pandas, networkx, Bio (specifically SeqIO and PDB)
* R version 3.2.5 or newer
* R packages: plyr, ggplot2, Biostrings


## How do I get set up?

1. Download the repository as zip-file from: https://github.com/raph333/Consensus-residue-contact-calculator/releases/  
2. Unpack the zip-file:
```bash
unzip Consensus-residue-contact-calculator-x.x  # x.x being the current release number
```
3. Test with example data:
```bash
cd Consensus-residue-contact-calculator-x.x
bash runall.sh test_data/raw_pdb_files test_data/ras_reference_alignment.fa test_data/pdb_chain_pfam.csv PF00071 1g16
```
4. If the script has finished, inspect the folder 'results' to see what to expect from the analysis. In case of an error, please address follow the instructions provided by the error-message (e.g. install required software).  
5. Start the script 'runall' with your own data.   
For more detailed instructions, please run
```bash
bash runall.sh --help
```
or read the next two sections.


## Input data

For simplicity, it is recommended to create a directory 'data' in the directory of the software (Consensus-residue-contact-calculator-x.x). As an example, see the directory 'test_data'. Make sure your data-directory has the following content:

* **PDB-files**: Select protein structures of interest from a given protein family and download the corresponding PDB-files. Put them in a directory which contains no other files.
* **reference alignment**: Create an alignment of all sequences of the structures to be analysed. This alignment is used by the software to identify structurally equivalent residues. Make sure the names in the alignment are the same as the names in the names of the PDB-files described above. Also make sure to use the same exact same sequences as in the PDB-structure. If a sequence in the alignment has residues not present in the structure or vice versa, the identification of structurally equivalent residues cross-referencing with other structures will be incorrect.
If a structure has no corresponding sequence in the reference alignment, it is simply ignored in the calculation of the consensus network. While the lack of one or a few sequences in the reference alignment decreases the size of the dataset, it does not corrupt the output. The presence of additional sequences in the alignment (which do not correspond to a structure in the dataset) has no impact on the result.
* **pdb_chain_pfam.csv**: You can simply use the file from the directory 'test_data' or download the latest version from SIFTS-database: from https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html

### How to create a reference alignment

The reference alignment is required in order to identify structurally equivalent residues across different protein structures. Two residues of different structures are referred to as structurally equivalent if superimposition of the two structures also superimposes the two residues.

There are two ways of creating a reference alignment:

1. Structure alignment (e.g. with MUSTANG). Obviously, the exact same structures used for residue contact calculation have to be used for the structural alignment. Therefore, it is recommended to first run 'process_pdb.py' on the downloaded 'raw' PDB-files and use the processed PDB-files for the structural alignment (as those same files are also used for residue contact calculation).
2. Use sequence alignment (e.g. with MUSCLE or ClustalW). Again, the sequences have to be exactly the same as in the structures used for residue contact calculation. Use the script 'extract_sequences_from_structures.py' to obtain the sequences from the structures. Then, these sequences can be aligned with a standard tool such as MUSCLE or ClustalW to give the reference alignment.

If the sequences are too diverse to be aligned in an unambiguous manner, it is recommended to use structure alignment instead as structure is more conserved than sequence. For highly similar sequences, structure and sequence alignment will give the same results. In these cases, simply choose the more convenient option.

## Usage

The data provided by the user (see above) is processed in a step-wise manner by python and R-scripts in the 'scripts' directory. The bash script 'runall.sh' executes all these scripts consecutively to arrive at the results. In order to do so, runall.sh requires five arguments.

Please make sure to provide the arguments in this order:

1. path to directory containing (only) the PDB-files to be analysed (see previous section)  
2. path to reference alignment in fasta format (see previous section)  
3. path to file for Pfam-domain and PDB-chain cross-referencing from the SIFTS-database (see previous section)  
4. Pfam-domain-ID of interest (e.g. PF00071): Argument is required to automatically identify the relevant part of the structures. For instance, in complex structures, only the chain containing the Pfam-domain of interest is used for analysis. Also note that only one chain per PDB-file is used. For instance, if a PDB-file contains multiple chains which contain the Pfam-domain of interest, only the first (alphabetically) is used.  
5. reference structure PDB-ID (e.g. 1g16): The positions in the reference alignment are used as a common residue numbering system. However, in most cases researchers have a particular structure of interest. For this reason, the software also provides the PDB-residue-numbers of the equivalent residues in the structure of interest - referred to as 'reference structure'. Simply provide the PDB-ID of your most interesting structure in the data set.  

Example:
```bash
bash runall.sh path/to/pdb_files_directory path/to/reference_alignment.fa path/to/pdb_chain_pfam.csv PF00071 1g16
```

Note: Two residues are considered to form a contact if any two atoms (excluding hydrogen atoms) are within 5 Angstrom of each other. This distance cutuff is defined in the script runall.sh. However, you can easily set the cutoff according to your preferences (the relevant line in the script is highlighted by a comment in capital letters).

## Results

The script runall.sh creates a directory 'results'. All intermediate and final results are saved in this directory.

### Final results:

* **mapped_networks.csv**: All residue contact networks of all input protein structure are saved in this file. In addition to the PDB-residue numbers, reference alignment position are provided for all residues. Moreover, the structurally equivalent residues in the reference PDB-structure are also provided.
* **consensus-network.csv** (main result): The file provides two alignment positions if at least one structure has a contact between the two residues. The number of structures which have an equivalent residue contact is given in the column 'contact_num'. The column 'conservation' divides the previous column by the total number of structures in the dataset. Thus, 'conservation' shows the fraction of structures that have a contact between two given residues. The value is very small the contact is only present in a single structure. A contact present in all structures of the dataset has the value 1. The last two column gives the PDB-number of the equivalent residues in the reference PDB-structure. (For instance, this information can be used to visualize the highly conserved residue contacts ('conservation' 1 or close to 1) on the reference structure.)
* **report** (analysis.HTML, analysis.PDF): An automatically created report is provided in both HTML (the images are in the folder 'figure') and PDF format. The report allows to quickly check the success of the analysis and explains the results. Moreover, the script 'analysis.Rmd' can serve as a convenient starting point for further analysis.

### Intermediate results (not important unless you modify the software):

* processed_pdb_files: Directory with processed pdb-files (only one chain per file, no heteroatoms, no hydrogen atoms etc.). See docstring of script process_pdb.py for more information.
* selected_chains_info.csv: Log file which specifies which chain has been selected for analysis from each PDB-file.
* raw_networks.csv: File containing all residue contact networks. See docstring of calculate_networks.py for more information.
* mapping.csv: File for cross-referencing PDB-residue-numbers and alignment positions in all structures. See docstring of map_networks.py for more information.
* analysis.md: Markdown-file for automatic creation of an HTML-report.


## Runtime

The runtime increases linearly with the number of protein structures in the dataset. The runtime of the entire analysis can be expected to be about three minutes for every 50 protein structures in the dataset. For large protein structures (more than 200 amino acids per relevant chain), the runtime will be longer.


## Can I use only parts of the software?

Yes. The software consists of several python and R scripts which conduct separate tasks. All these scripts are executed consecutively by the bash script 'runall.sh' for the user's convenience. However, the scripts can also be used individually. The purpose and input requirements of each script is explained in comments therein.

For instance, it is possible calculate the residue contact networks using another tool and use the rest of the software to calculate a consensus residue contact network.


## Who do I talk to?

For questions regarding the usage of the software, bug-reports, suggestions or feedback please don't hesitate to write an e-mail to raphael1peer@gmail.com.

---
output: html_document
---

# Consensus residue contact calculator


## What does this software do?

The program calculates residue contact networks for a number of related protein structures from the protein data base (PDB) provided by the user. Furthermore, for every residue contact, the software determines which fraction of structures have an equivalent contact. In order to identify common residue contacts, the protein structures provided by the user have to be related (The software will finish either way but won't find any common residue contacts if unrelated protein structures are provided). In order to identify structurally equivalent residues, the user has to provide an alignment of the sequences of the protein structures.

The user can chose a particular structure of interest. The results show, for every residue contact in the structure, the fraction of related protein structures which have an equivalent contact. This information can be used to gauge the structural importance of residue contacts thus aiding interpretation of protein structures. For instance, a residue contact present in all members of a given protein family can be expected to be important for the common fold of this protein family. Moreover, the results allow to distinguish between features unique to a given protein structure of interest and features shared across the whole protein family.

Finally, the software automatically generates a report which visualizes the results and can serve as a base for further analysis.

For further information about the method of consensus residue contact networks, please refer to the documentation and the papers cited therein.  **insert reference(s)**


## Software requirements

The program was developed and has been tested on Linux. It should work on any Linux-distribution.

Before running the program, please make sure you have the following dependencies installed on your computer:

* Python version 2.7 or newer
* IPython
* Python modules: os, sys, argparse, numpy, pandas, networkx, Bio (specifically SeqIO and PDB)
* R version 3.2.5 or newer
* R packages: plyr, ggplot2, Biostrings


## How do I get set up?

1) Download the repository as zip-file from: https://github.com/raph333/Consensus-residue-contact-calculator/releases/
2) Unpack the zip-file:  
unzip Consensus-residue-contact-calculator-x.x  # x.x being the current release number
3) Test with example data:  
cd Consensus-residue-contact-calculator-x.x  # x.x being the current release number  
bash runall.sh test_data/raw_pdb_files test_data/ras_reference_alignment.fa 1g16 PF00071
4) If the script has finished, inspect the folder 'results' to see what to expect from the analysis. In case of an error, please address follow the instructions provided by the error-message (e.g. install required software).
5) Start the script 'runall' with your own data.   
For more detailed instructions, please run  
bash runall.sh --help  
or read the next two sections.


## Input data

For simplicity, it is recommended to create a dircetory 'data' in the directory of the software (Consensus-residue-contact-calculator-x.x). As an example, see the directory 'test_data'. Make sure your data-directory has the following content:

* PDB-files: Select protein structures of interest from a given protein family and download the corresponding PDB-files. Put them in a directory which contains no other files.
* reference alignment: Create an alignment of all sequences of the structures to be analyised. This alignment is used by the software to identify structurally equivalent residues. Therefore, for proteins with low sequence similarity, it is recommended to use structural alignment (and get the sequence alignment form the structural alignment) rather than sequence alignment. However, most protein families have sufficient sequence similarity to use sequence alignment. Make sure the names in the alignment are the same as the names in the names of the PDB-files described above. Also make sure to use the same exact same sequences as in the PDB-structure. If a sequence in the alignment has residues not present in the structure or vice versa, the identification of structurally equivalent residues cross-referencing with other structures will be incorrect.
* pdb_chain_pfam.csv: You can simply use the file from the directory 'test_data' or download the latest version from SIFTS-database: from https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html

## Usage

The data provided by the user (see above) is processed in a step-wise manner by python and R-scripts in the 'scripts' directory. The bash script 'runall.sh' exectutes all these scripts consecutively to arrive at the results. In order to do so, runall.sh requires five arguments.

Please make sure to provide the arguments in this order:  
1) path to directory containing (only) the PDB-files to be analysed (see previous section)  
2) path to reference alignment in fasta format (see previous section)  
3) path to file for Pfam-domain and PDB-chain cross-referencing from the SIFTS-database (see previous section)  
4) Pfam-domain-ID of interest (e.g. PF00071): Argument is required to automatically identify the relevant part of the structures. For instance, in complex structures, only the chain containing the Pfam-domain of interest is used for analysis. Also note that only one chain per PDB-file is used. For instance, if a PDB-file contains multiple chains which contain the Pfam-domain of interest, only the first (alphabetically) is used.  
5) Reference structure PDB-ID (e.g. 1g16): The positions in the reference alignment are used as a common residue numbering system. However, in most cases researchers have a particular structure of interest. For this reason, the software also provides the PDB-residue-numbers of the equivalent residues in the stucture of interest - referred to as 'reference structure'. Simply provide the PDB-ID of your most interesting structure in the data set.  

Note: Two residues are considered to form a contact if any two atoms (excluding hydrogen atoms) are within 5 Angstrom of each other. This distance cutuff is defined in the script runall.sh. However, you can easily set the cutoff according to your prefereces (the relevant line in the script is highlighted by a comment in capital letters).


## Who do I talk to?

For questions regarding the usage of the software, bug-reports or suggestions, please don't hesitate write an e-mail to raphae1peer@gmail.com.

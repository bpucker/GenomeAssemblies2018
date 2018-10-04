# GenomeAssemblies2018
all scripts involved in the genome assembly and investigation process



# analyze_all_gene_predictions.py

python analyze_all_gene_predictions.py \
--in <FULL_PATH_TO_INPUT_DIRECTORY> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>


Description:

This scripts compares different gene predictions within the given folder. Expected are AUGUSTUS gene prediction result files (GFF3, amino acid sequences, mRNA sequences). Sets of representative sequences are identified based on maximizing the length of the encoded peptide. General statistics about the gene prediction are calculated e.g. number of genes and average exon number. 

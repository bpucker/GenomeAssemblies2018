# GenomeAssemblies2018

Scripts involved in the genome assembly and investigation process.



# analyze_all_gene_predictions.py

python analyze_all_gene_predictions.py \
--in <FULL_PATH_TO_INPUT_DIRECTORY> \
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>


Description:

This scripts compares different gene predictions within the given folder. Expected are AUGUSTUS gene prediction result files (GFF3, amino acid sequences, mRNA sequences). Sets of representative sequences are identified based on maximizing the length of the encoded peptide. General statistics about the gene prediction are calculated e.g. number of genes and average exon number. 



# assembly_polishing_pipeline.py

python assembly_polishing_pipeline.py \
--result_dir <FULL_PATH_TO_FINAL_RESULT_DIR> \
--cluster_dir <FULL_PATH_TO_TMP_CLUSTER_DIR> \
--assembly <FULL_PATH_TO_ASSEMBLY_FILE> \
--read1 <FULL_PATH_TO_FASTQ_FILE> \
--read2 <FULL_PATH_TO_FASTQ_FILE> \
--fastq_dir <FULL_PATH_TO_DIRECTORY_WITH_RNA-Seq_READS>


This pipeline takes a raw genome sequence assembly (FASTA) and runs all processing steps to polish it. Required are numerous tools, supportive scripts, and resources:


1) adapter_clipping.py

2) adapter_file containing all illumina adapter sequences (FASTA)

3) purging_short_sequences.py

4) remove_non_plant_sequences.py
	
5) FASTA file with all genome assemblies of closely related species (white list)

6) BWA_MEM_wrapper.py
	
7) REAPR

8) sort_contigs_on_ref.py

9) genome sequence of reference genome for (super-)scaffolding (FASTA)

10) organel_cleaning.py
	
11) RepeatMasker

12) run_INFERNAL_on_cluster.py

13) reads2counts_PE_mod.py

14) AUGUSTUS_hint_wrapper.py

15) FASTA file with representative peptide sequences of a reference



# Reference
Boas Pucker, Tao Feng, Samuel Brockington. Next generation sequencing to investigate genomic diversity in Caryophyllales. bioRxiv 646133; doi: https://doi.org/10.1101/646133 

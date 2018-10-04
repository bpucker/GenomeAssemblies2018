### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python assembly_polishing_pipeline.py\n
					--result_dir <FULL_PATH_TO_FINAL_RESULT_DIR>
					--cluster_dir <FULL_PATH_TO_TMP_CLUSTER_DIR>
					--assembly <FULL_PATH_TO_ASSEMBLY_FILE>
					--read1 <FULL_PATH_TO_FASTQ_FILE>
					--read2 <FULL_PATH_TO_FASTQ_FILE>
					--fastq_dir <FULL_PATH_TO_DIRECTORY_WITH_RNA-Seq_READS>
					
					bug reports and feature requests: bpucker@cebited.uni-bielefeld.de
					"""


import os, sys, time
from datetime import datetime

# --- end of imports --- #


def remove_adapter_sequences( input_file, adapter_clipping_dir, adapter_file, script, state ):
	"""! @brief remove Illumina adapters by clipping off all sequences with certain similarity """
	
	cmd = [ "python ", script, " --assembly_file ", input_file, " --output_dir ", adapter_clipping_dir, " --adapter_file ",  adapter_file ]
	if state:
		os.popen( "".join( cmd ) )
	return "".join( cmd )


def bwa_mapping( ref_file, out_dir, fastq1, fastq2, script, state ):
	"""! @brief use BWA mem to map all reads to the assembly """
	
	cmd = [ "python ", script, " --output_dir ", out_dir, " --read1 ", fastq1, " --read2 ", fastq2, " --do_not_split --reference ", ref_file ]
	if state:
		os.popen( "".join( cmd ) )
	return "".join( cmd )


def remove_short_and_low_complexitiy_seqs( input_file, output_file, script, state ):
	"""! @brief remove all very short sequences and sequences with low complexity """
	
	cmd = [ "python ", script, " --in ", input_file, " --out ", output_file ]
	if state:
		os.popen( "".join( cmd ) )
	return "".join( cmd )


def remove_non_plant_sequence( input_file, ref_file, cluster_dir, final_dir, script, state ):
	"""! @brief run BLASTn vs. plant database and nt to remove non plant sequences from assembly """
	
	cmd = [ "python ", script, " --assembly_file ", input_file, " --plant_ref_file ", ref_file, " --tmp_cluster_dir ", cluster_dir, " --final_result_dir ", final_dir, " --active", " > ", final_dir + "non_plant_removal.doc" ]
	if state:
		os.popen( "".join( cmd ) )
	return "".join( cmd )


def run_reapr( reapr, assembly_file, bam_file, output_dir, reapr_result_file, state ):
	"""! @brief run REAPR to evaluate and correct the assembly """
	
	os.chdir( output_dir )
	cmd = [ reapr, "pipeline", assembly_file, bam_file, output_dir+"sub_dir", ">", output_dir+"reapr.log" ]
	if state:
		os.popen( " ".join( cmd ) )
		while not os.path.exists( reapr_result_file ):
			time.sleep( 60 )
		time.sleep( 90 )
	return " ".join( cmd )


def remove_organel_seqs( assembly_file, output_dir, output_file, script, state ):
	"""! @brief remove organel sequences """
	
	cmd = [ "python ", script, " --prefix ", output_dir, ' --in ', assembly_file, ' --out ', output_file ]
	if state:
		os.popen( "".join( cmd ) )
	return "".join( cmd )


def sort_on_refbeet( input_file, ref_beet, output_dir, name_mapping_table, script, state ):
	"""! @brief sort contigs on RefBeet 1.5 reference """
	
	cmd = [ "python ", script, " --contig_file ", input_file, " --ref_file ", ref_beet, " --output_dir ", output_dir, " > ", name_mapping_table ]
	if state:
		os.popen( "".join( cmd ) )
	return "".join( cmd )


def repeat_masking( repeat_masker, output_dir, input_file, state ):
	"""! @brief run RepeatMasker to masker low complexity sequences prior to gene prediction """
	
	cmd = "".join( [ repeat_masker, " -engine crossmatch -parallel 4 -s -no_is -noint -species caryophyllales -gccalc -gff -dir ", output_dir, " ", input_file ] )
	if state:
		os.popen( cmd )
	return cmd


def run_INFERNAL_on_cluster( infernal_script, output_dir, input_file, state ):
	"""! @brief run INFERNAL on cluster to process multiple sequences in parallel """
	
	cmd ="".join( [ "python ", infernal_script, " --assembly ", input_file, " --output_dir ", output_dir, " --para_jobs 50", " > ", output_dir+"my_stdout.doc" ] )
	if state:
		os.popen( cmd )
	return cmd


def run_STAR( STAR_script, tmp_cluster_dir, output_dir, fastq_dir, assembly_file, state ):
	"""! @brief run STAR to generate gene prediction hints """
	
	os.chdir( tmp_cluster_dir )
	cmd = "".join( [ "python ", STAR_script, " --fastq_file_dir ", fastq_dir, " --tmp_cluster_dir ", tmp_cluster_dir, " --result_dir ", output_dir, " --ref_genome_file ", assembly_file, " > ", output_dir+"my_stdout.doc" ] )
	if state:
		os.popen( cmd )
	return cmd


def start_AUGUSTUS_gene_prediction( AUGUSTUS_script, assembly_file, bam_file, bv_prot_file, output_dir, state ):
	"""! @brief run AUGUSTUS gene prediction with protein and RNA-Seq hints """
	
	cmd = "".join( [ "python ", AUGUSTUS_script, " --genome ", assembly_file, " --bam ", bam_file, " --prot ", bv_prot_file, " --out_dir ", output_dir, " --bam_is_sorted", " > ", output_dir+"my_stdout.doc" ] )
	if state:
		os.popen( cmd )
	return cmd


def load_config_file( config_file ):
	"""! @brief load config file to find out which parts need to be run """
	
	config = { 	'adapter_clipping': True,
						'bwa': True,
						'reapr': True,
						'short_low_complex': True,
						'remove_non_plant': True,
						'sort_refbeet': True,
						'repeat_masking': True,
						'small_rna_pred': True,
						'STAR': True,
						'gene_prediction': True,
						'organel': True
					}
	
	with open( config_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[1] == "done":
					del config[ parts[0] ]
					config.update( { parts[0]: False } )
			line = f.readline()
	return config


def main( arguments ):
	"""! @brief runs the complete pipeline """
	
	final_prefix = arguments[ arguments.index( '--result_dir' )+1 ]
	cluster_prefix = arguments[ arguments.index( '--cluster_dir' )+1 ]
	fastq_dir = arguments[ arguments.index( '--fastq_dir' )+1 ]
	config_file = arguments[ arguments.index( '--config' )+1 ]
	
	input_assembly_file = arguments[ arguments.index( '--assembly' )+1 ]
	fastq1 = arguments[ arguments.index( '--read1' )+1 ]
	fastq2 = arguments[ arguments.index( '--read2' )+1 ]
	
	if final_prefix[-1] != '/':
		final_prefix += "/"
	if cluster_prefix[-1] != '/':
		cluster_prefix += "/"
	
	if not os.path.exists( final_prefix ):
		os.makedirs( final_prefix )
	if not os.path.exists( cluster_prefix ):
		os.makedirs( cluster_prefix )
	
	# --- script locations --- #	#=> could be replaced by config file
	adapter_clipper_script = "/vol/gf-arabseq/project_Ath-Nd1/members/bpucker/ca_stuff/genome_assemblies/assembly_polishing/adapter_clipping.py"
	adapter_file = "/vol/gf-arabseq/project_Ath-Nd1/members/bpucker/ca_stuff/genome_assemblies/assembly_polishing/illumina_adaptors_collapsed.fasta"
	remove_short_and_low_complex_seq_script = "/vol/gf-arabseq/project_Ath-Nd1/members/bpucker/ca_stuff/genome_assemblies/assembly_polishing/purging_short_sequences.py"
	remove_non_plant_seqs_script = "/vol/gf-arabseq/project_Ath-Nd1/members/bpucker/ca_stuff/genome_assemblies/assembly_polishing/remove_non_plant_sequences.py"
	pan_genome_file = "/vol/cluster-data/bpucker/ca/assembly_polishing/ref_genomes/caryophyllales_pan_genome.fasta"
	bwa_wrapper_script = "/vol/cluster-data/bpucker/bin/scripts/BWA_MEM_wrapper.py"
	reapr = "/vol/biotools/bin/reapr"
	sort_on_ref_script = "/vol/gf-arabseq/project_Ath-Nd1/members/bpucker/ca_stuff/genome_assemblies/assembly_polishing/sort_contigs_on_refbeet.py"
	refbeet_file = "/prj/gf-gabibeet/data/Assemblies/RefBeet/RefBeet-1.5/RefBeet-1.5.v20160810.fa"
	organel_cleaner_script = "/vol/gf-arabseq/project_Ath-Nd1/members/bpucker/ca_stuff/genome_assemblies/assembly_polishing/organel_cleaning.py"
	
	repeat_masker = "/vol/cluster-data/bpucker/bin/RepeatMasker/RepeatMasker"
	infernal_script = "/vol/gf-arabseq/project_Ath-Nd1/members/bpucker/ca_stuff/genome_assemblies/assembly_polishing/run_INFERNAL_on_cluster.py"
	STAR_script = "/vol/cluster-data/bpucker/bv_gene_prediction/reads2counts_PE_mod.py"
	AUGUSTUS_script = "/prj/gf-resmabs/members/bpucker/scripts/AUGUSTUS_hint_wrapper.py"
	bv_prot_file = "/prj/gf-gabibeet/members/Boas/20160825_protein_comparison/BeetSet-2.genes.1408.pep.fasta"
	
	
	# --- load config file --- #
	config = load_config_file( config_file )
	
	# --- starting processing and documentaiton --- #	
	doc_file = final_prefix + "DOCUMENTATION.log"
	with open( doc_file, "w", 0 ) as doc:
		t1 = datetime.now()
		doc.write( str(t1) + '\nremoving illumina adapters ...\n' )
		
		# ---- clipping of remaining illumina adapters --- #
		adapter_clipping_dir = final_prefix + "adapter_clipping/"
		if not os.path.exists( adapter_clipping_dir ):
			os.makedirs( adapter_clipping_dir )
		cmd = remove_adapter_sequences( input_assembly_file, adapter_clipping_dir, adapter_file, adapter_clipper_script, config['adapter_clipping'] )
		doc.write( cmd + '\n\n' )
		adapter_free_file = adapter_clipping_dir + "clean_assembly.fasta"
		t2 = datetime.now()
		doc.write( str( t2-t1 ) + '\n\nBWA mapping ...\n' )
		
		
		# --- BWA mapping for REAPR --- #
		out_dir = cluster_prefix + "bwa_mapping/"
		if not os.path.exists( out_dir ):
			os.makedirs( out_dir )
		cmd = bwa_mapping( adapter_free_file, out_dir, fastq1, fastq2, bwa_wrapper_script, config['bwa'] )
		doc.write( cmd + '\n\n' )
		bwa_result_file = out_dir + "final_bam_file_sorted_duplicates_removed_sorted_rg.bam.gz"
		t3 = datetime.now()
		doc.write( str( t3-t2 ) + '\n\nREAPR ...\n' )
		
		
		# --- run REAPR --- #
		output_dir = cluster_prefix + "reapr/"
		if not os.path.exists( output_dir ):
			os.makedirs( output_dir )
		reapr_result_file = output_dir + "sub_dir/04.break.broken_assembly.fa"
		cmd = run_reapr( reapr, adapter_free_file, bwa_result_file, output_dir, reapr_result_file, config['reapr'] )
		doc.write( cmd + '\n\n' )
		
		t4 = datetime.now()
		doc.write( str( t4-t3 ) + '\n\nShort and low complexity sequence removal ...\n' )
		
		
		# --- remove short and low complexity sequences --- #
		short_and_low_complex_dir = final_prefix + "short_and_low_complex/"
		if not os.path.exists( short_and_low_complex_dir ):
			os.makedirs( short_and_low_complex_dir )
		short_and_low_complex_clean_file = short_and_low_complex_dir + "short_and_low_complex_clean.fasta"
		cmd = remove_short_and_low_complexitiy_seqs( reapr_result_file, short_and_low_complex_clean_file, remove_short_and_low_complex_seq_script, config['short_low_complex'] )	
		doc.write( cmd + '\n\n' )
		t5 = datetime.now()
		doc.write( str( t5-t4 ) + '\n\nRemoval of non-plant sequences ...\n' )
		
		
		# --- remove non-plant sequences --- #
		cluster_dir = cluster_prefix + "non_plant_removal/"
		if not os.path.exists( cluster_dir ):
			os.makedirs( cluster_dir )
		final_dir = final_prefix + "non_plant_removal/"
		if not os.path.exists( final_dir ):
			os.makedirs( final_dir )
		cmd = remove_non_plant_sequence( short_and_low_complex_clean_file, pan_genome_file, cluster_dir, final_dir, remove_non_plant_seqs_script, config['remove_non_plant']  )
		clean_assembly_file = final_dir + "cleaned_assembly_file.fasta"
		doc.write( cmd + '\n\n' )
		t6 = datetime.now()
		doc.write( str( t6-t5 ) + '\n\nRemoving organel seqs ...\n' )
		
		
		# --- remove organel sequences --- #
		output_dir = final_prefix + "organel_cleaning/"
		if not os.path.exists( output_dir ):
			os.makedirs( output_dir )
		cleaned_assembly_file = output_dir + "organel_cleaned_assembly.fasta"
		cmd = remove_organel_seqs( clean_assembly_file, output_dir, cleaned_assembly_file, organel_cleaner_script, config['organel'] )
		doc.write( cmd + '\n\n' )
		t6b = datetime.now()
		doc.write( str( t6b-t6 ) + '\n\nSorting on RefBeet 1.5 ...\n' )
		
		
		# --- sort on RefBeet 1.5 --- #
		output_dir = final_prefix + "sorting_on_refbeet/"
		if not os.path.exists( output_dir ):
			os.makedirs( output_dir )
		name_mapping_table = output_dir + "contig_name_mapping_table.txt"
		cmd= sort_on_refbeet( cleaned_assembly_file, refbeet_file, output_dir, name_mapping_table, sort_on_ref_script, config['sort_refbeet'])
		ref_sorted_contig_file = output_dir + "sorted_contigs.fasta"
		doc.write( cmd + '\n\n' )
		t7 = datetime.now()
		doc.write( str( t7-t6b ) + '\n\n' )
		

		# --- repeat masking --- #
		output_dir = final_prefix + "repeat_masking/"
		if not os.path.exists( output_dir ):
				os.makedirs( output_dir )
		input_file = output_dir + ref_sorted_contig_file.split('/')[-1]
		os.popen( "cp " + ref_sorted_contig_file + " " + input_file )
		cmd = repeat_masking( repeat_masker, output_dir, input_file, config['repeat_masking'] )
		repeat_masked_seq_file = input_file + ".masked"
		doc.write( cmd + '\n\n' )
		t8 = datetime.now()
		doc.write( str( t8-t7 ) + '\n\n' )
	
	
		# --- small RNA prediction => use & to just start it --- #
		output_dir = cluster_prefix + "infernal/"
		if not os.path.exists( output_dir ):
				os.makedirs( output_dir )
		cmd = run_INFERNAL_on_cluster( infernal_script, output_dir, repeat_masked_seq_file, config['small_rna_pred'] )
		infernal_output = output_dir + "FINAL_INFERNAL_RESULT.txt"
		#copy result file
		doc.write( cmd + '\n\n' )
		t9 = datetime.now()
		doc.write( str( t9-t8 ) + '\n\n' )
	
	
		# --- generation of RNA-Seq hints --- #
		output_dir = final_prefix + "rna_seq_hints/"
		if not os.path.exists( output_dir ):
				os.makedirs( output_dir )
		tmp_cluster_dir = cluster_prefix + "rna_seq_hints/"
		if not os.path.exists( tmp_cluster_dir ):
				os.makedirs( tmp_cluster_dir )
		cmd = run_STAR( STAR_script, tmp_cluster_dir, output_dir, fastq_dir, ref_sorted_contig_file, config['STAR'] )
		rna_seq_mapping_file = tmp_cluster_dir + "Aligned.sortedByCoord.out.bam"
		#grep log file
		doc.write( cmd + '\n\n' )
		t10 = datetime.now()
		doc.write( str( t10-t9 ) + '\n\n' )
	
		# --- AUGUSTUS hin-based gene prediction --- #
		output_dir = cluster_prefix + "augustus_gene_prediction/"
		if not os.path.exists( output_dir ):
				os.makedirs( output_dir )
		assembly_file = output_dir + repeat_masked_seq_file.split('/')[-1]
		os.popen( "cp " + ref_sorted_contig_file + " " + input_file )
		start_AUGUSTUS_gene_prediction( AUGUSTUS_script, assembly_file, rna_seq_mapping_file, bv_prot_file, output_dir, config['gene_prediction'] )
		#gene predictions runs for a long time on cluster => manual transfer of results required!!
		doc.write( cmd + '\n\n' )
		t11 = datetime.now()
		doc.write( str( t11-t10 ) + '\n\n' )
	
	
		doc.write( 'Pipeline completed. Do not forget to save intermediate files and reports.\n' )



if __name__ == '__main__':
	
	if '--result_dir' in sys.argv and '--cluster_dir' in sys.argv and '--assembly' in sys.argv and '--read1' in sys.argv and '--read2' in sys.argv and '--fastq_dir' in sys.argv and '--config' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

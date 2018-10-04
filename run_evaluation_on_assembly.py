### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python run_evaluation_on_assembly.py\n
					--cluster_dir <FULL_PATH_TO_OUTPUT_DIR>
					--in <INPUT_ASSEMBLY_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


import sys, os, re, glob
from operator import itemgetter

# --- end of imports --- #

def run_augustus( augustus_path, input_file, output_file ):
	"""! @brief run AUGUSTUS3.2.2 gene prediction of given input file """
	
	cmd = [ augustus_path, "--species=arabidopsis", "--gff3=on", "--UTR=on", "--uniqueGeneId=true", "--codingseq=on", input_file, ">", output_file ]
	os.popen( " ".join( cmd ) )


def extract_predicted_sequences( gff3_file, assembly_file, augustus_seqs_ex_script, working_dir ):
	"""! @brief extract predicted sequences """
	
	os.chdir( working_dir )
	cmd = [ augustus_seqs_ex_script, "--seqfile=" + assembly_file, gff3_file ]
	os.popen( " ".join( cmd ) )


def load_results_from_BLAST_result_file( BLAST_result_file, cutoff=0.9999 ):
	"""! @brief load data from BLAST result file """
	
	data = {}
	
	with open( BLAST_result_file, "r" ) as f:
		line = f.readline()
		prev_query = line.split('\t')[0]
		hits = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_query:
				sorted_hits = sorted( hits, key=itemgetter( 'score' ) )
				if len( sorted_hits ) > 1:
					if ( sorted_hits[-2]['score'] / sorted_hits[-1]['score'] ) < cutoff:
						data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
				else:
					data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
				hits = []
				prev_query = parts[0]
			hits.append( { 'query': parts[0], 'subject': parts[1], 'score': float( parts[-1] ) } )
			line = f.readline()
		sorted_hits = sorted( hits, key=itemgetter( 'score' ) )
		if len( sorted_hits ) > 1:
			if ( sorted_hits[-2]['score'] / sorted_hits[-1]['score'] ) > cutoff:
				data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
		else:
			data.update( { sorted_hits[-1]['query']: sorted_hits[-1]['subject'] } )
	return data


def compare_datasets( data1, data2, outputfile ):
	"""! @brief compares datasets and identifies bidirectional best hits """
	
	seq_IDs_of_interest = []
	
	counter = 0
	keys = data1.keys()
	with open( outputfile, "w" ) as out:
		out.write( "seq_IDs_of_file1\tseq_IDs_of_file2\n" )
		for key in keys:	#key=candidate gene
			try:
				value = data1[ key ]	#value=contig_ID
				try:
					other_value = data2[ value ]	#other_value=candidate_gene_ID
					if key == other_value:
						counter += 1
						out.write( key + '\t' + value + '\n' )
						seq_IDs_of_interest.append( value )
				except:
					pass
			except:
				pass
	return seq_IDs_of_interest


def load_multiple_fasta_file( fasta_file ):
	"""!@brief load content of multiple fasta file """
	
	content = {}
	with open( fasta_file, "r" ) as f:
		header = f.readline().strip()[1:].split(' ')[0]
		line = f.readline()
		seq = ""
		while line:
			if line[0] == '>':
				content.update( { header: seq } )
				header = line.strip()[1:].split(' ')[0]
				seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		content.update( { header: seq } )
	return content


def RBH_identification( prefix, seq_file1, seq_file2, RBH_file ):
	"""! @brief identifies RBHs between given data sets """
	
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	
	seq_file1_db = prefix + "seq_file1_db"
	seq_file2_db= prefix + "seq_file2_db"
	
	seq_file1_blast_result_file = prefix + "seq_file1_blast_result_file.txt"
	seq_file2_blast_result_file = prefix + "seq_file2_blast_result_file.txt"
	
	# --- identify RBHs --- #
	db_command1 = "makeblastdb -in " + seq_file1 + " -out " + seq_file1_db + " -dbtype 'prot' -parse_seqids"
	db_command2 = "makeblastdb -in " + seq_file2 + " -out " + seq_file2_db + " -dbtype 'prot' -parse_seqids"
	
	os.popen( db_command1 )
	os.popen( db_command2 )
	
	BLAST_command1 = "blastp -query " + seq_file1 + " -db " + seq_file2_db + " -out " + seq_file1_blast_result_file + " -outfmt 6 -evalue 0.0001 -num_threads 8"
	os.popen( BLAST_command1 )
	
	BLAST_command2 = "blastp -query " + seq_file2 + " -db " + seq_file1_db + " -out " + seq_file2_blast_result_file + " -outfmt 6 -evalue 0.0001 -num_threads 8"
	os.popen( BLAST_command2 )
	
	data1 = load_results_from_BLAST_result_file( seq_file1_blast_result_file )
	data2 = load_results_from_BLAST_result_file( seq_file2_blast_result_file )
	seq_IDs_of_interest = compare_datasets( data1, data2, RBH_file )


def run_BUSCO( input_file, prefix, busco_path, augustus_path, augustus_config_path ):
	"""! @brief run BUSCO in genome mode on given assembly file """
	
	os.chdir( prefix )
	os.environ["PATH"] = augustus_path + ":" + os.environ["PATH"]
	print os.environ["PATH"]
	os.environ["AUGUSTUS_CONFIG_PATH"] = augustus_config_path
	print os.environ["AUGUSTUS_CONFIG_PATH"]
	cmd = "python " + busco_path + " --in " + input_file + " --out busco_run > " + prefix +"log.txt"
	os.popen( cmd )


def count_genes_in_gff3_file( gff3_file ):
	"""! @brief count number of genes in given gff3 file """
	
	gene_lengths = []
	
	with open( gff3_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					gene_lengths.append( int( parts[4] ) - int( parts[3] ) )
			line = f.readline()
	return gene_lengths


def count_RBHs( RBH_file ):
	"""! @brief load all RBHs and return number """
	
	with open( RBH_file, "r" ) as f:
		f.readline()	#header
		counter = 0
		line = f.readline()
		while line:
			counter += 1
			line = f.readline()
	return counter


def calculate_formal_contig_stats( filename ):
	"""! @brief calculates some formal stats of the given multiple fasta file (assembly)
		@param filename (string) full path to a assembly output file (multiple fasta file)
		@return (dictionary) contains all formal stats of the analyzed assembly
	"""
	
	print "calculation of formal assembly stats ... please wait!"
	number_of_bases_without_N = 0	#counts all bases without N
	number_of_gc = 0		#counts occurences of G or C in sequence
	contig_lengths = []		#lengths of all contigs in the assembly; used for calculation of min, max and mean
	
	with open( filename, 'r' ) as f:
		first_line = f.readline()
		line = f.readline()
		sequence = ""
		counter = 1
		while line:
			if line[0] == '>':	#new header => evaluate current sequence and set back to empty string
				for base in sequence.upper():
					if base == 'G' or base == 'C':
						number_of_gc += 1
						number_of_bases_without_N += 1
					elif base == 'A' or base == 'T':
						number_of_bases_without_N += 1
				contig_lengths.append( len( sequence ) )
				sequence = ""
			else:
				sequence += line.strip()
			line = f.readline()
			counter += 1
			if counter % 1000 == 0:
				print str( counter/1000 ) + ' x1000 lines processed'
		#place block from new header here again (for last sequence in file)
		for base in sequence.upper():
			if base == 'G' or base == 'C':
				number_of_gc += 1
				number_of_bases_without_N += 1
			elif base == 'A' or base == 'T':
				number_of_bases_without_N += 1
		contig_lengths.append( len( sequence ) )
	
	# --- calculate remaining stats --- #
	number_of_contigs = len( contig_lengths )	#counts number of contigs / scaffolds in this assembly
	total_number_of_bases = sum( contig_lengths )	#counts all bases in the assembyl
	mean_contig_length = total_number_of_bases / number_of_contigs	#average contig lengths
	minimal_contig_length = min( contig_lengths )
	maximal_contig_length = max( contig_lengths )
	

	# --- sort list of contig length decreasing --- #
	sorted_contig_lengths = sorted( contig_lengths )[::-1]	#invert to get it decreasing
	N25 = False
	N50 = False
	N75 = False
	N90 = False
	
	cum_length = total_number_of_bases
	
	for contig_length in sorted_contig_lengths:
		cum_length -= contig_length
		if cum_length <= 0.1 * total_number_of_bases:
			if not N90:
				N90 = contig_length
		elif cum_length <= 0.25 * total_number_of_bases:
			if not N75:
				N75 = contig_length
		elif cum_length <= 0.5 * total_number_of_bases:
			if not N50:
				N50 = contig_length
		elif cum_length <= 0.75 * total_number_of_bases:
			if not N25:
				N25 = contig_length
	
	stats = { 	'number_of_contigs': number_of_contigs,
						'mean_contig_length': mean_contig_length,
						'minimal_contig_length': minimal_contig_length,
						'maximal_contig_length': maximal_contig_length,
						'total_number_of_bases': total_number_of_bases,
						'number_of_bases_without_N': number_of_bases_without_N,
						'gc_content': float( number_of_gc ) /number_of_bases_without_N,
						'N25': N25,
						'N50': N50,
						'N75': N75,
						'N90': N90
					 }	
	return stats


def main( parameters ):
	"""! @brief run all genome analysis methods """
	
	working_dir = parameters[ parameters.index('--cluster_dir')+1 ]
	input_file = parameters[ parameters.index('--in')+1 ]
	active =True
	if '--inactive' in parameters:
		active = False
	
	beetset2 = "representative_sequences.pep.fasta"
	augustus = "augustus-3.2"
	augustus_seqs_ex_script = "getAnnoFasta.pl"
	
	busco = "BUSCO/scripts/run_BUSCO.py"
	augustus_path = "augustus-3.2.2/bin/"
	augustus_config_path = "augustus-3.2.2/config/"
	
	
	if working_dir[-1] != '/':
		working_dir += "/"
	if not os.path.exists( working_dir ):
		os.makedirs( working_dir )
	
	info_file = input_file.replace( '.fasta', '.info' )
	with open( info_file, "w" ) as info_out:
		
		# --- transfer file to cluster --- #
		assembly_file = working_dir + input_file.split('/')[-1]
		if active:
			os.popen( "cp " + input_file + " " + assembly_file )
		
		#run AUGUSTUS
		gff3_file = assembly_file.replace( ".fasta", ".gff3" )
		if active:
			run_augustus( augustus, assembly_file, gff3_file )
		
		# --- extract sequences --- #
		if active:
			extract_predicted_sequences( gff3_file, assembly_file, augustus_seqs_ex_script, working_dir )
		pep_file = assembly_file.replace( ".fasta", "3.aa" )
		
		# --- identify RBHs vs. BeetSet2 --- #
		rbh_dir = working_dir + "rbh_identification/"
		RBH_file = working_dir + "RBHs_vs_BeetSet2.txt"
		if active:
			RBH_identification( rbh_dir, pep_file, beetset2, RBH_file )
		
		# --- run BUSCO v3 --- #
		prefix = working_dir + "busco/"
		if not os.path.exists( prefix ):
			os.makedirs( prefix )
		if active:
			run_BUSCO( input_file, prefix, busco, augustus_path, augustus_config_path )
		
		# --- collect all generated information in one file --- #
		stats = calculate_formal_contig_stats( input_file )
		info_out.write( "number of contigs\t" + str( stats["number_of_contigs"] ) + '\n' )
		info_out.write( "maximal contig length\t" + str( stats["maximal_contig_length"] ) + '\n' )
		info_out.write( "N50\t" + str( stats["N50"] ) + '\n' )
		info_out.write( "N90\t" + str( stats["N90"] ) + '\n' )
		info_out.write( "assembly size (withoutN)\t" + str( stats[ 'number_of_bases_without_N' ] ) + '\n' )
		info_out.write( "GC content\t" + str( stats[ 'gc_content' ] ) + '\n' )			
		
		gene_lengths = count_genes_in_gff3_file( gff3_file )
		info_out.write( "number of predicted genes\t" + str( len( gene_lengths ) ) + '\n' )
		info_out.write( "average gene lengths\t" + str( sum( gene_lengths ) / float( len( gene_lengths ) ) ) + '\n' )
		
		RBH_file = working_dir+ "RBHs_vs_BeetSet2.txt"
		number_of_RBHs = count_RBHs( RBH_file )
		info_out.write( "number of RBHs against BeetSet2\t" + str( number_of_RBHs ) + '\n' )
		
		busco_summary_file = glob.glob( working_dir + "busco/run_busco_run/short_summary_busco_run.txt" )[0]	#this error should be fixed now
		with open( busco_summary_file, "r" ) as f:
			content = f.read()
			BUSCO_string = re.findall( "C:\d+\.\d+%\[S:\d+\.\d+%,D:\d+\.\d+%\],F:\d+\.\d+%,M:\d+\.\d+%,n:\d+", content )[0]
		info_out.write( "BUSCO result string\t" + BUSCO_string + "\n" )
		

if __name__ == '__main__':
	
	if '--cluster_dir' in sys.argv and '--in' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

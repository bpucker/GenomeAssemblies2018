### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

import sys, glob, re, os, time, datetime, shutil

# --- end of imports --- #

__usage__ = """
		python run_blastn_on_cluster.py\n
		--assembly <FULL_PATH_TO_INPUT_FILE (FASTA)>\n
		--output_dir <FULL_PATH_TO_OUTPUT_DIR>\n
				
		feature requests and bug reports:
		bpucker@cebitec.uni-bielefeld.de
"""


def submit_jobs_to_cluster( prefix, query_file_names, para_jobs, infernal_path, infernal_db ):
	"""! @brief submit BLAST jobs for each file to cluster """
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, file_name in enumerate( query_file_names ):
		ID = "R_" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		out_file = prefix + ID + '.out'
		err_file = prefix + ID + '.err'
		result_file = prefix + ID + ".infernal"
		
		cmd = infernal_path + " --cpu 0 " + infernal_db + " " + file_name + " > " + result_file
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-N",
																ID,
																"-l vf=1G",
																"-l arch=lx-amd64",
																"-P fair_share",
																"-o",
																out_file,
																"-e",
																err_file
															] ) + '\n'
							)
		os.popen( "chmod +x " + sh_file )
		os.popen( sh_file )
		time.sleep(1)
		#os.remove( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "R_" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 10 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "R_" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 120 )


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split('.')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split('.')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def produce_multiple_query_files( query_file, cutoff ):
	"""! @brief produce multiple query files """
	
	prefix = query_file + str( datetime.datetime.now() )[-5:]  + "/"
	os.makedirs( prefix )
	
	sequences = load_sequences( query_file )
	
	query_file_names = []
	
	len_counter = 0
	name_counter = 1
	query_file = prefix + "0".zfill(4) + ".fasta"
	query_file_names.append( query_file )
	out = open( query_file, "w" )
	for idx, seq_id in enumerate( sorted( sequences.keys() ) ):
		if len_counter >= cutoff:
			len_counter = 0
			out.close()
			query_file = prefix + str( name_counter ).zfill(4) + ".fasta"
			query_file_names.append( query_file )
			out = open( query_file, "w" )
			name_counter += 1
		out.write( '>' + seq_id + '\n' + sequences[ seq_id ] + '\n' )
		len_counter += len( sequences[ seq_id ] )
	out.close()
	return prefix, query_file_names


def load_contig_order( query_file ):
	"""! @brief load order of contigs in query file """
	
	contig_order = []
	with open( query_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] == '>':
				contig_order.append( line.strip()[1:].split(' ')[0] )
			line = f.readline()
	return contig_order


def load_INFERNAL_results( result_file ):
	"""! @brief load all results from given INTERNAL result file """
	
	results_per_contig = {}
	with open( result_file, "r" ) as f:
		line = f.readline()
		ID = ""
		while line:
			try:
				if line[:len( "Query:" )] == "Query:":
					ID = re.findall( "scaffold\d+", line )[0]
				elif len( ID ) > 0 and line[0] != '#' and line.strip()[0] == "(":	#try to match table lines
					if "!" in line and not ".." in line:	#some additional criteria
						mod_line = line.strip().replace( "    ", " " ).replace( "   ", " " ).replace( "  ", " " ).replace( "  ", " " ).replace( "  ", " " )
						parts = mod_line.split(' ')
						if len( parts[6] ) > 0 and len( parts[7] ) > 0:
							start, end = sorted( map( int, parts[6:8] ) )
							entry = { 'chr': ID,
											'start':  start,
											'end': end,
											'comment': parts[5],
											'score': parts[3],
											'orientation': parts[8] 
										}
							try:
								content = results_per_contig[ ID ]
								del results_per_contig[ ID ]
								content.append( entry )
								results_per_contig.update( { ID: content } )
							except KeyError:
								results_per_contig.update( { ID: [ entry ] } )
			except IndexError:
				pass
			line = f.readline()
	return results_per_contig


def final_processing( prefix, query_file ):
	"""! @brief blt processing of BLAT results for identification of best hit """
	
	final_result_file = "/".join( query_file.split('/')[:-1] ) + "FINAL_INFERNAL_RESULT.txt"
	
	result_file_names = glob.glob( prefix + "*.infernal" ) + glob.glob( prefix + "*/*.infernal" ) 
	
	contig_order = load_contig_order( query_file )	#get order of initial sequences
	
	results_per_sequence = {}
	for result_file in result_file_names:
		results = load_INFERNAL_results( result_file )	#load all results from result files
		results_per_sequence.update( results )
	
	# --- write ordered results into final output file (GFF3) --- #
	with open( final_result_file, "w" ) as out:
		for contig in contig_order:
			try:
				results = results_per_sequence[ contig ]
			except KeyError:
				results = []
			for each in results:
				new_line = [ each['chr'], "INFERNAL", "ncRNA", each['start'], each['end'], each['score'], each['orientation'], ".", each['comment'] ]
				out.write( "\t".join( map( str, new_line ) ) + '\n' )
	
	#shutil.rmtree( prefix )	#remove all temp files


def main( arguments ):
	"""! @brief check inputs and call functions """
	
	if '--para_jobs' in arguments:
		para_jobs = int( arguments[ arguments.index( '--para_jobs' ) + 1 ] )
	else:
		para_jobs = 50
	
	query_file = arguments[ arguments.index( '--assembly' ) + 1 ]
	
	if '--splitting_cutoff' in arguments:
		cutoff = int( arguments[ arguments.index( '--splitting_cutoff' ) + 1 ] )
	else:
		cutoff=500000
	
	infernal_path = "infernal-1.1.2-linux-intel-gcc/binaries/cmscan"
	infernal_db = "infernal-1.1.2-linux-intel-gcc/Rfam/Rfam.cm"
	
	
	output_dir = arguments[ arguments.index( '--output_dir' ) + 1 ]
	if output_dir[-1] != '/':
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	os.popen( "cp " + query_file + " " + output_dir )
	query_file = output_dir + query_file.split('/')[-1]
	
	prefix, query_file_names = produce_multiple_query_files( query_file, cutoff )	#seq length increads compared to BLAT
	
	submit_jobs_to_cluster( prefix, query_file_names, para_jobs, infernal_path, infernal_db )
	
	final_processing( prefix, query_file )


if __name__ == '__main__':
	
	if '--assembly' in sys.argv and '--output_dir' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

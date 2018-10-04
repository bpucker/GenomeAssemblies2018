### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python adapter_clipping.py\n
					--assembly_file <FULL_PATH_TO_FASTA_FILE>
					--output_dir <FULL_PATH_TO_OUTPUT_DIRECTORY>
					
					--adapter_file <FULL_PATH_TO_FASTA_FILE> (optional if hard coded)
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import os, sys

# --- end of imports --- #


def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( multiple_fasta_file ) as f:
		header = f.readline()[1:].strip().split(' ')[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split(' ')[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )
	return sequences


def adapter_trimming( input_file, assembly ):
	"""! @brief collect information for clipping of illumina adapters """
	
	# --- load information for adapter clippping --- #
	black_list = {}
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				current_clipping = black_list[ parts[1] ]
				current_clipping.append( sorted( [ int( parts[8] ), int( parts[9] ) ] ) )
				del black_list[ parts[1] ]
				black_list.update( { parts[1]: current_clipping } )
			except KeyError:
				black_list.update( { parts[1]: [ sorted( [ int( parts[8] ), int( parts[9] ) ] ) ] } )
			line = f.readline()
	
	# --- do adapter clipping --- #
	lost_sequences_counter = 0
	clipped_sequences_counter = 0
	for key in black_list.keys():
		seq = assembly[ key ]
		positions = black_list[ key ]
		del assembly[ key ]
		start, end = positions[0]
		for each in positions:
			if each[0] < start:
				start = 0 + each[0]
			if each[1] > end:
				end = 0 + each[1]
		
		if start > 200 and len( seq ) - end > 200:	#take longest remaining seq if adapter is located in the middle
			if start > len( seq ) - end:
				assembly.update( { key: seq[ :start ] } )
			else:
				assembly.update( { key: seq[ end: ] } )
			clipped_sequences_counter += 1
		elif start > 200:
			assembly.update( { key: seq[ :start ] } )
			clipped_sequences_counter += 1
		elif len( seq ) - end > 200:
			assembly.update( { key: seq[ end: ] } )
			clipped_sequences_counter += 1
		else:
			print "sequence lost due to length cutoff: " + key
			lost_sequences_counter += 1
	
	print "number of clipped sequences: " + str( clipped_sequences_counter )
	print "number of lost sequences (due to length cutoff): " + str( lost_sequences_counter )
	return assembly


def main( arguments ):
	"""! @brief run everything """
	
	assembly_file = arguments[ arguments.index( '--assembly_file' )+1 ]
	output_dir = arguments[ arguments.index( '--output_dir' )+1 ]
	
	if output_dir[-1] != '/':
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	if '--adapter_file' in arguments:
		adapter_file = arguments[ arguments.index( '--adapter_file' )+1 ]
	else:
		adapter_file = "Illumina_adapters.fa"	### ADD YOUR FILE HERE ###
		
	blast_result_file = output_dir + "illumina_adapters_vs_assembly.txt"
	blast_db = output_dir + "blast_db"
	clean_assembly_file = output_dir + "clean_assembly.fasta"
	
	print "running BLASTn ... "
	os.popen( "makeblastdb -in " + assembly_file + " -out " + blast_db + " -dbtype nucl" )
	os.popen( "blastn -query " + adapter_file + " -db " + blast_db + " -out " + blast_result_file + " -outfmt 6 -evalue 0.01 -word_size 4 -num_threads 4" )
	
	assembly = load_sequences( assembly_file )
	clean_assembly = adapter_trimming( blast_result_file, assembly )
	
	with open( clean_assembly_file, "w" ) as out:
		for key in sorted( clean_assembly.keys() ):
			out.write( '>' + key + '\n' + clean_assembly[ key ] + '\n' )
	print "all done!"


if __name__ == '__main__':
	if '--assembly_file' in sys.argv and '--output_dir' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

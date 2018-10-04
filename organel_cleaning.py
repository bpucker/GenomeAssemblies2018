import os, sys

# --- end of imports --- #

def load_all_seqs_from_multiple_fasta_file( filename ):
	"""! @brief load all sequences from multiple fasta file """
	
	data = {}
	
	with open( filename, "r" ) as f:
	 	header = f.readline().strip()[1:]
		line = f.readline()
		seq = []
		while line:
			if line[0] == '>':
				data.update( { header: "".join( seq ) } )
				header = line.strip()[1:]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		data.update( { header: "".join( seq ) } )
	return data


def load_blast_results( result_file, similarity, length ):
	"""! @brief load all BLAST results for further filtering """
	
	hits = []
	with open( result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[2] ) > similarity:
				if int( parts[3] ) > length:
					hits.append( { 'id': parts[0], 'len': abs( int( parts[6] )-int( parts[7] ) ) } )
			line = f.readline()
	return hits


def main( arguments ):
	"""! @brief controls everything """
	
	
	prefix = arguments[ arguments.index( '--prefix' )+1 ]
	input_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	
	bv_organell_seq_file = "organel_sequences.fasta"	### REPLACE THIS FILE NAME ###
	similarity = 85
	length = 500
	length_fraction = 0.8
	
	blast_db = prefix + "blastdb"
	os.popen( "makeblastdb -in " + bv_organell_seq_file + " -out " + blast_db + " -dbtype nucl" )
	
	blast_results_file = prefix + "blast_results.txt"
	os.popen( "blastn -query " + input_file + " -db " + blast_db + " -out " + blast_results_file + " -outfmt 6 -evalue 0.0000000001 -num_threads 8" )
	
	blast_results = load_blast_results( blast_results_file, similarity, length )
	seqs = load_all_seqs_from_multiple_fasta_file( input_file )
	
	for hit in blast_results:
		try:
			if hit['len'] > length_fraction*len( seqs[ hit['id'] ] ):
				del seqs[ hit['id'] ]
		except KeyError:
			pass
	with open( output_file, "w" ) as out:
		for key in seqs.keys():
			out.write( '>' + key + '\n' + seqs[ key ] + '\n' )


if __name__ == '__main__':
	
	main( sys.argv )

### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python sort_contigs_on_ref.py\n
					--contig_file <FULL_PATH_TO_FILE>
					--ref_file <FULL_PATH_TO_FILE>
					--output_dir <FULL_PATH_TO_DIR>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
"""


import os, sys
from operator import itemgetter

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip().split( " " )[0]
		seq = ""
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: seq } )
					header = line.strip()[1:].split( " " )[0]
					seq = ""
			else:
				seq += line.strip()
			line = f.readline()
		sequences.update( { header: seq } )	
	return sequences


def revcomp( seq ):
	"""! @brief construct reverse complement of sequence """
	
	new_seq = []
	
	bases = { 'a':'t', 't':'a', 'c':'g', 'g':'c' }
	for nt in seq.lower():
		try:
			new_seq.append( bases[nt] )
		except:
			new_seq.append( 'n' )
	return ''.join( new_seq[::-1] ).upper()


def blastn_and_sorting( contigs_file, ref_seq_file, prefix, contig_order_file, contigs ):
	""""! @brief run BLASTn and identify optimal order and orientation of contigs """
	
	cmd = "makeblastdb -in " + ref_seq_file + " -out " + prefix+"ref_blastn_db -dbtype 'nucl'"
	os.popen( cmd )
	
	blastn_output_file = prefix + "blastn_output_file.txt"
	cmd = "blastn -query " + contigs_file + " -db " + prefix+"ref_blastn_db -out " + blastn_output_file + " -outfmt 6 -evalue 0.00000001 -num_threads 8" 
	os.popen( cmd )
	
	best_blastn_hits = {}
	
	with open( blastn_output_file, "r" ) as f:
		line = f.readline()
		prev_id = ""
		while line:
			parts = line.strip().split('\t')
			if parts[0] != prev_id:
				if int( parts[8] ) < int( parts[9] ):
					orientation = "fw"
				else:
					orientation = "rv"
				try:
					value = best_blastn_hits[ parts[0] ]['score']
					if value < float( parts[-1] ):
						del best_blastn_hits[ parts[0] ]
						best_blastn_hits.update( { parts[0]: { 'id': parts[0], 'chr': parts[1], 'pos': (int( parts[8] )+int( parts[9] ))/2, 'orientation': orientation, 'score': float( parts[-1] ) } } )
				except KeyError:
					best_blastn_hits.update( { parts[0]: { 'id': parts[0], 'chr': parts[1], 'pos': (int( parts[8] )+int( parts[9] ))/2, 'orientation': orientation, 'score': float( parts[-1] ) } } )
			line = f.readline()
	
	best_blastn_hits = best_blastn_hits.values()	
	sorted_hits = sorted( best_blastn_hits, key=itemgetter( 'chr', 'pos' ) )
	
	with open( contig_order_file, "w" ) as out:
		for hit in sorted_hits:
			try:
				out.write( '\t'.join( map( str, [ hit['id'], hit['chr'], hit['pos'], hit['orientation'], len( contigs[ hit['id'] ] ) ] ) ) + '\n' )
			except:
				print hit['id']


def produce_final_contig_order( contig_order_file, assembly_file, final_contig_file ):
	"""! @brief load contig order information from file """
	
	contigs = load_sequences( assembly_file )
	
	counter = 1
	with open( final_contig_file, "w" ) as out:
		with open( contig_order_file, "r" ) as f:
			line = f.readline()
			prev_chr_name = line.split('\t')[1]
			
			while line:
				parts = line.strip().split('\t')
				if parts[ 3 ] == "fw":
					out.write( '>scaffold' + str( counter ).zfill(6) + '\n' + contigs[ parts[0] ] + '\n' )
					print parts[0] + ">>" + 'scaffold' + str( counter ).zfill(6)
					del contigs[ parts[0] ]
					
				elif parts[3] == "rv":
					out.write( '>scaffold' + str( counter ).zfill(6) + '\n' + revcomp( contigs[ parts[0] ] ) + '\n' )
					print parts[0] + ">>" + 'scaffold' + str( counter ).zfill(6)
					del contigs[ parts[0] ]
					
				else:
					print "ERROR: SEQ ORIENTATION UNKNOWN!"
				line = f.readline()
				counter += 1
		
		# --- writing all remaining contigs sorted by length into output file --- #
		sorting_list = []
		for key in contigs:
			sorting_list.append( { 'id': key, 'len': len( contigs[ key ] ) } )
		for entry in sorted( sorting_list, key=itemgetter('len') )[::-1]:
			out.write( '>scaffold' + str( counter ).zfill(6) + '\n' + contigs[ entry['id'] ] + '\n' )
			print entry['id']+ ">>" + 'scaffold' + str( counter ).zfill(6)
			counter += 1


def main( parameters ):
	"""! @brief runs all parts of reference-based pseudochromosome construction """
	
	contig_file = parameters[ parameters.index( '--contig_file' )+1 ]
	ref_seq_file = parameters[ parameters.index( '--ref_file' )+1 ]
	
	prefix = parameters[ parameters.index( '--output_dir' )+1 ]
	
	contig_order_file = prefix + "contig_order.txt"
	final_contig_file = prefix + "sorted_contigs.fasta"
	
	contigs = load_sequences( contig_file )
	blastn_and_sorting( contig_file, ref_seq_file, prefix, contig_order_file, contigs )	
	produce_final_contig_order( contig_order_file, contig_file, final_contig_file )


if __name__ == '__main__':
	
	if '--contig_file' in sys.argv and '--ref_file' in sys.argv and '--output_dir' in sys.argv:
		main( sys.argv )	
	else:
		sys.exit( __usage__ )

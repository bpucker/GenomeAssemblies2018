### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python purging_short_sequences.py\n
					--in <ASSEMBLY_FILE>
					--out <CLEANED_ASSEMBLY_FILE>
					"""

import os, sys
import matplotlib.pyplot as plt

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


def main( parameters ):
	"""! @brief run everything """
	
	input_assembly_file = parameters[ parameters.index( '--in' ) +1 ]
	output_assembly_file = parameters[ parameters.index( '--out' ) +1 ]
	
	contig_min_cutoff = 500
	max_N_percentage = 0.5
	
	contigs = load_sequences( input_assembly_file )
	
	
	long_counter = 0
	blast_hit_counter = 0
	dropped_counter = 0
	very_long_counter = 0
	
	with open( output_assembly_file, "w" ) as out:
		for key in sorted( contigs.keys() ):
			if len( contigs[ key ] ) >= contig_min_cutoff and contigs[ key ].count('N') < max_N_percentage*len( contigs[ key ] ):
				out.write( '>' + key + '\n' + contigs[ key ] + '\n' )
				long_counter += 1
			else:
				dropped_counter += 1
	
	print "sequence kept due to length: " + str( long_counter )
	print "sequence kept due to annotation: " + str( blast_hit_counter )
	print "dropped sequence: " + str( dropped_counter )
	
	fig_file = output_assembly_file + ".png"
	fig, ax = plt.subplots()
	
	N_content = []
	for key in contigs.keys():
		N_content.append( contigs[ key ].count('N') / float( len( contigs[ key ] ) ) )
	
	ax.hist( N_content, bins=100 )
	ax.set_xlabel( "% of N in sequence" )
	ax.set_ylabel( "number of sequences" )
	
	fig.savefig( fig_file, dpi=300 )


if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

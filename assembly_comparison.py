### Boas Pucker ###

import re, sys, glob
import matplotlib.pyplot as plt
from operator import itemgetter
import numpy as np

# --- end of imports --- #

def construct_feature_response_curve( assembly_data, output_figure, estimated_genome_size ):
	"""! @brief construct figure for visualization of assembly stats """
	
	fig, ax1 = plt.subplots( figsize=( 12, 5 ) )
	ax2 = ax1.twinx()
	ax3 = ax1.twinx()
	ax4 = ax1.twinx()
	ax5 = ax1.twinx()
	ax6 = ax1.twinx()
	
	ax1_values = []
	ax2_values = []
	ax3_values = []
	ax4_values = []
	ax5_values = []
	ax6_values = []
	labels = []
	
	
	for i, assembly in enumerate( assembly_data ):
		labels.append( assembly['id'] )
		try:
			ax6_values.append( assembly['gene_len'] )
		except KeyError:
			ax6_values.append( 0 )
		try:
			ax5_values.append( assembly['RBHs'] )
		except KeyError:
			ax5_values.append( 0 )
		try:
			ax4_values.append( assembly['BUSCOs'] )
		except KeyError:
			ax4_values.append( 0 )
		try:
			ax3_values.append( assembly['N50'] )
		except KeyError:
			ax3_values.append( 0 )
		try:
			ax2_values.append( assembly['number_of_contigs']/1000 )
		except KeyError:
			ax2_values.append( 0 )
		try:
			ax1_values.append(assembly['size']/1000000.0  )
		except KeyError:
			ax1_values.append( 0 )
	
	# --- plot feature response curve --- #
	ax6.scatter( range( len( ax6_values ) ), ax6_values, color="grey", s=5 )
	ax5.scatter( range( len( ax5_values ) ), ax5_values, color="green", s=5 )
	ax4.scatter( range( len( ax4_values ) ), ax4_values, color="purple", s=5 )
	ax3.scatter( range( len( ax3_values ) ), ax3_values, color="orange", s=5 )
	ax2.scatter( range( len( ax2_values ) ), ax2_values, color="blue", s=5 )
	ax1.scatter( range( len( ax1_values ) ), ax1_values, color="red", s=5 )
	
	# --- add expected genome size --- #
	ax1.plot( [ 0, len( ax4_values ) ], [ estimated_genome_size, estimated_genome_size ], color="red", linestyle="--", linewidth=1 )
	ax1.text( 0, estimated_genome_size, "estimated genome size", fontsize=5, color="red", va="top" )
	
	#labels.append( "" )
	ax1.set_ylabel( "assembly size [Mbp]", color="red" )
	ax2.set_ylabel( "number of contigs [*1000]", color="blue" )
	ax3.set_ylabel( "N50 [bp]", color="orange" )
	ax4.set_ylabel( "complete BUSCOs [%]", color="purple" )
	ax5.set_ylabel( "RBHs vs. BeetSet2", color="green" )
	ax6.set_ylabel( "average gene length", color="grey" )
	
	ax1.set_xlabel( "assemblies with different k-mer sizes" )
	
	rspine = ax3.spines['right']
	rspine.set_position(('axes', 1.12 ) )
	
	lspine = ax4.spines['right']
	lspine.set_position( ('axes', 1.25 ) )
	
	lspine = ax5.spines['right']
	lspine.set_position( ('axes', 1.4 ) )
	
	lspine = ax6.spines['right']
	lspine.set_position( ('axes', 1.6 ) )
	
	
	ax1.set_xlim( 0, len( labels ) )
	start, end = ax1.get_xlim()
	ax1.xaxis.set_ticks(np.arange(start, end, 1))
	ax1.set_xlim( -0.5, len( labels ) )
	
	ax1.set_xticklabels( labels, rotation=45 )
	ax1.set_ylim( int( min( ax1_values )*0.9 ), int( max( ax1_values )*1.1 ) )
	ax2.set_ylim( int( min( ax2_values )*0.9 ), int( max( ax2_values )*1.1 ) )
	ax3.set_ylim( int( min( ax3_values )*0.9 ), int( max( ax3_values )*1.1 ) )
	ax4.set_ylim( int( min( ax4_values )*0.9 ), int( max( ax4_values )*1.1 ) )
	ax5.set_ylim( int( min( ax5_values )*0.9 ), int( max( ax5_values )*1.1 ) )
	ax6.set_ylim( int( min( ax6_values )*0.9 ), int( max( ax6_values )*1.1 ) )
	
	plt.subplots_adjust( left=0.05, right=0.6, top=0.97, bottom=0.15 )
	
	fig.savefig( output_figure, dpi=1200 )
	plt.close('all')


def load_data_from_info_file( info_file ):
	"""! @brief load all data from info file"""
	
	collected_infos = { 'id': info_file.split('/')[-1].split('.')[0], 'species': info_file.split('/')[-1].split('.')[0][:3] }
	
	with open( info_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if parts[0] == "number of contigs":
				collected_infos.update( { 'number_of_contigs': int( parts[1] ) } )
			elif parts[0] == "maximal contig length":
				collected_infos.update( { 'max_contig_len': int( parts[1] ) } )
			elif parts[0] == "N50":
				collected_infos.update( { 'N50': int( parts[1] ) } )
			elif parts[0] == "N90":
				collected_infos.update( { 'N90': int( parts[1] ) } )
			elif parts[0] == "assembly size (withoutN)":
				collected_infos.update( { 'size': int( parts[1] ) } )
			elif parts[0] == "GC content":
				collected_infos.update( { 'GC': float( parts[1] ) } )
			elif parts[0] == "number of predicted genes":
				collected_infos.update( { 'gene_number': int( parts[1] ) } )
			elif parts[0] == "average gene lengths":
				collected_infos.update( { 'gene_len': float( parts[1] ) } )
			elif parts[0] == "number of RBHs against BeetSet2":
				collected_infos.update( { 'RBHs': int( parts[1] ) } )
			elif parts[0] == "BUSCO result string":
				collected_infos.update( { 'BUSCOs': float( parts[1].split(':')[1].split('%')[0] ) } )
				collected_infos.update( { 'BUSCO_info': parts[1] } )
			line = f.readline()
	return collected_infos


def main( ):
	"""! @brief run all parts of this analysis """
	
	input_directory = "/genome_assemblies_of_interest/"
	output_directory = "/genome_assembly_evaluation_summary/"
	
	stat_key_names = [ 'id', "number_of_contigs", "max_contig_len", "N50", "N90", "size", "GC", "gene_number", "gene_len", "RBHs", "BUSCOs", "BUSCO_info" ]
	
	estimated_genome_sizes = { 'kew': 620, 'lim': 730, 'mac': 500, 'phe': 265, 'sim': 765, 'bet': 550, 'cor': 280, 'spe': 245 }	#Mbp	#'mic': 500, 
	
	assembly_stats = []
	for filename in glob.glob( input_directory + "*.info" ):
		assembly_stats.append( load_data_from_info_file( filename ) )
	
	species = []
	for each in assembly_stats:
		species.append( each['species'] )
	species = list( set( species ) )
	
	print "number of different species: " + str( len( species ) )
	for spec in species:
		output_stats_file = output_directory + spec + ".txt"
		output_pic_file = output_directory + spec + ".png"
		relevant_data_sets = []
		for entry in assembly_stats:
			if entry['species'] == spec:
				relevant_data_sets.append( entry )
		relevant_data_sets = sorted( relevant_data_sets, key=itemgetter('id') )
		
		# --- wite statistics into output file --- #
		with open( output_stats_file, "w" ) as out:
			for key in stat_key_names:
				new_line = [ key ]
				for entry in relevant_data_sets:
					try:
						new_line.append( entry[ key ] )
					except:
						new_line.append( "N/A" )
				out.write( "\t".join( map( str, new_line ) ) + '\n' )
		
		# --- construct figure --- #
		estimated_genome_size = estimated_genome_sizes[ spec ]	#in Mbp
		construct_feature_response_curve( relevant_data_sets, output_pic_file, estimated_genome_size )


if __name__ == '__main__':
	
	main()
	
	print "all done!"

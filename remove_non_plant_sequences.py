### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

import sys, glob, re, os, time, datetime, shutil, urllib2
from operator import itemgetter

# --- end of imports --- #

__usage__ = """ python run_blastn_on_cluster.py\n
							--assembly_file <INPUT_FILE (FASTA)>\n
							--plant_ref_file <INPUT_FILE (FASTA)>\n
							--tmp_cluser_dir <OUTPUT_DIRECTORY>\n
							--final_result_dir <OUTPUT_DIRECTORY>\n
							--active <activates execution of BLASTs>
					
							feature requests and bug reports:
							bpucker@cebitec.uni-bielefeld.de
						"""


def submit_jobs_to_cluster( prefix, query_file_names, reference_blastn_db, para_jobs ):
	"""! @brief submit BLAST jobs for each file to cluster """
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, file_name in enumerate( query_file_names ):
		ID = "B_" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = prefix + ID + '.sh'
		out_file = prefix + ID + '.out'
		err_file = prefix + ID + '.err'
		
		cmd = "blastn -query " + file_name + " -db " + reference_blastn_db + " -out " +  '.'.join( file_name.split('.')[:-1] ) + ".txt " 
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + "-outfmt 6 -evalue 0.00001 -max_target_seqs 1" +'"',
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
		os.remove( sh_file )
		waiting_status = True
		while waiting_status:
			qstat = os.popen( "qstat" )
			content = qstat.read()
			qstat_IDs = re.findall( "B_" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )
					
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 3 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "B_" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				for each in content.split('\n')[2:-1]:
					if ID in each.split()[2] and not 'd' in each.split()[4]:
						waiting_status = True
		time.sleep( 10 )


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


def final_processing( query_file_names, prefix, query_file, final_result_file ):
	"""! @brief blt processing of BLAT results for identification of best hit """
	
	result_file_names = []
	for filename in query_file_names:
		result_file_names.append( '.fasta'.join( filename.split('.fasta')[:-1] ) + '.txt' )
	
	cmd1 = "cat " + " ".join( result_file_names ) + " > " + final_result_file
	os.popen( cmd1 )
	
	#shutil.rmtree( prefix )
	return final_result_file


def run_blastn_on_cluster( query_file, reference_blastn_db, final_result_file, para_jobs ):
	"""! @brief check inputs and call functions """
	
	cutoff=1000000
	
	prefix, query_file_names = produce_multiple_query_files( query_file, cutoff )
	
	submit_jobs_to_cluster( prefix, query_file_names, reference_blastn_db, para_jobs )
	
	final_processing( query_file_names, prefix, query_file, final_result_file )


def load_blast_results( result_file ):
	"""! @brief load all hit taxonomic IDs and score for all queries """
	
	results = {}
	
	with open( result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				hits = results[ parts[ 0 ] ]
				hits.append( { 'id': parts[1], 'score': float( parts[-1] ) } )		#NEED TO BE CHANGED TO -1
				del results[ parts[ 0 ] ]
				results.update( { parts[0]: hits } )
			except KeyError:
				results.update( { parts[0]: [ { 'id': parts[1], 'score': float( parts[-1] ) } ] } )	#NEED TO BE CHANGED TO -1
			line = f.readline()
	return results


def load_high_quality_blast_results( result_file ):
	"""! @brief load all hit taxonomic IDs and score for all queries """
	
	results = {}
	
	with open( result_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if float( parts[2] ) > 80 and int( parts[3] ) > 100:
				try:
					hits = results[ parts[ 0 ] ]
				except KeyError:
					results.update( { parts[0]: float( parts[-4] ) } )
			line = f.readline()
	return results


def load_taxonomic_data( taxonomic_table_file ):
	
	taxonomic_data = {}
	
	with open( taxonomic_table_file, "r" ) as f:
		header = f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			
			line = f.readline()


def get_information_from_NCBI( gis, gi_look_up_table_file ):
	"""! @brief get all sequences corresponding to IDs provided in a text file from NCBI """
	
	gi_look_up_table = {}
	with open( gi_look_up_table_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				gi_look_up_table.update( { parts[0]: parts[1] } )
			except:
				print line
			line = f.readline()
	
	accession_numbers = []
	for gi in gis:
		try:
			gi_look_up_table[ gi ]
		except KeyError:
			accession_numbers.append( gi )
	
	url_prefix = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id="
	with open( gi_look_up_table_file, "a", 0 ) as out:
		for idx, accession_ID in enumerate( accession_numbers ):
			#print str( idx+1 ) + "/" + str( len( accession_numbers ) )
			try:
				parts = urllib2.urlopen( url_prefix + accession_ID ).read().strip().split('\n')
				out.write( accession_ID + "\t" + parts[0][1:].lower() + '\n' )
			except:
				pass	#print "ERROR: " + accession_ID
			if idx % 5 == 0:
				time.sleep( 1 )
	return gi_look_up_table


def load_sequenced_species_names( input_file ):
	"""! @brief load all species names into one string """
	
	genomes = []
	
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			genomes.append( line.split('\t')[0] )			
			line = f.readline()
	#print "number of loaded genome names: " + str( len( genomes ) )
	return " ".join( genomes ).lower()


def extract_annotation_from_string( input_string ):
	""""! @brief extracts species name from given annotation string """
	
	if "predicted: " in input_string:
		return " ".join( input_string.split('predicted: ')[1].split(' ')[:2] )	
	else:
		return " ".join( input_string.split(' ')[1:3] )


def filter_contamination_sequences( blast_results, all_seqs,  plant_table_file, prokaryote_table_file, gi_look_up_table_file ):
	"""! @brief identify plant sequences and contaminations based on information about sequenced genomes """
	
	collected_gis = []
	for key in blast_results.keys():
		for entry in blast_results[ key ]:
			try:
				collected_gis.append( re.findall( "gi\|\d+", entry['id'] )[0] )
			except IndexError:
				print "ERROR: no gi detected - " + each			
	collected_gis = list( set( collected_gis ) )
	print "number of collected gis: " + str( len( collected_gis ) )
	gi_look_up_table = get_information_from_NCBI( collected_gis, gi_look_up_table_file )
	
	prokaryote_genomes = load_sequenced_species_names( prokaryote_table_file )
	plant_genomes = load_sequenced_species_names( plant_table_file )
	
	###
	plant_genomes += " vitis vinifera malus x	103 erythranthe guttatus vernicia fordii gossypium harknessii vernicia montana tradescantia ohiensis gossypioides kirkii hesperelaea palmeri "
	plant_genomes += "populus tremula croton texensis croton stellatopilosus olea europaea m.truncatula dna millettia pinnata camellia sinensis poncirus trifoliata heuchera parviflora croton stellatopilosus"
	plant_genomes += " populus tomentosa picea sitchensis mimulus guttatus populus est boea hygrometrica malus hupehensis geranium maderense salix purpurea pisum sativum croton sublyratus "
	plant_genomes += "croton insularis croton bonplandianus vitis hybrid solanum demissum allium cepa bambusa oldhamii plantago major salix suchowensis batis maritima croton floribundus "
	plant_genomes += " sinapis arvensis populus davidiana musa abb medicago sativa glycine tomentella vitis amurensis stevia rebaudiana populus balsamifera viscum scurruloideum soybean clone "
	plant_genomes += " salvia miltiorrhiza vitis labrusca capsicum frutescens grapevine fanleaf hirtella racemosa populus tremuloides prunus pseudocerasus humulus japonicus populus x petunia x "
	plant_genomes += " oryza granulata [camellia sinensis]camellia lupinus albus camellia saluenensis morus alba vitis berlandieri croton viminalis croton megalobotrys croton rosmarinoides croton michauxii "
	plant_genomes += " croton ekmanii "
	###
	
	critical_annotations = []
	
	plant_seq_ids = {}
	contamination_seq_ids = {}
	for key in blast_results.keys():
		plant_status = 0
		contamination_status = 0
		for entry in blast_results[ key ]:
			try:
				try:
					annotation = extract_annotation_from_string( gi_look_up_table[ re.findall( "gi\|\d+", entry['id'] )[0] ] )
					if annotation in plant_genomes:
						plant_status += entry['score']
					elif annotation in prokaryote_genomes:
						contamination_status += entry['score']
					#else:
						#try:
							#annotation_black_list[ annotation ]
							#contamination_status += 1000
						#except KeyError:
							#critical_annotations.append( annotation )
				except IndexError:
					print "ERROR: no gi number detected in - " + entry['id']
			except KeyError:
				print "ERROR: gi number not found - " + entry['id']
		if plant_status > 0 and plant_status >= contamination_status:
			plant_seq_ids.update( { key: None } )
		elif contamination_status > 0:
			contamination_seq_ids.update( { key: None } )
	
	sorted_critical_annotations = []
	for each in list( set( critical_annotations ) ):
		sorted_critical_annotations.append( { 'id': each, 'value': critical_annotations.count( each ) } )
	sorted_critical_annotations = sorted( sorted_critical_annotations, key=itemgetter('value') )[::-1]
	for each in sorted_critical_annotations:
		print each['id'] + '\t' + str( each['value'] )
	
	return plant_seq_ids, contamination_seq_ids


def main( arguments ):
	"""! @brief controls all functions of this workflow """
	
	assembly_file = arguments[ arguments.index( '--assembly_file' )+1 ]
	prefix = arguments[ arguments.index( '--tmp_cluster_dir' )+1 ]
	
	
	plant_ref_file = arguments[ arguments.index( '--plant_ref_file' )+1 ]
	output_dir = arguments[ arguments.index( '--final_result_dir' )+1 ]
	
	
	black_list = ""
	
	para_jobs = 200
	
	if '--active' in arguments:
		active = True
	else:
		active = False
	
	if output_dir[-1] != "/":
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	# --- constant file locations --- #
	plant_table_file = "genomes_euks.txt"	#needs to be a list of all sequenced eukaryotic species
	prokaryote_table_file = "genomes_proks.txt"	#needs to be a list of all sequences prokaryotic species
	gi_look_up_table_file = "gi_look_up_table.txt"	#should be changed to a writeable directory
	
	
	# --- output files --- #
	cleaned_assembly_file = output_dir + "cleaned_assembly_file.fasta"
	contamination_seq_file = output_dir + "contamination_seq_file.fasta"
	seqs_kept_without_blast_hits_file = output_dir + "seqs_kept_without_blast_hits_file.fasta"
	
	plant_seqs = {}
	contamination_seqs = {}
	uncertain_seqs = {}
	
	os.popen( "cp " + assembly_file + " " + prefix )
	assembly_file = prefix + assembly_file.split('/')[-1]
	
	# --- run BLASTn against close relative to identify all high quality contigs --- #
	if len( plant_ref_file ) > 3:
		plant_ref_db = prefix + "plant_ref_db"
		if active:
			os.popen( "makeblastdb -in " + plant_ref_file + " -out " +plant_ref_db + " -dbtype nucl"  )
		remaining_query_file = prefix + "remaining_query_file.fasta"
		final_result_file = output_dir + "BLASTn_vs_plant_ref.txt"
		if active:
			run_blastn_on_cluster( assembly_file, plant_ref_db, final_result_file, para_jobs )
		high_qual_blast_results = load_high_quality_blast_results( final_result_file )
		all_seqs = load_sequences( assembly_file )
		plant_seq_counter = 0
		uncertain_seq_counter = 0
		with open( remaining_query_file, "w" ) as out:
			for key in all_seqs.keys():
				try:
					high_qual_blast_results[ key ]
					plant_seq_counter += 1
					plant_seqs.update( { key: all_seqs[ key ] } )
				except KeyError:
					out.write( '>' + key + '\n' + all_seqs[ key ] + '\n' )
					uncertain_seq_counter += 1
		print "number of reference-based identified plant contigs: " + str( plant_seq_counter )
		print "number of remaining uncertain contigs: " + str( uncertain_seq_counter )
	else:
		remaining_query_file = assembly_file
	
	
	# --- run BLASTn vs nt to classify the remaining contigs --- #
	final_blastn_result_file = output_dir + "BLASTn_vs_nt.txt"
	if active:
		run_blastn_on_cluster( remaining_query_file, "nt", final_blastn_result_file, para_jobs )	
	blast_results = load_blast_results( final_blastn_result_file )
	print "number of BLAST hits: " + str( len( blast_results.keys() ) )
	all_seqs = load_sequences( remaining_query_file )
	print "number of sequences without hit: " + str( len( all_seqs.keys() ) - len( blast_results.keys() ) )
	
	plant_seq_ids, contamination_seq_ids = filter_contamination_sequences( blast_results, all_seqs,  plant_table_file, prokaryote_table_file, gi_look_up_table_file )
	print "number of plant sequence IDs: " + str( len( plant_seq_ids.keys() ) )
	print "number of contamination sequence IDs: " + str( len( contamination_seq_ids.keys() ) )
	for key in all_seqs.keys():
		try:
			plant_seq_ids[ key ]
			plant_seqs.update( { key: all_seqs[ key ] } )
		except KeyError:
			try:
				contamination_seq_ids[ key ]
				contamination_seqs.update( { key: all_seqs[ key ] } )
			except KeyError:
				uncertain_seqs.update( { key: all_seqs[ key ] } )
	
	
	# --- check with black list --- #
	if len( black_list ) > 3:
		black_ids = []
		with open( black_list, "r" ) as f:
			line = f.readline()
			while line:
				black_ids.append( line.strip() )
				line = f.readline()
		print "numbers of IDs on black list: " + str( len( black_ids ) )
		for ID in black_ids:
			try:
				del uncertain_seqs[ ID ]
				contamination_seqs.update( { ID: all_seqs[ ID ] } )
			except KeyError:
				del plant_seqs[ ID ]
				contamination_seqs.update( { ID: all_seqs[ ID ] } )
	
	
	# --- generate output files --- #
	with open( cleaned_assembly_file, "w" ) as out:
		for key in sorted( plant_seqs.keys() ):
			out.write( '>' + key + '\n' + plant_seqs[ key ] + '\n' )
		for key in sorted( uncertain_seqs.keys() ):
			out.write( '>' + key + '\n' + uncertain_seqs[ key ] + '\n' )
	print "final number of sequences in assembly: " + str( len( plant_seqs.keys() + uncertain_seqs.keys() ) )
	
	with open( contamination_seq_file, "w" ) as out:
		for key in sorted( contamination_seqs.keys() ):
			out.write( '>' + key + '\n' + contamination_seqs[ key ] + '\n' )
	
	with open( seqs_kept_without_blast_hits_file, "w" ) as out:
		for key in sorted( uncertain_seqs.keys() ):
			out.write( '>' + key + '\n' + uncertain_seqs[ key ] + '\n' )


if __name__ == '__main__':
	
	if '--assembly_file' in sys.argv and '--plant_ref_file' in sys.argv and '--tmp_cluster_dir' in sys.argv and '--final_result_dir' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
	
	print "all done!"

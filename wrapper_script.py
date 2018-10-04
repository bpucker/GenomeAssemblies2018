### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

import sys, glob, re, os, time, datetime, shutil

# --- end of imports --- #

def read_information_file( info_file ):
	"""! @brief load all information from given file """
	
	information = {}
	
	with open( info_file, "r" ) as f:
		line = f.readline()
		while line:
			#x = line[5]	#just to check file for content
			line = f.readline()
	return information


def submit_jobs_to_cluster( assembly_file_names, cluster_dir_names, para_jobs, script_name ):
	"""! @brief submit BLAST jobs for each file to cluster """
	
	IDs_to_check = []
	batch_ID = str( datetime.datetime.now() )[-3:]
	for idx, file_name in enumerate( assembly_file_names ):
		ID = "E_" + batch_ID + '_' + str( idx ).zfill(4)
		IDs_to_check.append( ID )
		sh_file = cluster_dir_names[ idx ] + ID + '.sh'
		out_file = cluster_dir_names[ idx ] + ID + '.out'
		err_file = cluster_dir_names[ idx ] + ID + '.err'
		
		cmd = "python " + script_name + " --in " + assembly_file_names[ idx ] + " --cluster_dir " +  cluster_dir_names[ idx ]
		
		with open( sh_file, "w" ) as out:
				out.write( "#!/bin/bash\n" + " ".join( [ 	"echo " + '"',
																cmd + '"',
																"| qsub -cwd",
																"-pe multislot 8",
																"-N",
																ID,
																"-l vf=5G",
																"-l arch=lx-amd64",
																"-P fair_share",
																"-l idle=1",
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
			qstat_IDs = re.findall( "E_" + batch_ID + "_\d{4}", content )
			counter = 0
			for ID in qstat_IDs:
				if ID in IDs_to_check:
					counter += 1
			if counter < para_jobs:
				waiting_status = False
			else:
				time.sleep( 1 )
	
	waiting_status = True
	while waiting_status:
		qstat = os.popen( "qstat" )
		content = qstat.read()
		qstat_IDs = re.findall( "E_" + batch_ID + "_\d{4}", content )
		waiting_status = False
		for ID in IDs_to_check:
			if ID in qstat_IDs:
				waiting_status = True
		time.sleep( 120 )


if __name__ == '__main__':
	
	input_dir = "/genome_assemblies_of_interest/"	#INPUT DIRECTORY
	cluster_prefix = "/genome_assembly_evaluation/"	#OUTPUT DIRECTORY
	
	if not cluster_prefix[ -1 ] == '/':
		cluster_prefix += "/"
	if not os.path.exists( cluster_prefix ):
		os.makedirs( cluster_prefix )
	
	para_jobs = 50
	script_name = "run_evaluation_on_assembly.py"
	
	input_assembly_files = glob.glob( input_dir + "*.fasta" )
	active_files = []
	for filename in input_assembly_files:
		try:
			information = read_information_file( filename.replace( '.fasta', '.info' ) )
		except:
			active_files.append( filename )
		
	# --- run everything for all active files --- #
	cluster_dir_names = []
	for idx, filename in enumerate( active_files ):
		dir_name = cluster_prefix + filename.split('/')[-1].split('.')[0] + str( idx ).zfill(3)+ "/"
		if not os.path.exists( dir_name ):
			os.makedirs( dir_name )
		cluster_dir_names.append( dir_name )
	
	print "number of jobs to submit: " + str( len( active_files ) )
	for each in active_files:
		print each
	time.sleep( 10 )
	submit_jobs_to_cluster( active_files, cluster_dir_names, para_jobs, script_name )
	
	print "all done!"

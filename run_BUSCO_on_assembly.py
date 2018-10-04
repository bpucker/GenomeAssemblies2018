### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

__usage__ = """
					python run_BUSCO_on_assembly.py\n
					--cluster_dir <FULL_PATH_TO_OUTPUT_DIR>
					--in <INPUT_ASSEMBLY_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


import sys, os, re, glob
from operator import itemgetter

# --- end of imports --- #

def run_BUSCO( input_file, prefix, busco_path, augustus_path, augustus_config_path ):
	"""! @brief run BUSCO in genome mode on given assembly file """
	
	os.chdir( prefix )
	os.environ["PATH"] = augustus_path + ":" + os.environ["PATH"]
	print os.environ["PATH"]
	os.environ["AUGUSTUS_CONFIG_PATH"] = augustus_config_path
	print os.environ["AUGUSTUS_CONFIG_PATH"]
	cmd = "python " + busco_path + " --in " + input_file + " --out busco_run > " + prefix +"log.txt"
	os.popen( cmd )


def main( parameters ):
	"""! @brief run all genome analysis methods """
	
	working_dir = parameters[ parameters.index('--cluster_dir')+1 ]
	input_file = parameters[ parameters.index('--in')+1 ]
	
	augustus = "augustus-3.2"
	augustus_seqs_ex_script = "getAnnoFasta.pl"
	
	busco = "BUSCO/scripts/run_BUSCO.py"
	augustus_path = "augustus-3.2.2/bin/"
	augustus_config_path = "augustus-3.2.2/config/"
	
	
	if working_dir[-1] != '/':
		working_dir += "/"
	if not os.path.exists( working_dir ):
		os.makedirs( working_dir )
	
	
	# --- transfer file to cluster --- #
	assembly_file = working_dir + input_file.split('/')[-1]
	os.popen( "cp " + input_file + " " + assembly_file )
		
	# --- run BUSCO v3 --- #
	prefix = working_dir + "busco/"
	if not os.path.exists( prefix ):
		os.makedirs( prefix )
	run_BUSCO( input_file, prefix, busco, augustus_path, augustus_config_path )
		
		
if __name__ == '__main__':
	
	if '--cluster_dir' in sys.argv and '--in' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )

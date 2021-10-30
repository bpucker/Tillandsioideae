### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###

__usage__ = """
					python3 collect_single_copy_seqs_per_spec.py
					--in <BUSCO_RESULT_FOLDER>
					--out <OUTPUT_FOLDER>
					"""


import os, sys, subprocess, glob

# --- end of imports --- #

def main( arguments ):
	
	BUSCO_input_folder = arguments[ arguments.index('--in')+1 ]
	BUSCO_output_folder = arguments[ arguments.index('--out')+1 ]
	
	if BUSCO_input_folder[-1] != "/":
		BUSCO_input_folder += "/"
	if BUSCO_output_folder[-1] != "/":
		BUSCO_output_folder += "/"
	
	busco_result_files = glob.glob( BUSCO_input_folder + "*/busco/run_busco_run/full_table_busco_run.tsv" )
	
	if not os.path.exists( BUSCO_output_folder ):
		os.makedirs( BUSCO_output_folder )
	
	print ( "number of BUSCO result files: " + str( len( busco_result_files ) ) )
	for filename in busco_result_files:
		p = subprocess.Popen( args="cp " + filename + " " + BUSCO_output_folder + filename.split('/')[-4] + ".busco_full_table.tsv", shell=True )
		p.communicate()
			


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

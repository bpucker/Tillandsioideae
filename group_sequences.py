### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.1 ###

__usage__ = """
					python3 group_sequences.py
					--busco <BUSCO_RESULT_FOLDER>
					--pep <PEP_FILE_FOLDER>
					--cds <CDS_FILE_FOLDER>
					--out <OUTPUT_FOLDER>
					
					optional:
					--minspec <MINIMAL_NUMBER_OF_REPRESENDED_SPECS_PER_GENE>
					"""


import os, sys, subprocess, glob, re

# --- end of imports --- #


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences



def main( arguments ):
	"""! @brief run everything """
	
	busco_result_folder = arguments[ arguments.index('--busco')+1 ]
	pep_sequence_folder = arguments[ arguments.index('--pep')+1 ]
	cds_sequence_folder = arguments[ arguments.index('--cds')+1 ]
	aln_input_file_folder = arguments[ arguments.index('--out')+1 ]
	
	if '--minspec' in arguments:
		min_spec_cutoff = int( arguments[ arguments.index('--minspec')+1 ] )
	else:
		min_spec_cutoff = 100	#minimal number of species having a homologous sequence in the assembly
	selection_string = "cluster"	#restrict to one type of dataset
	
	
	cds_aln_input = aln_input_file_folder + "cds/"
	pep_aln_input = aln_input_file_folder + "pep/"
	
	if not os.path.exists( cds_aln_input ):
		os.makedirs( cds_aln_input )
	if not os.path.exists( pep_aln_input ):
		os.makedirs( pep_aln_input )
	
	# --- load IDs --- #
	busco_result_tables = glob.glob( busco_result_folder + "*.busco_full_table.tsv" )
	seq_IDs_for_tree = {}	#{ Busco_ID1: { 'spec1': seqID, 'spec2': seqID, ... }, Busco_ID2: { 'spec1': seqID, 'spec2': seqID, ... }, ... }
	for result_file in busco_result_tables:
		ID = result_file.split('/')[-1].split('.')[0]
		if selection_string in ID or "_" not in ID:
			candidates = {}
			with open( result_file, "r" ) as f:
				line = f.readline()
				while line:
					if line[0] != "#":
						parts = line.strip().split('\t')
						if parts[1] == "Complete":
							try:
								candidates[ parts[0] ].append( parts[2] )
							except KeyError:
								candidates.update( { parts[0]: [ parts[2] ] } )
					line = f.readline()
			for key in candidates.keys():
				if len( candidates[ key ] ) == 1:
					try:
						seq_IDs_for_tree[ key ].update( { ID: candidates[ key ][0] } )
					except KeyError:
						seq_IDs_for_tree.update( { key: { ID: candidates[ key ][0] } } )
	
	
	# --- load all peptides --- #
	peptide_collection = {}
	pep_files = glob.glob( pep_sequence_folder + "*.pep.fasta" )
	for pep_file in pep_files:
		ID = pep_file.split('/')[-1].split('.')[0]
		if selection_string in ID or "_" not in ID:
			peps = load_sequences( pep_file )
			peptide_collection.update( { ID: peps } )
	
	# --- load all CDS --- #
	cds_collection = {}
	cds_files = glob.glob( cds_sequence_folder + "*.cds.fasta" )
	for cds_file in cds_files:
		ID = cds_file.split('/')[-1].split('.')[0]
		if selection_string in ID or "_" not in ID:
			cdss = load_sequences( cds_file )
			cds_collection.update( { ID: cdss } )
	
	
	# --- construct alignment input files --- #
	for key in seq_IDs_for_tree.keys():
		if len( seq_IDs_for_tree[ key ] ) > min_spec_cutoff:
			pep_output_file = pep_aln_input + key + ".pep.fasta"
			with open( pep_output_file, "w" ) as out:
				for spec in seq_IDs_for_tree[ key ].keys():
					seq_id = seq_IDs_for_tree[ key ][ spec ]
					seq = peptide_collection[ spec ][ seq_id ]
					out.write( '>' + spec + "@" + seq_id + "\n" + seq + "\n" )
			cds_output_file = cds_aln_input + key + ".cds.fasta"
			with open( cds_output_file, "w" ) as out:
				for spec in seq_IDs_for_tree[ key ].keys():
					seq_id = seq_IDs_for_tree[ key ][ spec ]
					seq = cds_collection[ spec ][ seq_id ]
					out.write( '>' + spec + "@" + seq_id + "\n" + seq + "\n" )


if '--busco' in sys.argv and '--pep' in sys.argv and '--cds' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

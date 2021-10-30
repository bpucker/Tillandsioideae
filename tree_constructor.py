### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###
### v0.11 ###

__usage__ = """
					python3 tree_constructor.py
					--aln <ALIGNMENT_FILE>
					--out <OUTPUT_FOLDER>
					"""


import os, sys, subprocess, glob, re

# --- end of imports --- #


def load_alignment( aln_file, tmp_mapping ):
	"""! @brief load alignment and replace query IDs by real sequence names """
	
	sequences = {}
	
	with open( aln_file ) as f:
		header = f.readline()[1:].strip()
		try:
			header = tmp_mapping[ header ]
		except KeyError:
			pass
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					try:
						header = tmp_mapping[ header ]
					except KeyError:
						pass
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def alignment_trimming( aln_file, cln_aln_file, occupancy ):
	"""! @brief remove all alignment columns with insufficient occupancy """
	
	alignment = load_alignment( aln_file, {} )
	# --- if there is an alignment (expected case) 
	if len( list(alignment.keys()) ) > 0:
		# --- identify valid residues in aligned sequences (columns with sufficient occupancy) --- #
		valid_index = []
		for idx, aa in enumerate( list(alignment.values())[0] ):
			counter = 0
			for key in list(alignment.keys()):
				if alignment[ key ][ idx ] != "-":
					counter += 1
			if counter / float( len( list(alignment.keys()) ) ) > occupancy:
				valid_index.append( idx )
		
		# --- generate new sequences --- #
		with open( cln_aln_file, "w" ) as out:
			for key in list(alignment.keys()):
				seq = alignment[ key ]
				new_seq = []
				for idx in valid_index:
					new_seq.append( seq[ idx ] )
				out.write( ">" + key + '\n' + "".join( new_seq ) + '\n' )
	# --- just in case the alignment file is empyt (is this possible?) ---#
	else:
		with open( cln_aln_file, "w" ) as out:
			out.write( "" )


def special_alignment_load( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if "@" in header:
			header = header.split('@')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
				seq = "".join( seq ).replace("\x00", "-")
				if not "\x00" in seq:
					sequences.update( { header: seq } )
				else:
					print ( header)
					print ( [ seq ])
				header = line.strip()[1:]
				if "@" in header:
					header = header.split('@')[0]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		seq = "".join( seq ).replace("\x00", "-")
		if not "\x00" in seq:
			sequences.update( { header: seq } )	
		else:
			print ( header)
			print ( [ seq ] )
	return sequences


def load_mapping_table( name_mapping_file ):
	"""! @brief load name mapping table """
	
	name_mapping_table = {}
	with open( name_mapping_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			name_mapping_table.update( { parts[0]: parts[1] + "_" + parts[2] } )
			line = f.readline()
	return name_mapping_table


def main( arguments ):
	"""! @brief run everything """
	
	aln_input_file_folder = arguments[ arguments.index('--aln')+1 ]
	output_folder = arguments[ arguments.index('--out')+1 ]
	
	if '--mode' in arguments:
		mode = arguments[ arguments.index('--mode')+1 ]
		if mode not in [ "fasttree", "raxml" ]:
			mode = "fasttree"
	else:
		mode = "fasttree"
	
	if '--raxml' in arguments:
		raxml = arguments[ arguments.index('--raxml')+1 ]
	else:
		raxml = "raxml"
	
	if "--fasttree" in arguments:
		fasttree = arguments[ arguments.index('--fasttree')+1 ]
	else:
		fasttree = "fasttree"
	
	mafft = "mafft"
	pxaa2cdn = "pxaa2cdn"
	
	if '--mapping' in arguments:
		name_mapping_file = arguments[ arguments.index('--mapping')+1 ]
	else:
		name_mapping_file = ""
	
	cpu = 4
	num_of_alignments = 10
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	files_per_gene = {}
	pep_files = glob.glob( aln_input_file_folder + "pep/*.fasta" )
	for pep in pep_files:
		files_per_gene.update( { pep.split('/')[-1].split('.')[0]: { 'pep': pep } } )
	cds_files = glob.glob( aln_input_file_folder + "cds/*.fasta" )
	for cds in cds_files:
		try:
			files_per_gene[ cds.split('/')[-1].split('.')[0] ].update( { 'cds': cds } )
		except KeyError:
			print ( "ERROR: corresponding pep file missing - " + cds )
	
	aln_collection = {}	#alignments: { 'gene1': { 'spec1': "AC---GT", 'spec2': "--GG----" }, 'gene2': { ... }, ... }
	species_list = []
	aln_len_per_gene = {}	#length of clean alignment per gene
	for gene in files_per_gene.keys():
		try:
			pep_file = files_per_gene[ gene ]['pep']
			cds_file = files_per_gene[ gene ]['cds']
			
			# --- construct inial alignment --- #
			aln_file = output_folder + gene + ".fasta.aln"
			if not os.path.isfile( aln_file ):
				p = subprocess.Popen( args= mafft + " --quiet " + pep_file + " > " + aln_file, shell=True )
				p.communicate()
						
			
			# --- translate into nucleotides --- #
			cds_alignment_file = aln_file + ".cds_aln"
			if not os.path.isfile( cds_alignment_file ):
				p = subprocess.Popen( args=" ".join( [ pxaa2cdn, "-a", aln_file, "-n", cds_file, "-o", cds_alignment_file ] ), shell=True )
				p.communicate()
			
			# --- trim alignment based on occupancy --- #
			cln_aln_file = cds_alignment_file + ".cln"
			if not os.path.isfile( cln_aln_file ):
				alignment_trimming( cds_alignment_file, cln_aln_file, occupancy=0.1 )
			
			alignment = special_alignment_load( cln_aln_file )
			aln_collection.update( { gene:  alignment } )
			species_list += alignment.keys()
			aln_len_per_gene.update( { gene: len( list( alignment.values() )[0] ) } )
			
		except KeyError:
			print ( "ERROR: corresponding cds file missing" -  files_per_gene[ gene ]['pep'] )
	
	
	# --- merge all clean alignments and fill species gaps --- #
	species_list = list( sorted( list( set( species_list ) ) ) )	#sorted list of all species of interest
	genes = list( sorted( aln_collection.keys() ) )		#sorted list of all genes for tree construction
	data_per_species = {}
	for spec in species_list:
		data_per_species.update( { spec: [] } )
	
	for gene in genes[:min( [ num_of_alignments, len( genes )-1 ] )]:	#only use given number of alignments for tree construction
		for spec in species_list:
			try:
				seq = aln_collection[ gene ][ spec ].upper()
			except KeyError:
				seq = "-" * aln_len_per_gene[ gene ]
			data_per_species[ spec ].append( seq.strip() )
	
	# --- construct tree --- #
	final_alignment_output_file = output_folder + "final_aln_file.fasta.aln"
	with open( final_alignment_output_file, "w" ) as out:
		for spec in species_list:
			out.write( '>' + spec.replace( "_cluster", "" ).replace( "_trinity", "" ).replace( "_unigene", "" ) + "\n" + "".join( data_per_species[ spec ] ) + "\n" )
	
	cln_final_alignment_output_file = final_alignment_output_file + ".cln"
	if not os.path.isfile( cln_final_alignment_output_file ):
		alignment_trimming( final_alignment_output_file, cln_final_alignment_output_file, occupancy=0.1 )
	
	if mode == "raxml":	#RAxML
		prefix = output_folder + "RAxML_tree"
		fin_tree_file = prefix + ".raxml.bestTree"
		if not os.path.isfile( fin_tree_file ):
			p = subprocess.Popen( args= " ".join( [ raxml, "--all --threads " + str( cpu ) + " --model GTR+F+G --msa", cln_final_alignment_output_file, "--prefix", prefix ] ), shell=True )
			p.communicate()
	else:	#FastTree2
		fin_tree_file = output_folder + "FastTree_tree.tre"
		if not os.path.isfile( fin_tree_file ):
			p = subprocess.Popen( args= " ".join( [ fasttree, "-gtr -nt -nopr -nosupport <", cln_final_alignment_output_file, ">", fin_tree_file ] ), shell=True )
			p.communicate()
	
	# --- load mapping table and replace names in tree files --- #
	if len( name_mapping_file ) > 0:
		renamed_fin_tree_file = fin_tree_file.replace( ".tre", ".renamed.tre" )
		if not os.path.isfile( renamed_fin_tree_file ):
			mapping_table = load_mapping_table( name_mapping_file )
			with open( fin_tree_file, "r" ) as f:
				tree = f.read()
			IDs_to_replace = list( sorted( re.findall( "[a-zA-Z]{1,3}[a-zA-Z0-9]{1,5}", tree ) ) )[::-1]
			print ( "number of leaves in tree: " + str( len( IDs_to_replace ) ) )
			for ID in IDs_to_replace:
				try:
					tree = tree.replace( ID, mapping_table[ ID ] )
				except KeyError:
					print ( "name replacement in tree failed: " + ID )
			with open( renamed_fin_tree_file, "w" ) as out:
				out.write( tree )


if '--aln' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

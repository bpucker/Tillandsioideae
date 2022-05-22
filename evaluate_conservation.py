### Boas Pucker ###
### b.pucker@tu-braunschweig.de ###

__version__ = "v0.1"

__usage__ = """
					python3 MYB_annotator.py
					--ids <ID_FILE>
					--orthogroups <ORTHOGROUP_TSV_FILE>
					--reference <REFERENCE_NAME>
					--results <RESULT_FILE>
					
					bug reports and feature requests: b.pucker@tu-braunschweig.de
					"""


import os, sys
import numpy as np

# --- end of imports --- #

def quantify_conservation( dataset ):
	"""! @brief quantify conservation of given gene """
	
	samples = list( dataset.keys() )
	counter = 0.0
	for sample in samples:
		if len( dataset[ sample ] ) > 0:
			counter += 1
	return ( counter / len( samples ) )


def check_presence( dataset ):
	"""! @brief check presence of ortholog across all samples """
	
	samples = list( dataset.keys() )
	names = []
	for sample in samples:
		if len( dataset[ sample ] ) > 0:
			names.append( sample )
	return names
	

def main( arguments ):
	"""! @brief run everything """
	
	ref_enriched_geneID_file = arguments[ arguments.index('--ids')+1 ]	#text file with one ID per line
	orthogroup_file = arguments[ arguments.index('--orthogroups')+1 ]	#Orthogroups.tsv
	reference = arguments[ arguments.index('--reference')+1 ]	#LYUTD
	result_file = arguments[ arguments.index('--results')+1 ]

	# --- load genes of interest --- #
	genes_of_interest = {}
	with open( ref_enriched_geneID_file, "r" ) as f:
		line = f.readline()
		while line:
			genes_of_interest.update( { line.strip(): None } )
			line = f.readline()

	# --- load data into dictionary with ref genes as keys --- #
	orthologs_per_ref_gene = {}
	with open( orthogroup_file, "r" ) as f:
		raw_headers = f.readline().strip().replace(".pep", "").split('\t')[1:]
		headers = []
		for each in raw_headers:
			if "_cluster" in each:
				headers.append( each.replace( "_cluster", "" ) )
			else:
				headers.append( each )
		ref_index = headers.index( reference )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			try:
				ref_genes = parts[ ref_index+1 ]
				if len( ref_genes ) > 3:
					if ", " in ref_genes:
						ref_genes = ref_genes.split(', ')
					else:
						ref_genes = [ ref_genes ]
					data = {}
					for idx, sample in enumerate( headers ):
						try:
							genes = parts[ idx +1 ]
							if len( genes ) > 3:
								if ", " in genes:
									genes = genes.split(', ')
								else:
									genes = [ genes ]
							else:
								genes = []
						except IndexError:
							genes = []
						data.update( { sample: genes } )
					for rg in ref_genes:
						orthologs_per_ref_gene.update( { rg: data } )
				else:	#no ref gene in this orthogroup (empty string)
					pass
			except IndexError:	#no ref gene in this orthogroup (line shorter due to strip() )
				pass
			line = f.readline()


	# --- check conservation per gene --- #
	conservation_genes_of_interest = []
	conservation_background = []
	for gene in list( orthologs_per_ref_gene.keys() ):
		cons = quantify_conservation( orthologs_per_ref_gene[ gene ] )
		try:
			genes_of_interest[ gene ]
			conservation_genes_of_interest.append( cons )
		except KeyError:
			conservation_background.append( cons )


	print( "median conservation genes of interest: " + str( np.median( conservation_genes_of_interest ) ) )
	print( "median conservation background genes: " + str( np.median( conservation_background ) ) )

	print( "mean conservation genes of interest: " + str( np.mean( conservation_genes_of_interest ) ) )
	print( "mean conservation background genes: " + str( np.mean( conservation_background ) ) )


	# --- check presence of orthologs across all samples --- #
	with open( result_file, "w" ) as out:
		for gene in list( orthologs_per_ref_gene.keys() ):
			names = check_presence( orthologs_per_ref_gene[ gene ] )
			out.write( gene + "\t" + ";".join( names ) + "\n" )


if '--ids' in sys.argv and '--orthogroups' in sys.argv and '--reference' in sys.argv and '--results' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )

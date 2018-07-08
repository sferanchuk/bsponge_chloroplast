#!/usr/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
from Bio import Phylo


selgenes = [ "rrs", "atpB" ]
sgstat = { "rrs" : [ 1499, 18 ], "psaB" : [ 1790, 21 ], "tufA" : [ 1225, 13 ], "atpB" : [3510, 66 ] , "rpoB" : [ 3765, 56 ] }
#rnconv = { "p9" : "tufA", "p34" : "rpl19", "p4" : "psaB", "p6" : "atpB", "p2" : "rpoB" }
fig, ax = plt.subplots( 1, len( selgenes ), figsize=(40,20), dpi=80 )

for k in range( len( selgenes ) ):
	gene = selgenes[k]
	tree = Phylo.read( "fragments_" + gene + ".nwk", 'newick' )
	matplotlib.rcParams["font.size"] = 28
	matplotlib.rcParams["lines.linewidth"] = 3
	leafs = tree.get_terminals( order="postorder" )

	bspnum = 0
	spid = "bsponge"
	labelconv = { "p_salinaru" : "Picocystis salinarum", "myrmecia" : "Myrmecia israeliensis", "botyococcu" : "Botryococcus braunii", "coccomyxa" : "Coccomyxa subellipsoidea", "c_parasiti" : "Choricystis parasitica", "bsponge" : "Sponge symbiont", "h_reticala" : "Hydrodictyon reticulatum", "m_jurisii" : "Mychonastes jurisii", "chlorella" : "Chlorella vulgaris" }
	sleafs = [ str( leaf ) for leaf in leafs ]
	spnum = sleafs.index( spid )
	spdepth = tree.distance( leafs[ spnum ] )
	scmap = {}
	for sleaf in sleafs:
		if sleaf in labelconv:
			scmap[ labelconv[ sleaf ] ] = "#5d3239"
	scmap[ labelconv[ spid ] ] = "#251923"
	def labelfunc( csn ):
		sn = str( csn )
		if sn != 'Clade':
			if sn in labelconv:
				return labelconv[ sn ]
			return sn
		return ''

	Phylo.draw( tree, axes=ax[ k ], label_func = labelfunc, label_colors = scmap  )	
	ax[ k ].set_title( "Gene: " + gene, loc = "left", size=32 )
	ax[ k ].set_xlim( [ -0.1, 0.42 ] )
	ax[ k ].set_xticks( [ 0, 0.1, 0.2, 0.3 ] )
	ax[ k ].set_xticklabels( [ "0.0", "0.1", "0.2", "0.3" ], size=22 )
	ax[ k ].set_yticks( [] )
	ax[ k ].set_ylabel( "" )
	ax[ k ].set_xlabel( "", size=20 )
	ax[ k ].set_ylim( [ 0, len( sleafs ) + 1 ] )
	bwidth = float( sgstat[ gene ][ 1 ] ) / sgstat[ gene ][ 0 ]
	ax[ k ].barh( spnum + 1, [ bwidth ], 0.4, spdepth - bwidth, color = "#f8a28d", edgecolor = "#251923", lw = 2 )
plt.savefig( "trees_chart.png" )
plt.savefig( "trees_chart.eps" )
	


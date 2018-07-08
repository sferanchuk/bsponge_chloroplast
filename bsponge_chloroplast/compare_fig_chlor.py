
from reportlab.lib import colors
from reportlab.lib.colors import red, grey, orange, green, brown
from reportlab.lib.colors import blue, lightblue, purple

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

seq1 = SeqIO.parse( "bchchlor_cds.gb", "gb")
seq2 = SeqIO.parse( "chori_chlor.gb", "gb")

name = "chlor_diag"
gd_diagram = GenomeDiagram.Diagram(name)
feature_sets = {}

gdtrack1 = gd_diagram.new_track( 1, name="n1")
gdtrack2 = gd_diagram.new_track( 2, name="n2")
gdtrack1.height = 0.2
gdtrack1.end = 55838
gdtrack1.name = "Sponge symbiont"
gdtrack1.axis_labels = True
gdtrack2.height = 0.2
gdtrack2.end = 94206

rotate = 66000 #79000
revrotate = gdtrack2.end - rotate

feature_sets[ "r1" ] = gdtrack1.new_set()
feature_sets[ "r2" ] = gdtrack2.new_set()

selgenes = [ "rrs", "psaB", "tufA", "atpB", "rpoB" ]
	
rnconv = { "p9" : "tufA", "p34" : "rpl19", "p4" : "psaB", "p6" : "atpB", "p2" : "rpoB" }
rbars1 = {}
if True:
	for seq_record in seq1 :
		for seq_feature in seq_record.features :
			if seq_feature.type in [ "CDS", "cdna", "rRNA", "tRNA" ] :
				#assert len(seq_feature.qualifiers['translation'])==1
				fname = seq_feature.qualifiers['ID'][0]
				if fname in rnconv:
					fname = rnconv[ fname ]
				flabel = True if fname in selgenes else False
				fcolor = colors.Color( 0.36, 0.19, 0.22, 0.65 ) if seq_feature.type in [ "cdna", "tRNA", "rRNA" ] else colors.Color( 0.96, 0.63, 0.55, 0.95)
				fangle = 0 if seq_feature.location.strand == 1 else 180
				rbars1[ fname ] = feature_sets[ "r1" ].add_feature( seq_feature, sigil="BOX", color=fcolor, label=flabel, name=fname, label_position="start", label_size=14, label_angle=fangle)


rbars2= {}
if True:
	for seq_record in seq2 :
		for seq_feature in seq_record.features :
			if seq_feature.type in [ "CDS", "rRNA", "tRNA" ] :
				#assert len(seq_feature.qualifiers['translation'])==1
				if seq_feature.location.start > rotate:
					seq_feature.location = FeatureLocation( seq_feature.location.start - rotate, seq_feature.location.end - rotate, seq_feature.location.strand )
				else:
					seq_feature.location = FeatureLocation( seq_feature.location.start + revrotate, seq_feature.location.end + revrotate, seq_feature.location.strand )
				fname = seq_feature.qualifiers['gene'][0]
				flabel = True if fname in selgenes else False
				fcolor = colors.Color( 0.36, 0.19, 0.22, 0.65 ) if seq_feature.type in [ "rRNA", "tRNA" ] else colors.Color( 0.46, 0.52, 0.83, 0.95 )
				fangle = 0 if seq_feature.location.strand == 1 else 180
				rbars2[ fname ] = feature_sets[ "r2" ].add_feature( seq_feature, sigil="BOX", color=fcolor, label=flabel, name=fname, label_position="start", label_size=14, label_angle=fangle )
				print fname

color = colors.Color( 0.8, 0.8, 0.8, 0.9 )
border = colors.lightgrey

with open( "ch_comp.out" ) as f:
	for line in f:
		ls = line.split()
		id1 = ls[0][7:]
		if id1 in rnconv:
			id1 = rnconv[ id1 ]
		id2 = ls[1]
		print( id1, id2 )
		if id1 in rbars1 and id2 in rbars2:
			gd_diagram.cross_track_links.append( CrossLink( rbars1[ id1 ], rbars2[ id2 ], color, border))
			
for pair in [ ( "rrs", "rrs" ), ( "rrl", "rrl" ) ]:
	if pair[0] in rbars1 and pair[1] in rbars2:
			gd_diagram.cross_track_links.append( CrossLink( rbars1[ pair[0] ], rbars2[ pair[1] ], color, border))
	
print border

gd_diagram.draw(format="linear", pagesize='A4', fragments=1 )

#gd_diagram.write(name + ".pdf", "PDF")
gd_diagram.write( name + "_linear.png", "PNG")
gd_diagram.write( name + "_linear.eps", "EPS")
#gd_diagram.write(name + ".svg", "SVG")

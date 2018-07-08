

from Bio import SeqIO
from Bio.Seq import Seq
import os

rnconv = { "p9" : "tufA", "p34" : "rpl19", "p4" : "psaB", "p6" : "atpB", "p2" : "rpoB", "rrnS" : "rrs", "rrn16" : "rrs", "rrn16S" : "rrs" }
snconv = { "bchchlor_cds" : "bsponge", "chori_chlor" : "c.parasitica" }

for gene in [ "rrl", "rrf" ]: # [ "rrs", "atpB" ]: #"rpl19", "tufA", "psbA", "psaB", "atpB", "rpoB" ]: 
	fprefix = "fragments_" + gene
	output_handle = open( fprefix + ".fa", "w")
	for species in [ "botyococcus",  "chlorella",  "chori_chlor",  "coccomyxa", "myrmecia", "bchchlor_cds", "h.reticalatum",  "m.jurisii", "p.salinarum" ]:
		input_handle  = open( species + ".gb", "r")
		for seq_record in SeqIO.parse(input_handle, "genbank") :
			for seq_feature in seq_record.features :
				if seq_feature.type in [ "CDS", "rRNA", "tRNA" ] :
					if not 'gene' in seq_feature.qualifiers:
						if 'ID' in seq_feature.qualifiers:
							cgene = seq_feature.qualifiers['ID'][0]
						else:
							continue
					else:
						cgene = seq_feature.qualifiers['gene'][0]
					if cgene in rnconv:
						cgene = rnconv[ cgene ]
					if cgene != gene: 
						continue
					loc = seq_feature.location
					seq = Seq( str( seq_record.seq[ loc.start - 1 : loc.end ] ) )
					if loc.strand == -1:
						seq = seq.reverse_complement() 
					sname = species if not species in snconv else snconv[ species ]
					output_handle.write(">%s gene=%s\n%s\n" % ( sname, gene, str( seq )  ) )
					break
		input_handle.close()
	output_handle.close()
	os.system( "mafft --maxiterate 1000 --genafpair %s.fa >%s-align.fa" % ( fprefix, fprefix ) )
	os.system( "java -jar ###/readseq.jar -a -format=12 -inform=8 %s-align.fa -output=tmp.phylip" % fprefix )
	os.system( "fastme -m N -d J -i tmp.phylip -o %s.nwk" % fprefix ) 
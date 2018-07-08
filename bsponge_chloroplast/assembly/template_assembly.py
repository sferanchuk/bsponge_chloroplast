#!/usr/bin/python2

import sys
import os

templates = [ "chori_chlor" ]

samples = [ "1747", "1791", "1820" ]

scaf_path = "###"

riname = "stemplates2"
rname = riname + ".fa"
if not os.path.isfile( rname ):
	for sample in samples:
		spath = scaf_path + "/S.%s/S.%s.scafSeq.fa" % ( sample, sample )
		tspath = "sc%s.fa" % sample
		if not os.path.isfile( tspath ):
			os.system( "ln -s %s %s" % ( spath, tspath ) )
		if not os.path.isfile( tspath + ".nsq" ):
			os.system( "makeblastdb -dbtype nucl -in %s -parse_seqids" % tspath )
		idname = "st%s.ids" % sample
		idfile = open( idname, "w" )
		sids = set( [] )
		for template in templates:
			onames = [ sample + "_" + template + ".out", sample + "_" + template + ".nout" ]
			for m in range( 2 ):
				oname = onames[m]
				if not os.path.isfile( oname ):
					if m == 0:
						os.system( "tblastx -db %s -query %s.fa -evalue 1e-40 -out %s -outfmt \"6 sseqid qseqid evalue\" -max_target_seqs 2000" % ( tspath, template, oname ) )
					else:
						os.system( "blastn -task blastn -db %s -query %s.fa -evalue 1e-8 -out %s -outfmt \"6 sseqid qseqid evalue\" -max_target_seqs 2000" % ( tspath, template, oname ) )
				lsid = []
				with open( oname ) as f:
					for line in f:
						ls = line.split()
						evalue = float( ls[-1] )
						if evalue < 1e-70 or True:
							lsid.append( ls[0] )
				csids = set( lsid )
				sids |= csids
		for sid in sids:
			if len( sid ) > 4:
				idfile.write( sid + "\n" )
		idfile.close()
		os.system( "blastdbcmd -db " + tspath + " -outfmt %f -entry_batch " + idname + " >> " + rname )

if not os.path.isfile( riname + ".1.bt2" ):
	os.system( "bowtie2-build " + rname + " " + riname )
rdname = "ra2.fq"
if not os.path.isfile( rdname ):
	os.system( "rm b*.bam" )
	for sample in samples:
		rspath = "d" + sample 
		if not os.path.isfile( "b" + sample + ".bam" ):
			os.system( "bowtie2 -x %s -1 %s_1.fq -2 %s_2.fq -p 6 -k 1 --no-unal | samtools view -Sbh - >b%s.bam" % ( riname, rspath, rspath, sample ) )
		os.system( "samtools1.7 fastq -n b%s.bam > rd%s.fq" % ( sample, sample ) )
		os.system( "cat rd%s.fq >> %s" % ( sample, rdname ) )
		 
rtdname = "rta2.fq"
if not os.path.isfile( rtdname ):
	for template in templates:
		if not os.path.isfile( template + ".1.bt2" ):
			os.system( "bowtie2-build " + template + ".fa " + template )
		for sample in samples:
			rspath = "d" + sample
			if not os.path.isfile( "tb" + sample + ".bam" ):
				os.system( "bowtie2 -x %s -1 %s_1.fastq -2 %s_2.fastq -p 6 -k 1 --no-unal | samtools view -Sbh - >tb%s.bam" % ( template, rspath, rspath, sample ) )
			os.system( "/home/sferanchuk/soft/samtools-1.7/samtools fastq -n tb%s.bam > rtd%s.fq" % ( sample, sample ) )
			os.system( "cat rtd%s.fq >> %s" % ( sample, rtdname ) )




os.system( "cat ra2.fq rta2.fq >raa2.fq" )
rdname = "raa2.fq"

trpath = "###"

os.system( trpath + "/trinity-plugins/BIN/seqtk-trinity seq -A raa2.fq >raa2.fa" )
os.system( trpath + "/Inchworm/bin/inchworm --reads raa2.fa --run_inchworm -K 32 --DS --no_prune_error_kmers --min_assembly_coverage 1 >contigs.fa" )
	
			                                                                   
		
#!/usr/bin/python

import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import scipy.stats as stats
import math
import sys


rseqs = {}
if False:
	with open( "bch-cds.fa" ) as f:
		for line in f:
			if line[0] == ">":
				cseqid = line.strip()[1:]
				rseqs[ cseqid ] = ""
			else:
				rseqs[ cseqid ] += line.strip()

rrawseq = []
with open( "bchchlor.fa" ) as f:
	for line in f:
		if line[0] != '>':
			rrawseq += line.strip()
			
rsbegs = {}
rsends = {}
with open( "bchchlor.cds" ) as f:
	for line in f:
		ls = line.strip().split( "\t" )
		gname = ls[-1]
		pgpos = ls[-1].rfind( ".p" )
		if pgpos != -1:
			gname = ( ls[ -1 ][ pgpos + 1 : ] ).split( ";" )[0]
		gbeg = int( ls[3] )
		gend = int( ls[4] )
		rsbegs[ gbeg ] = gname
		rsends[ gbeg ] = gend
		rseqs[ gname ] = rrawseq[ gbeg : ( gend + 1 ) ]


rnconv = {}
ndic = {}		
with open( "ch_comp.out" ) as f:
	for line in f:
		ls = line.split()
		id1 = ls[0][7:]
		if id1 in rnconv:
			id1 = rnconv[ id1 ]
		id2 = ls[1]
		ndic[ id1 ] = id2

		
#genome = ""
#with open( "bchchlor.fa" ) as f:
#	for line in f:


gs = gridspec.GridSpec( 2, 2, height_ratios=[1, 4], width_ratios=[ 10, 1 ] )
fig = plt.figure( figsize=(40,13), dpi=80 )
ax = plt.subplot( gs[1,0] )
ax0 = plt.subplot( gs[0,0] )
ax1 = plt.subplot( gs[:,1] )

refs = sorted( rseqs.keys() )

#samples = [ "cg1s.bam", "cg2s.bam", "d1747s.bam", "d1791s.bam", "d1820s.bam", "r1747s.bam", "r1791s.bam" ]
samples = [ "d1747s.bam", "r1747s.bam", "d1791s.bam", "r1791s.bam", "d1820s.bam" ] 
snames = [ "DNA healthy", "RNA healthy", "DNA diseased", "RNA diseased", "DNA dead" ]
#pcolors = [ "#732c7b", "#bb3333", "#5555bb", "#bdaec6" ]
#pcolors = [ "#251923", "#5d3239", "#7687d5", "#f8a28d" ]
pcolors = [ "#5d3239", "#f8a28d", "#d1d9ff" ]
#epcolors = pcolors + [ "white" ]
aalet = [ "A", "C", "G", "T" ]
ierr = 0.01

def size_func( val ):
	if val == 0:
		return 0
	#return min( 0.4, 0.03 * math.log( val ) + 0.1 )	
	powf = math.exp( 0.1 * math.log( val ) )
	return min( 0.4, 0.17 * powf )	


def size_func_total( val ):
	if val == 0:
		return 0
	#val += 300
	powf = math.exp( 0.25 * math.log( val ) )
	print ( val, math.log( val ), math.sqrt( val ), powf )
	return min( 1.6, 0.2 * powf )	

tsize = [ 0 ] * len( samples )
tdistr = [ [ 0 ] * len( pcolors ) for k in xrange( len( samples ) ) ]

rskeys = sorted( rsbegs.keys() )
rsnames = []
ctg = "a1;411"
for rc in range( len( rskeys ) ):
	gbeg = rskeys[ rc ]
	ref = rsbegs[ gbeg ]
	gend = rsends[ gbeg ]
	if ref in ndic:
		rsnames.append( ndic[ ref ] )
	else:
		rsnames.append( ref )
						
	rseq = rseqs[ ref ]
	allspos = set( [] )
	alltpos = set( [] )
	allcdistr = []
	allsradius = []
	for k in range( len( samples ) ):
		sf = pysam.AlignmentFile( samples[k] )
		sreads = set( [] )
		spos = set( [] )
		tpos = set( [] )
		
		if len( rseq ) < gend - gbeg:
			print "mismatch %s %d %d %d" % ( ref, len( rseq ), gbeg, gend )
			gend = gbeg + len( rseq )
		nreads = 0
		ncreads = 0
		squal = 0.
		lsum = 0.
		for read in sf.fetch( ctg,  gbeg, gend ):
			#if snames[k][0] != "d" and not ( read.flag & 2 ):
			#	continue
			#if read.reference_id != 1:
			#	continue
			ap = read.get_aligned_pairs( matches_only=True )
			if len( ap ) == 0:
				continue
			#print ap
			s_eq = 0
			s_tot = 0
			#rseq = read.get_reference_sequence()
			qseq = read.query_alignment_sequence
			qqual = read.query_qualities
			qbeg = read.query_alignment_start
			rplist = read.get_reference_positions() 
			#print qqual
			cspos = set( [] )
			ctpos = set( [] )
			nreads += 1
			lsum += float( len( ap ) )
			cgbeg = ap[0][1] 
			cgend = ap[-1][1]
			crseq = rrawseq[ cgbeg  : cgend + 1 ]
			for cp in ap:
				if not qseq[ cp[0] ] in aalet:
					continue
				if qqual[ cp[0] ] < 30:
					continue
				if cp[1] < gbeg or cp[1] > gend:
					continue
				ctpos.add( cp[1] )
				s_tot += 1
				#if rseq[ cp[1] - gbeg ] == qseq[ cp[0] ]:
				if rrawseq[ cp[1] ] == qseq[ cp[0] ]:
					s_eq += 1
				else:
					cspos.add( cp[1] - gbeg )
			if s_tot > 0:
				qual = float( s_tot - s_eq ) / s_tot
				ncreads += 1
				squal += qual
				if qual < 0.01:
					sreads.add( read.query_name )
					spos |= cspos
					tpos |= ctpos
		alltpos |= tpos
		print "sample %s prot %s spos %d tpos %d sreads %d total ( %d %d %g %g %d )" % ( samples[k], ref, len( spos ), len(  tpos ), len( sreads ), nreads, ncreads, squal / max( ncreads, 1 ),  lsum / max( ncreads, 1 ), len( crseq ) )
		sf = pysam.AlignmentFile( samples[k] )
		pdistr = {}
		for read in sf.fetch( ctg, gbeg, gend ):
			if not read.query_name in sreads:
				continue
			ap = read.get_aligned_pairs( matches_only=True )
			qseq = read.query_alignment_sequence
			qqual = read.query_qualities
			for cp in ap:
				if not qseq[ cp[0] ] in aalet:
					continue
				if qqual[ cp[0] ] < 30:
					continue
				rind = cp[1] - gbeg
				if not rind in spos:
					continue
				nlet = aalet.index( qseq[ cp[0] ] )
				if not rind in pdistr:
					pdistr[ rind ] = [ 0 ] * 4
				pdistr[ rind ][ nlet ] += 1
		spoly = 0
		smut = 0
		spolymut = 0
		for pos in pdistr:
			nr = aalet.index( rseq[ pos ] )
			scdistr = sorted( pdistr[ pos ], reverse=True )
			nsuccess = scdistr[0]
			ntot = sum( scdistr )
			if pdistr[ pos ][ nr ] != nsuccess:
				smut += 1
			pvsum = 0
			pvcur = 0
			for ne in range( ntot - nsuccess + 1 ):
				lpv = math.log( 1 - ierr ) * ( ntot - ne ) + math.log( ierr ) * ( ne ) + math.lgamma( ntot + 1 ) - math.lgamma( ntot + 1 - ne ) - math.lgamma( ne + 1 )
				pvcur = math.exp( lpv )
				pvsum += pvcur
			if pvsum > 0.95:
				spoly += 1
				allspos.add( pos )
				if pdistr[ pos ][ nr ] != nsuccess:
					spolymut += 1
		print "poly %d mut %d polymut %d all %d tot %d tsize %d " % ( spoly, smut, spolymut, len( spos ), len( tpos ), len( rseq ) )
		mcdistr = [ 0, 0 ]
		if len( sreads ) > 0:
			#ax.pie( [ smut, spoly - smut, len( spos ) - spoly, len( tpos ) - len( spos ), len( rseq ) - len( tpos ) ], colors = [ "red", "blue", "cyan", "gray", "lightgray" ], radius = min( 0.3, 0.02 * math.log( len( sreads ) ) + 0.07 ), center = ( rc, 2 * k ) ) 
			cdistr = [ spolymut, spoly - spolymut ]
			for dc in range( len( cdistr ) ):
				tsize[ k ] += len( sreads )
				mcdistr[ dc ] = cdistr[dc] * float( len( rseq ) ) / len( tpos )
				tdistr[ k ][ dc ] += mcdistr[ dc ]
			print ( cdistr, mcdistr )
		allcdistr.append( mcdistr )
		allsradius.append( len( sreads ) )
	print "prot %s allspos %d ratio %g" % ( ref, len( allspos ), float( len( allspos ) ) / len( rseq ) )
	for k in range( len( samples ) ):
		if len( alltpos ) > 0:
			v2 = len( allspos ) * float( len( rseq ) ) / len( alltpos ) - sum( allcdistr[k] )
			ax.pie( allcdistr[k] + [ v2 ], colors = pcolors, radius = size_func( allsradius[k] ), center = ( rc, len( samples ) - k - 1 ) )
			tdistr[ k ][ 2 ] += v2
					
	#ax[ rc ].hist( distr, density = True, range = ranges[ rc ] )
	#ax[ rc ].set_title( refs[rc] )
ax.set_xlim( [ -1, len( rskeys ) ] )
ax.set_ylim( [ -1, len( snames ) ] )
ax.set_xticks( range( len( rskeys ) ) )
ax.set_xticklabels( rsnames, rotation='vertical', size = 20 )
ax.set_yticks( range( len( snames ) ) )
ax.set_yticklabels( reversed( snames ), size = 24 )
ax.set_aspect( "equal" )
scstep = len( rskeys ) / float( len( samples ) )
for sc in range( len( samples ) ):
	avgsize = tsize[ sc ] / len( rskeys )
	pradius = size_func_total( avgsize )
	ax0.pie( tdistr[ sc ], colors = pcolors, radius = pradius, center = ( sc * scstep, 1 ) )
	ax0.text( sc * scstep + pradius + 0.3, 1, snames[sc], ha="left", size=24 )
ax0.set_xlim( [ -1, len( rskeys ) ] )
ax0.set_ylim( [ 0.3, 1.7 ] )
ax0.set_xticks( () )
ax0.set_yticks( () )
ax0.set_aspect( "equal" )
pclabels = [ "Mutated", "Polymorphic", "Consensus" ]
scalevalues = [ 5, 20, 100, 500, 10000, 50000 ]	
if False:
	ccstep = len( rskeys ) / float( len( pcolors ) + len( scalevalues ) )
	for cc in range( len( pcolors ) + len( scalevalues ) ):
		ax1.pie( [ 1 ], colors = [ pcolors[cc] if cc < len( pcolors ) else pcolors[0] ], radius = 0.4 if cc < len( pcolors ) else size_func( scalevalues[ cc - len( pcolors ) ] ) , center = ( 1, cc * ccstep ) )
		ax1.text( 1, cc * ccstep + 0.7, pclabels[cc] if cc < len( pcolors ) else str( scalevalues[ cc - len( pcolors ) ] ), ha="left", size=24 )
for rcc in range( len( pcolors ) ):
	cc = len( pcolors ) - rcc - 1
	ax1.pie( [ 1 ], colors = [ pcolors[rcc] ], radius = 0.1, center = ( 1, cc ) )
	ax1.text( 1, cc + 0.15, pclabels[rcc], ha="center", size=24 )
ax1.set_ylim( [ -1, 3 ] )
ax1.set_xlim( [ 0.7, 1.3 ] )
ax1.set_xticks( () )
ax1.set_yticks( () )
ax1.set_aspect( "equal" )



	


plt.savefig( "poly_pie.png" )
plt.savefig( "poly_pie.eps" )



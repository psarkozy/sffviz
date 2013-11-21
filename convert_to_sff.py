#!/usr/bin/python2.6
import cProfile
import time
import pysam
from Bio import SeqIO
import optparse
import numpy as np
import scipy.misc.pilutil as smp
import sys
import copy

from numpy import *
from matplotlib.colors import LogNorm
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
#from pylab import plot,savefig
import matplotlib.pyplot as plt
from pylab import *
db=[]
try:
	dbfname='sampleid_barcode_refseq_mapping.tsv'
	db=[x.strip('\n').split('	') for x in open(dbfname).readlines()] #create our DB:
except:
	print 'cant open db', dbfname
	pass
#	116	116	117	117
# Barcode	BRCA1	BRCA2	BRCA1	BRCA2
# 13			NGS107	NGS90
# 14	NGS109			NGS92
# 25	NGS111	NGS94		


t=0
samba= 'TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGA'
class Myrec:
	def __init__(self):
		self.id = None
		self.seq = None
		self.letter_annotations = {}
		self.annotations = {}

def print_flowgram(record):
	members=['flow_values',	'flow_index',	'flow_chars',	'clip_adapter_right', 'clip_qual_right', 'clip_qual_left', 'clip_adapter_left',	'flow_key']
	print 'id	',record.id
	for m in members:
		try:
			
			print  m+'	',record.annotations[m]
		except:
			print 'member not in record',m, record
			raise
	print 'len(seq)	',len(record.seq)
	print 'seq	',record.seq
	print 'phred_quality type',type(record.letter_annotations['phred_quality'])
	print 'phred_quality len',len(record.letter_annotations['phred_quality'])
	print 'phred_quality',record.letter_annotations['phred_quality']


def print_alignment(a):
	print a.qname,'	',#query name, in this case, the read ID
	print a.qstart,'	',# start index of the aligned query portion of the sequence (0-based, inclusive)
	print a.qend,'	',# end index of the aligned query portion of the sequence (0-based, exclusive)
	print a.aend,'	',# aligned reference position of the read on the reference genome.  aend points to one past the last aligned residue. Returns None if not available.
	print a.alen,'	',# aligned length of the read on the reference genome. Returns None if not available.
	print a.cigar,'	',# CIGAR string, 
	print a.flag,'	',# flag, 16 is reversed read!
	print a.mapq,'	',# mapping quality
	print a.pos,'	',# 0-based leftmost coordinate
	print a.positions, '	',# a list of reference positions that this read aligns to.
	print a.qual,'	',#  read sequence base qualities, including soft clipped bases (None if not present).
	print a.query,'	',# aligned portion of the read and excludes any flanking bases that were soft clipped (None if not present). SAM/BAM files may included extra flanking bases sequences that were not part of the alignment. These bases may be the result of the Smith-Waterman or other algorithms, which may not require alignments that begin at the first residue or end at the last. In addition, extra sequencing adapters, multiplex identifiers, and low-quality bases that were not considered for alignment may have been retained.
	print a.rlen,'	', # length of the read (read only). Returns 0 if not given.
	print a.tlen,'	', # the insert size
	print a.seq,'	', # read sequence bases, including soft clipped bases (None if not present).
	print a.tags

	

if sys.argv[1]=='bam2sff' and len(sys.argv)>2:
	barcodes={}
	for line in open('bar-code.txt').readlines():
		bc=line.strip().split('	')
		barcodes[bc[0]]=bc[1]
	
	print sys.argv
	fname=sys.argv[2]
	for bc in barcodes.iterkeys():
		barseq=barcodes[bc]
		print barseq
		records=[]
		mainrecord=0
		for r in SeqIO.parse('singleread.sff', "sff"):
			mainrecord=r
		if '.bam' in fname or '.sam' in fname:
			print fname
			bamf=pysam.Samfile(fname,'rb')
			print 'loaded',fname
			cnt=0
			first=True
			newsff=open(fname.partition('.')[0]+'_'+bc+'_newsff.sff','wb')
			fastq=open(fname.partition('.')[0]+'_'+bc+'_newsff.fastq','w')
			fw=SeqIO.SffIO.SffWriter(newsff)
			#fw.write_file([mainrecord])
			for read in bamf.fetch(until_eof=True):
				if len(read.query)<9  or read.query[0:5] != barseq:
					#print barseq,read.query[4:9] 
					continue
				#print read.query
				cnt+=1
				# print_alignment(read)
				# print_flowgram(mainrecord)
				#print mainrecord
				# newrec=copy.copy(mainrecord)
				newrec=Myrec()
				newrec.id=read.qname#.replace('MB1YX','P2EWR')
				newrec.seq='TCAG'+read.query	
				newrec.annotations['flow_chars']=samba
				
				tags=read.tags
				# print tags
				pos=0
				seq='TCAG'+read.query
				fz=False
				for tag in tags: 
					if tag[0]=='FZ':
						newflowgram=[]
						newflowindex=[]#[1,2,3,2] #cause of samba and the added TCAG
						fz=tag[1]
						lastindex=-1
						seqindex=0
				if fz:
					for k in range(len(samba)):
						flowchar=samba[k]
						mer=0
						
						if seqindex < len(seq):
							
							if seq[seqindex]==flowchar: #we have at least 1 base at this flow
								mer=1
								# print k, lastindex,flowchar
								newflowindex.append(k-lastindex)
								lastindex=k
								while(seqindex+mer<len(seq) and seq[seqindex+mer]==flowchar):
									mer+=1
									newflowindex.append(0)
								seqindex+=mer
								#print mer
								newflowgram.append(100*(mer-1)+50+fz[k])
							else:
								newflowgram.append(max(0,fz[k]-50))
						else:
							newflowgram.append(0)
						# print flowchar,k,mer
					# for flowv in newflowgram: #ZM
						# newv=int((flowv*100)/256)
						# newv=flowv
						# if newv<1000 and newv>0:
							# flowvspos[newv][pos]+=1
						# pos+=1
				else:
					print 'NO FZ tag in read',read.qname
				newrec.annotations['flow_values']=newflowgram
				newrec.annotations['flow_index']=newflowindex
				newqual=[32,32,32,32] #for our TCAG key
				for c in read.qual:
					newqual.append(max(0,ord(c)-33))
				newrec.letter_annotations['phred_quality']=newqual
				newrec.annotations['clip_adapter_left']=0
				newrec.annotations['clip_adapter_right']=0
				newrec.annotations['clip_qual_left']=13
				newrec.annotations['clip_qual_right']=len(newrec.seq)-4
				newrec.annotations['flow_key']='TCAG'
				fastq.write('@'+newrec.id+'\n')
				fastq.write(newrec.seq[newrec.annotations['clip_qual_left']:newrec.annotations['clip_qual_right']]+'\n+\n')
				fastq.write(read.qual[newrec.annotations['clip_qual_left']:newrec.annotations['clip_qual_right']])
				fastq.write('\n')
				
				
				records.append(newrec)
				# print_flowgram(newrec)
				# if first:
					# first=False
					# fw.write_file([newrec])
				# else:
					# fw.write_record(newrec)
				if cnt%1000==0:
					print 'wrote ',cnt,'reads'
				if cnt>600000:
					break
			
			fw.write_file(records)
	
print sys.argv
#!/usr/bin/python2.6
import pysam
import math
import sys
from numpy import zeros
from numpy import *
import time
from matplotlib.colors import LogNorm

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

#from pylab import plot,savefig
import matplotlib.pyplot as plt
from pylab import *
from Bio import SeqIO

#-------PARAMS-----
maxflowsize=1000
maxreads=3000000
#------------------

print 'lauching'


histogram=[[],[],[],[]]  #tacg , RGBK
badextensions=[]
flowsum=[]
flowcnt=[]
phreds=ones((64,maxflowsize))
for i in range(maxflowsize):

	flowsum.append(0)
	flowcnt.append(1)
	badextensions.append(0)
	for j in range(4):
		histogram[j].append(0)
phreds_hist=[]
flows_hist=[]
records=[]
#records= SeqIO.parse("ecoli_10x_0.sff", "sff")
fname=sys.argv[1]
clip=False
makephredvsflow=False
makebadextensions=False
makeflowhist=False
makeflowpos=False
makeflowmerpos=False
makemerflow=False
bases=['A','C','G','T']
flowchars='' 
flowcharmap={}
def makeFlowCharMap(chars):
	global flowcharmap
	baseindex={'A':0,'C':0,'G':0,'T':0}
	pos=0
	for char in chars:
		flowcharmap[char][pos]=baseindex[char]
		baseindex[char]+=1
		pos+=1
	for base in bases:
		flowvspos[base]=ones((maxflowsize,len(flowcharmap[char])))
flowvspos={}
for base in bases:
	flowvspos[base]=[[],[]]
	flowcharmap[base]={}
flowmerpos=[]
maxmers=20
def makeFlowMerPos(flowseq):
	global flowmerpos
	for mer in range(maxmers):
		flowmerpos.append(ones((maxflowsize,len(flowseq))))
		
numflows=100
if len(sys.argv)>2 and sys.argv[2]=='clip':
	clip=True
	print 'Clipping on'
	
if len(sys.argv)>2 and sys.argv[2]=='phredflow':
	makephredvsflow=True
	print 'makephredvsflow on'
	
if len(sys.argv)>2 and sys.argv[2]=='badextensions':
	makebadextensions=True
	print 'makebadextensions on'
if len(sys.argv)>2 and sys.argv[2]=='flowhist':
	makeflowhist=True
	print 'makeflowhist on'
if len(sys.argv)>2 and sys.argv[2]=='flowvspos':
	makeflowpos=True
	print 'makeflowhist on'
if len(sys.argv)>2 and sys.argv[2]=='flowmerpos':
	makeflowmerpos=True
	print 'makeflowhist on'
def supersumer(flow_indices):#read.annotations['flow_index'])
	ss=[flow_indices[0]]
	for i in xrange(1,len(flow_indices)):
		ss.append(ss[-1]+flow_indices[i])
	return ss
	
i=0
tards=0
t=time.clock()
lengths=[]
for record in SeqIO.parse(fname, "sff"):
	i+=1
	lengths.append(len(record.seq))
	if makeflowpos:
		if i<10:
			numflows=len(record.annotations['flow_chars'])
			flowchars=record.annotations['flow_chars']
			makeFlowCharMap(flowchars)
			#	print flowcharmap
		k=0
		for flowvalue in record.annotations['flow_values']:
			
			if flowvalue<maxflowsize and flowvalue:
				flowchar=flowchars[k]
				try:
					x=flowcharmap[flowchar][k]
					flowvspos[flowchar][flowvalue][x]+=1
					# flowvspos[flowchar][0].append()
					# flowvspos[flowchar][1].append(flow)
				except KeyError:
					print flowchars
					print flowchar,k,flowcharmap,flowvalue
					raise
			k+=1
		if i>maxreads:
			break
	if makeflowmerpos:
		if i==1:
			flowchars=record.annotations['flow_chars']
			numflows=len(record.annotations['flow_chars'])
			makeFlowMerPos(flowchars)
			print flowchars
		k=0
		supersum=[]
		fli=record.annotations['flow_index']
		supersum=supersumer(fli)#[sum(fli[0:x]) for x in range(len(fli))][1:]
		#print supersum
		#print supersumer(fli)
		#exit(1)
		fii=0
		#print supersum
		while (fii<len(supersum)):
			fi=supersum[fii]-1
			mer=1
			while(len(supersum)>fii+1 and supersum[fii]==supersum[fii+1]):
				mer+=1
				fii+=1
			flowvalue=record.annotations['flow_values'][int(fi)]
		#	print mer,fi,int(fi),fii,flowvalue
			#print record.annotations['flow_index']
			#print flowmerpos
			if flowvalue<maxflowsize:
				flowmerpos[mer][flowvalue][fi]+=1
			fii+=1
			
		if i>maxreads:
			break
	if makephredvsflow:
		k=0 
		flowindex=-1
		#flow_values=record.annotations['flow_values']
		for phred in record.letter_annotations['phred_quality']:
			flowindex+=record.annotations['flow_index'][k]
			#if record.annotations['flow_chars'][flowindex]==record.seq[k].upper():
			phreds[phred][min(maxflowsize-1,record.annotations['flow_values'][flowindex])]+=1
				#phreds_hist.append(phred)
				#flows_hist.append(min(maxflowsize-1,record.annotations['flow_values'][flowpos]))
			#else:
				#print 'non-matching bases in seq vs flowchars:',record.seq[k],record.annotations['flow_chars'][flowpos-1]
					
			k+=1 
			
	if makeflowhist:
		rclip=record.annotations['clip_qual_right']
		lclip=record.annotations['clip_qual_left']
		if rclip==0:
			rclip=len(record.seq)
		if clip:
			rclip=rclip/2 #OOOOOOOOO
		j=0
		flowindex=-1
		for flow in record.annotations['flow_values']:	
			flowindex+=record.annotations['flow_index'][j]
			j+=1	
			if clip and rclip<=flowindex:
				continue
			if j<maxflowsize:
				flowsum[j]+=flow		
				flowcnt[j]+=1		
			if flow>=maxflowsize:
				print 'flow value greater than ',maxflowsize,'!',flow,record.id
				continue
		
			histogram[j%4][flow]+=1
		fi=-1
	if makebadextensions:
		rclip=record.annotations['clip_qual_right']
		lclip=record.annotations['clip_qual_left']
		if rclip==0:
			rclip=len(record.seq)
		if clip:
			rclip=rclip/2 #OOOOOOOOO
		for k in range(len(record.annotations['flow_index'])): 
			if clip and (k>rclip or k<lclip):
				continue
				p=0
			index=record.annotations['flow_index'][k]
			
			if index>3 and fi>=0:
				if fi==-1:
					print record.seq,record.annotations['flow_key']
				badextensions[fi]+=1
				tards+=1
				# print 'tards in ',fi,'[',
				# for f in xrange(fi,fi+index+1):
					# print record.annotationpyps['flow_values'][f],
				# print ']'
			fi+=index	#print ',',index,
	if i%1000==0:
		print i, 'clipping at','total length=',len(record.annotations['flow_index']),'time=', time.clock()-t
		t=time.clock()
		
	if i>maxreads:
		break
	#records=[record]
print tards,'tards, in',i,'reads'
	#print record.id, len(record), record.seq[:20]+"..."
	#222018 tards in 100k reads
numreads=i
print 'read all'
# print records[0]
# print records[0].annotations['flow_values']
#make incomplete extension plot:
if makebadextensions:
	print 'badextensions[999]',badextensions[999]
	badextensions[999]=0
	plt.plot([a for a in range(1000)],badextensions[:1000],'r')
	plt.title('%i bad extensions vs flow number'%(tards)) 
	if clip:
		plt.savefig(fname+'_badextensions_clipped.png')#,dpi=200)
	else:
		plt.savefig(fname+'_badextensions.png')#,dpi=200)

if makeflowpos:
	for base in bases:
		print 'makeflowpos',base
		plt.clf()
		plt.subplot(111, axisbg='#000080')
		# plt.hist2d(flowvspos[base][0],flowvspos[base][1],bins=[len(flowcharmap[base]),1000],norm=LogNorm())
		plt.imshow(flowvspos[base],norm=LogNorm(), interpolation='nearest', origin='lower',aspect=0.13)
		plt.title(fname+' '+base+' base flow index vs flow value')
		plt.grid(b=True,axis='y', which='both',markevery=100)
		colorbar()
		plt.savefig(fname+'_'+base+'_imshow_flowvspos.png',dpi=200)
	outf=open(fname+'_flowvspos.tsv','w')
	for base in flowchars:
		outf.write(base+'	')
	outf.write('\n')
	fullimage=ones((maxflowsize,len(flowchars)))
		
	for flowv in range(maxflowsize):
		outf.write(str(flowv))
		i=0
		for char in flowchars:
			
			x=flowcharmap[char][i]
			outf.write('	'+str(flowvspos[char][flowv][x]))
			fullimage[flowv][i]+=flowvspos[char][flowv][x]
			i+=1
		outf.write('\n')
	outf.close()
	plt.clf()
	plt.subplot(111, axisbg='#000080')
	# plt.hist2d(flowvspos[base][0],flowvspos[base][1],bins=[len(flowcharmap[base]),1000],norm=LogNorm())
	plt.imshow(fullimage,norm=LogNorm(), interpolation='nearest', origin='lower',aspect=float(len(flowchars)/float(maxflowsize)))
	plt.title(fname+' FULL base flow index vs flow value')
	plt.grid(b=True,axis='y', which='both',markevery=100)
	colorbar()
	plt.savefig(fname+'_FULL_imshow_flowvspos.png',dpi=200)

if makeflowmerpos:
	revlook=[]
	for i in range(len(flowchars)):
		found=False
		for k in range(i-1,max(0, i-8),-1):
			if flowchars[k] == flowchars[i]:
				revlook.append(i-k)
				found=True
				break
		if not found:
			revlook.append(4)
	print revlook
	curvedb={}
	colors=['','','m','r','y','c','g','b']
	for base in bases:
		plt.clf()
		for mer in range(1,4):
			for rv in range(2,8):
				xs=[]
				ys=[]
				for fi in range(350):#len(flowchars)):
					if flowchars[fi]!=base:
						continue
					if revlook[fi]==rv:
						fstart=max(0,mer*100-49)
						fend=mer*100+50
						avg=0
						cnt=0
						for r in range(fstart,fend):
							cnt+=flowmerpos[mer][r][fi]-1
							avg+=(flowmerpos[mer][r][fi]-1)*r
						
						avg=avg/cnt
						if cnt>0:
							xs.append(fi)
							ys.append(avg)
						#print avg, cnt
						# flowmercnt=sum(flowmerpos[mer][fstart:fend][fi])
						# print len(range(fstart,fend)), len(flowmerpos[mer][fstart:fend][fi]), range(fstart,fend), flowmerpos[mer][fstart:fend][fi]
						# ys.append( dot(range(fstart,fend), flowmerpos[mer][fstart:fend][fi])/flowmercnt)
						
				
				plt.plot(xs,ys, color=colors[rv])
				plt.plot(xs,ys, colors[rv]+'|')
				print base,mer,rv
		plt.title(fname+'base'+base+' base flow index vs flow value')
		plt.savefig(fname+'_'+str(base)+'_curves_dotted_flowmerpos_v2.png',dpi=200)

	for mer in range(maxmers):
		print 'makeflowpos',mer
		plt.clf()
		plt.subplot(111, axisbg='#000080')
		# plt.hist2d(flowvspos[base][0],flowvspos[base][1],bins=[len(flowcharmap[base]),1000],norm=LogNorm())
		plt.imshow(flowmerpos[mer],norm=LogNorm(), interpolation='nearest', origin='lower',aspect=0.5)
		plt.title(fname+' MER:'+str(mer)+' base flow index vs flow value')
		plt.grid(b=True,axis='y', which='both',markevery=100)
		colorbar()
		plt.savefig(fname+'_'+str(mer)+'_imshow_flowmerpos.png',dpi=200)
	# outf=open(fname+'_flowmerpos.tsv','w')
	# for mer in range(maxmers):
		# outf.write(str(mer)+'	')
	# outf.write('\n')
	# for flowv in range(maxflowsize):
		# outf.write(str(flowv))
		# i=0
		# for char in flowchars:
			# x=flowcharmap[char][i]
			# outf.write('	'+str(flowvspos[char][flowv][x]))
			# i+=1
		# outf.write('\n')
	# outf.close()
	
	
if makeflowhist:
	for i in range(maxflowsize):
		flowsum[i]/=float(flowcnt[i])
		for j in range(4):
			histogram[j][i]=math.log(float(histogram[j][i]+1))
	x=[a/100.0 for a in range(maxflowsize)]
	plt.clf()
	plt.grid(b=True,which='major', linestyle='-',color='k')
	plt.plot(x, histogram[0],'r')
	plt.plot(x, histogram[1],'g')
	plt.plot(x, histogram[2],'b')
	plt.plot(x,histogram[3],'k')

	## print histogram[0]
	if clip:
		plt.savefig(fname+'_histogramclipped_per2.png',dpi=200)
	else:
		plt.savefig(fname+'_histogram.png',dpi=200)
		
	plt.clf()
	plt.plot(range(maxflowsize),flowsum,'r')
	if clip:
		plt.savefig(fname+'_flowvalues_clipped.png',dpi=200)
	else:
		plt.savefig(fname+'_flowvalues.png',dpi=200)

plt.clf()
#plt.hexbin(phreds_hist,flows_hist,gridsize=600)
# hist2d(phreds_hist,flows_hist, bins=[64,maxflowsize])
## print histogram[0]
#plt.savefig(fname+'_flow_vs_phred_hist.png',dpi=200)
plt.clf()
if makephredvsflow:
	plt.clf()
	for x in range(64):
		for y in range( maxflowsize):
			phreds[x][y]=log(phreds[x][y])

	plt.imshow(phreds,aspect=4,origin='lower')
	# hist2d(phreds_hist,flows_hist, bins=[64,maxflowsize])
	## print histogram[0]

	plt.savefig(fname+'_flow_vs_phred_hist_imshow.png',dpi=200)
	plt.clf()

print 'Total read cnt=',len(lengths)
avglen=sum(lengths)/float(len(lengths))
variances=[(l-avglen)*(l-avglen) for l in lengths]
stdev=math.sqrt(sum(variances)/float(len(lengths)))
print 'avglen=',avglen
print 'stdev=',stdev

# ID: gi|215263233|emb|FM180568.1|
# Name: gi|215263233|emb|FM180568.1|
# Number of features: 0
#/flow_values=(91, 8, 77, 8, 9, 100, 9, 113, 8, 103, 10, 9, 106, 9, 81, 7, 100, 97, 97, 106, 7, 111, 7, 7, 90, 8, 10, 83, 99, 7, 8, 111, 8, 9, 94, 99, 8, 8, 185, 11, 14, 96, 8, 101, 205, 8, 101, 8, 86, 8, 7, 115, 211, 101, 204, 100, 8, 7, 93, 6, 5, 207, 195, 9, 200, 7, 9, 88, 8, 214, 8, 106, 9, 91, 98, 11, 196, 6, 99, 113, 193, 10, 8, 97, 7, 12, 310, 94, 6, 10, 190, 9, 11, 104, 7, 9, 109, 6, 92, 80, 86, 6, 82, 8, 8, 93, 7, 108, 7, 109, 12, 7, 93, 5, 5, 108, 6, 204, 12, 10, 191, 6, 12, 106, 7, 198, 13, 9, 192, 8, 113, 6, 122, 9, 104, 5, 106, 8, 8, 205, 7, 113, 8, 6, 186, 9, 101, 8, 8, 99, 102, 8, 5, 101, 119, 212, 99, 6, 99, 11, 6, 99, 10, 104, 9, 8, 99, 5, 7, 89, 9, 11, 202, 99, 10, 204, 9, 10, 104, 13, 8, 120, 6, 104, 9, 8, 90, 77, 114, 8, 87, 8, 9, 93, 5, 101, 8, 9, 98, 6, 11, 225, 101, 7, 113, 7, 100, 10, 9, 99, 211, 6, 8, 83, 10, 11, 110, 101, 7, 111, 10, 15, 80, 8, 106, 98, 205, 90, 104, 7, 98, 189, 7, 8, 93, 10, 5, 87, 9, 9, 107, 7, 7, 97, 12, 9, 106, 8, 8, 108, 10, 100, 7, 94, 11, 11, 111, 7, 8, 279, 93, 5, 7, 209, 7, 90, 8, 90, 10, 107, 6, 10, 102, 8, 7, 98, 4, 102, 8, 9, 216, 9, 300, 209, 7, 7, 104, 9, 97, 8, 88, 10, 188, 9, 107, 6, 7, 125, 9, 113, 8, 94, 7, 168, 6, 8, 98, 10, 9, 115, 11, 90, 9, 5, 183, 11, 14, 107, 9, 90, 7, 11, 127, 8, 8, 85, 189, 94, 8, 8, 115, 10, 8, 109, 7, 90, 7, 7, 120, 93, 6, 7, 291, 111, 10, 6, 99, 123, 8, 10, 206, 8, 7, 95, 6, 9, 124, 10, 101, 95, 7, 11, 213, 7, 9, 418, 12, 10, 102, 222, 8, 6, 84, 204, 9, 8, 109, 94, 9, 127, 114, 7, 274, 8, 7, 114, 93, 111, 7, 99, 6, 8, 97, 8, 109, 212, 214, 5, 7, 216, 8, 89, 8, 10, 96, 7, 92, 92, 6, 12, 114, 8, 11, 111, 7, 203, 105, 9, 7, 86, 5, 5, 108, 66, 9, 8, 95, 88, 9, 11, 226, 12, 7, 110, 6, 104, 87, 99, 82, 9, 14, 97, 7, 8, 86, 7, 101, 234, 11, 213, 222, 114, 8, 73, 100, 134, 9, 124, 212, 6, 8, 78, 8, 110, 11, 90, 11, 108, 6, 8, 82, 12, 8, 79, 7, 6, 99, 113, 117, 17, 231, 10, 111, 6, 6, 111, 219, 107, 7, 7, 358, 9, 11, 76, 9, 98, 6, 97, 6, 7, 286, 8, 112, 12, 10, 111, 84, 90, 95, 103, 7, 10, 210, 6, 7, 98, 4, 89, 15, 5, 80, 5, 119, 206, 7, 11, 85, 8, 97, 105, 97, 98, 15, 6, 111, 7, 182, 237, 8, 7, 94, 9, 13, 84, 6, 8, 134, 9, 8, 94, 7, 9, 129, 98, 7, 212, 89, 8, 9, 90, 10, 295, 231, 10, 100, 9, 10, 104, 7, 5, 104, 108, 7, 8, 109, 5, 11, 97, 6, 191, 5, 9, 122, 9, 94, 6, 10, 110, 12, 8, 116, 9, 92, 11, 6, 110, 13, 8, 92, 8, 9, 212, 84, 108, 114, 13, 81, 9, 10, 300, 6, 8, 108, 7, 9, 89, 9, 9, 206, 296, 9, 9, 111, 8, 6, 208, 193, 12, 293, 207, 10, 9, 70, 8, 6, 119, 7, 5, 308, 5, 18, 115, 7, 134, 5, 9, 110, 9, 8, 112, 8, 9, 122, 9, 11, 116, 13, 13, 126, 9, 255, 6, 11, 117, 15, 176, 9, 13, 112, 8, 110, 9, 89, 11, 6, 103, 5, 114, 94, 7, 7, 211, 8, 113, 8, 12, 116, 8, 107, 11, 9, 84, 4, 9, 103, 93, 4, 8, 97, 8, 128, 12, 86, 306, 17, 356, 13, 77, 6, 11, 116, 71, 83, 11, 9, 113, 8, 13, 72, 124, 92, 93, 7, 83, 11, 9, 121, 7, 119, 8, 8, 95, 15, 206, 6, 111, 126, 14, 97, 96, 180, 5, 6, 88, 5, 7, 187, 6, 16, 92, 93, 12, 4, 208, 86, 105, 8, 90, 10, 9, 109, 10, 97, 4, 7, 93, 5, 270, 8, 5, 190, 5, 8, 190, 6, 6, 95, 9, 97, 12, 87, 6, 11, 178, 87, 7, 90, 12, 325, 13, 97, 4, 87, 7, 118, 8, 8, 239, 6, 7)
#/flow_index=(1, 2, 3, 2, 2, 3, 2, 2, 1, 1, 1, 2, 3, 3, 1, 3, 3, 1, 3, 0, 3, 2, 1, 0, 2, 2, 3, 1, 0, 1, 1, 0, 1, 3, 3, 0, 1, 0, 2, 0, 3, 2, 0, 2, 2, 1, 2, 0, 2, 1, 1, 0, 3, 3, 0, 0, 1, 3, 0, 3, 3, 2, 1, 1, 2, 3, 2, 2, 3, 3, 2, 0, 3, 0, 3, 2, 0, 3, 0, 2, 2, 2, 2, 3, 0, 2, 3, 0, 2, 3, 1, 3, 1, 1, 0, 1, 2, 3, 2, 3, 3, 3, 0, 1, 2, 0, 3, 3, 2, 3, 1, 1, 2, 3, 2, 3, 3, 0, 1, 2, 2, 3, 1, 0, 3, 3, 1, 2, 3, 2, 1, 1, 0, 1, 1, 2, 1, 0, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 0, 0, 1, 3, 0, 2, 2, 2, 3, 3, 2, 3, 0, 2, 0, 0, 1, 0, 3, 2, 2, 2, 0, 2, 3, 2, 2, 2, 0, 3, 3, 2, 3, 0, 3, 2, 3, 3, 1, 0, 1, 3, 3, 2, 3, 1, 3, 0, 0, 1, 3, 1, 3, 0, 3, 3, 2, 1, 3, 0, 3, 0, 0, 0, 3, 1, 0, 3, 1, 0, 3, 1, 2, 1, 2, 0, 0, 3, 1, 1, 2, 3, 2, 1, 0, 1, 0, 3, 0, 2, 3, 2, 1, 3, 3, 2, 0, 1, 3, 3, 1, 3, 1, 3, 0, 3, 2, 1, 1, 1, 3, 3, 2, 1, 0, 2, 0, 1, 0, 1, 2, 1, 1, 2, 1, 0, 3, 2, 2, 2, 3, 3, 3, 1, 1, 2, 0, 2, 3, 1, 0, 1, 3, 0, 0, 0, 3, 2, 2, 3, 0, 0, 2, 3, 1, 1, 1, 1, 3, 0, 3, 2, 3, 2, 1, 0, 3, 2, 1, 1, 1, 3, 2, 0, 1, 0, 3, 3, 3, 3, 3, 1, 2, 0, 1, 3, 2, 0, 0, 1, 0, 2, 3, 3, 1, 3, 3, 2, 0, 3, 2, 3, 3, 2, 3, 3, 3, 0, 1, 1, 1, 2, 3, 0, 0, 3, 3, 3, 0, 1, 0, 0, 3, 3, 0, 1, 0, 2, 0, 0, 1, 0, 3, 3, 3, 0, 0, 3, 2, 3, 3, 3, 3, 3, 2, 0, 0, 3, 2, 0, 3, 2, 2, 3, 2, 1, 3, 0, 2, 3, 2, 3, 3, 1, 3, 2, 2, 1, 0, 0, 2, 0, 0, 0, 2, 3, 1, 1, 3, 3, 1, 1, 1, 2, 3, 2, 3, 2, 0, 2, 1, 2, 1, 1, 0, 3, 3, 0, 3, 1, 3, 0, 1, 1, 2, 3, 2, 3, 2, 0, 0, 3, 0, 3, 0, 3, 2, 2, 3, 0, 1, 2, 2, 0, 0, 2, 2, 2, 3, 0)
# /flow_chars=TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
# /clip_adapter_right=0
# /clip_qual_right=478
# /clip_qual_left=4
# /clip_adapter_left=0
# /flow_key=TCAG
# Per letter annotation for: phred_quality
# Seq('tcagATCTACGATGTGCGCCAGTTCTGTTACCGCAACCTTGAAGACTTCGTTGC...CAA', DNAAlphabet())


# gi|215263233|emb|FM180568.1|       16      0       0       100     [(0, 6), (2, 1), (0, 73)]       -1      -1      79      AGCTTTCATTCTGCCTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTG ((]]OO]]]]]W]]MM]===]]]]]]]]]]]]]]]]SS]WW)))))))]]]]]]]]]]]]]]]]]]]]]][[]]]]]]]      [('NM', 10)]
# next gi|215263233|emb|FM180568.1|       16      0       0       100     [(0, 86)]       -1      -1      86      AGCTTTTCATTCTGCCTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACC  :::]]]]]]]WW]]]FF]QQQ]AA]]]]]]]]]]]]]]]]SS)))))))]]]]]]]]]]]]]]]]]]]]]]]]]]]]GG]BBB]]]       [('NM', 10)]
# next gi|215263233|emb|FM180568.1|       16      0       0       100     [(0, 93)]       -1      -1      93      AGCTTTTCATTCTGCCTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTG   ---]]]]]]]KK]]]SS]+++]YY]]]]]]]]]]]]]VVYEE-------]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]MM]]WW]]]]]]]        [('NM', 10)]

# 1 QNAME String [!-?A-~]f1,255g Query template NAME
# 2 FLAG Int [0,216-1] bitwise FLAG
# 3 RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
# 4 POS Int [0,229-1] 1-based leftmost mapping POSition
# 5 MAPQ Int [0,28-1] MAPping Quality
# 6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
# 7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next segment
# 8 PNEXT Int [0,229-1] Position of the mate/next segment
# 9 TLEN Int [-229+1,229-1] observed Template LENgth
# 10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
# 11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
# QNAME	FLAG	RNAME	POS	MAPQ	CIGAR	RNEXT	PNEXT	TLEN	SEQ	QUAL	
# next1k: ECOLI14JDOYRCMZ6        0       0       8485    100     [(0, 44), (1, 1), (0, 61), (2, 1), (0, 121), (1, 1), (0, 9), (1, 1), (0, 60), (1, 1), (0, 150), (1, 1), (0, 50), (2, 1), (0, 19), (5, 1)]       -1      -1      519     TCTACAACGATGCAGGTATTAGCAACGATCGTATTCTGATCAAAACTGGCTTCTACCTGGCAGGGTATCCGTGCTGCAGAACAGCTGGAAAAAGAAGGCATCAACTTAACCTGACCCTGCTGTTCTCCTTCGCTCAGGCTCGTGCTTGTGCGGAAGCGGGCGTGTTCCTGATCTCGCCGTTTGTTGGCCGTATTCTTGACTGGTACAAAGCGAATACCGATAAGAAAAGAGTACGCTACCGGCAGAAGATCCGGGCGTGGTTTCTGTATCTGAAATCTACCAGTACTACAAAGAGCACCGGTTATGAAACCGTGGTTATGGGCGCAAGCTTCCGTAACATCGGCGAAATTCTGGAACTGGCAGGTTGCGACCGTCTGACCATCGCACCGGCACTGCTGAAAGAGCTGGCTGAGAGNGAAGGGGCTATCGAACGTAAACTGTCTTACACCCGGCGAAGTGAAAGCGCGTCCGGCGCGTATCACTGAGTCCGAGTTCCTGTGCAGGATAACCAGGATCCAA        ]]]]]]]]]]]]]]]]]]KK]]TT]]]PP]]]]]]]]]]]]]]]]????]]ZZ]99]]]XX]]]]]PPP]]]PP]]]]]]\]W]]]]]]]]]-----]66DD]]X]RR];;]]AA]]](((]]]]]VV]]??[[]]]]]]KK]]]]]]]QQ]YY]XXMM]]EEEUT]]KKUU]]]]][]YMM]BBB];;SS;;]JXDDT??W]]A,,U]]@@@]]]JJXQKKQ\]CC\++++]]]]]]]]]]//UU]]]OO]]]NN===]]]AA>>>]\]XM]XX[>>>W]]]@@]]]]A]VI===]E]]]::CC>>>W]555QQ]]GGCC]](((]]];;J]HHFF]]JJ]]YU//U]...77]TDDAAIWLL]]==BB[]]]::P]\]RVGGS]]ROQHHHHK]]]]VST,,,ETD]P;;CD[K8Z]!J55))))YV6]TG77GQ]///PWPVUEE#]P$$$44GF;;]LI)))SPGL[NBB;;RIMLQHXKHCDOWGZCCVZG55??.EH/EL??CT;;??<55GG   [('NM', 51)]


#-------------------------------------------------------------------------------------
#!/usr/bin/env python
# Usage: python histogram_from_bed.py bed histogramfigure nbinaries
#This program obtains a scores histogram for the  .bed file and saves this figure to the histogramfigure file name.
#-------------------------------FUNCTIONS-LIST----------------------------------------
#def process_inputfile(file1,fout):
#def int_main(file1,fout):
#---------------------------------FUNCTIONS-------------------------------------------



def createhistogram(scores,figurename,nbin,file1):
	import numpy
	import pylab
#number of binaries
#	nbin=20
	n,bins,patches=pylab.hist(scores,bins=nbin)
	pylab.title('scores histogram for'+file1)
	pylab.xlabel('background scores')
	pylab.savefig(figurename+'.png')



def process_inputfile(file1,figurename,nbin):
	g = open(file1,'r')
	lineas=g.readlines()
	nline=0;
	total_lineas=len(lineas)
	scores=[]
	while nline<total_lineas :
		line=lineas[nline]
		chrom, startp, endp,name,score,strand = line.split("\t")
		scores.append(int(score))
		nline=nline+1;
	g.close()
	createhistogram(scores,figurename,nbin,file1)



def int_main(file1,figurename,nbin):
	print "Input file - "+file1
	print "figure Name -"+figurename
	print "Number of binaries -"+nbin
	process_inputfile(file1,figurename,int(nbin))



#----------------------------------MAIN-------------------------------------------

import sys
int_main(sys.argv[1],sys.argv[2],sys.argv[3])


#-------------------------------------------------------------------------------------
#!/usr/bin/env python
# Usage: python stats_from_bed.py bed mean
#This program obtains a scores histogram for the  .bed file and saves this figure to the histogramfigure file name.
#-------------------------------FUNCTIONS-LIST----------------------------------------


def process_inputfile(file1):
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
	return scores


def getstatistics(scores,mean,alpha):
	import scipy
	from scipy import stats	
#	confidence intervals
	ascores=scipy.array(scores)
	sem=scipy.stats.sem(ascores)
	samplemean=ascores.mean()
	if alpha==0.01 :
		factor=2.58
	elif alpha==0.05 :
		factor=1.96
	ci=[samplemean-factor*sem,samplemean+factor*sem]
# 	pvalue - ttest	
	[t,pvalue]=scipy.stats.ttest_1samp(ascores,mean)
	return ci,pvalue/2

	
# ttest - 1 sample





def int_main(file1,mean,alpha):
	print "Input file - "+file1
	print "mean"+mean
	alpha=float(alpha)
	scores=process_inputfile(file1)
	mean=float(mean)
	if len(scores) >30 :
		ci,pvalue=getstatistics(scores,mean,alpha)
		if pvalue>alpha :
			print "Cannot reject the null hipothesis (mean = "+str(mean)+"), pvalue("+str(pvalue)+")"
			print "Mean confidence intervals (alpha="+str(alpha)+"), ci:"+str(ci) 
		else:
			print "Reject the null hypothesis (mean = "+str(mean)+"), pvalue("+str(pvalue)+")"
			print "Mean confidence intervals (alpha="+str(alpha)+"), ci:"+str(ci) 
	else :
		print "Not enough data to make an inference"
	
	
# For statistica test -> float(mean)
	



#----------------------------------MAIN-------------------------------------------

import sys
int_main(sys.argv[1],sys.argv[2],sys.argv[3])

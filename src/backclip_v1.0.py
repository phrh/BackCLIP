# BackCLIP - identifies the presence of common background in a PAR-CLIP dataset.
#
# Created by Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra, Santiago Herrera on 5 February 2015
# Copyright (c) 2015 Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra, Santiago Herrera. All rights reserved.
#
#
# BackCLIP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2.
#-------------------------------------------------------------------------------------
#!/usr/bin/env python
# USAGE: python backclip.py filebackground fileclusters fileparameters
# Use this program to have a quantitative measure of the amount of background present in a determined dataset.
# INPUT:         fileparameters -> file with parameters such as:
#		 	alpha=significance level(default value 0.01)
#		 	threshold=to define if the amount of common background is significant (default (maximum score)/2) 
#		 	histogram=give as output histogram of the intersection scores (default false)
#			fbackground=(.BED)file background and corresponding scores
#			fclusters=(.BED)file with clusters detected from PAR-CLIP dataset	
#		 	filename=in case histogram is true, the name of the output file (default fileclusters.histbackground)	
#			namebed=name bed with set of common background found in fileclusters
# OUTPUT: confidence interval (significance level 99%)
# 	  histogram score distribution in the intersection file (saved in filename)
#	  common background in file clusters (saved as namebed)
#-------------------------------FUNCTIONS-LIST----------------------------------------
#	getparameters(parameters)
#	getstatistics(scores,alpha,threshold)
#	findseq(chrom,startp,endp,seq,strand,fclu)
#	findscore(chrom,startp,endp,seq,strand,fint2)
#	clean_intersection(fintfil,fclu,fint2)
#	getthreshold(fback)
#	createhistogram(scores_hist,figurename,nbin,file1)
#	DeleteDuplicatesList(x)
#	deleteandmerge(bedint)
#	check_overlap(lineas,chrom,startp,endp,score,strand,nline)
#	process_inputfile(parameters)
#	int_main(parameters) 
#---------------------------------FUNCTIONS-------------------------------------------

def getparameters(parameters):
	from datetime import date
	fparam=open(parameters,'r')
	lineas=fparam.readlines()	
#default parameters in case are deleted from parameters file
	alpha=0.01
	histogram=False
	smotif=False
	filename='histbackground '+str(date.today().ctime())
	namebed=filename+'intersection.bed'
	threshold=0
	for line in lineas:
		paramname,value=line.split("=")
		if paramname.lower().strip() == "alpha" :
			alpha=float(value)
		elif paramname.lower().strip() == "histogram" :
			if value.lower().strip() == "true" :
				histogram=True
		elif paramname.lower().strip() == "filename" :		
			filename=value.strip()
		elif paramname.lower().strip() == "namebed" :
			namebed=value.strip()+'.intersection.bed'
		elif paramname.lower().strip() == "fbackground" :
			fbackground=value.strip()
		elif paramname.lower().strip() == "fclusters" :
			fclusters=value.strip()
		elif paramname.lower().strip() == "threshold" :
			threshold=int(value.strip())
	return alpha,histogram,filename,namebed,fbackground,fclusters,threshold


def getstatistics(scores,alpha,threshold):
	import scipy
	from scipy import stats	
#	confidence intervals
	if len(scores) >= 30:
		ascores=scipy.array(scores)
		sem=scipy.stats.sem(ascores)
		samplemean=ascores.mean()
		if alpha==0.01 :
			factor=2.58
			ci=[samplemean-factor*sem,samplemean+factor*sem]
		elif alpha==0.05 :
			factor=1.96
			ci=[samplemean-factor*sem,samplemean+factor*sem]
		else :
			print "sample size must be >30"
			ci=[]
		[t,pvalue]=scipy.stats.ttest_1samp(ascores,threshold)
	return ci,pvalue/2


def findseq(chrom,startp,endp,seq,strand,fclu):
	flag=0
	for interval in fclu:
		if seq in interval.name :
			cchrom, cstartp, cendp,cseq,cscore,cstrand = interval.chrom,interval.start,interval.end,interval.name,interval.score,interval.strand
			if (cstartp<=startp) and (endp<=cendp) and (cchrom==chrom): 
				newseq=cseq[startp-cstartp:endp-cstartp+1]
				flag=1			
	return newseq


def findscore(chrom,startp,endp,seq,strand,fint2):
	newscore=0
	for interval in fint2 :
		if (startp==interval.start) and (endp==interval.end) and (chrom==interval.chrom): 
			newscore=interval.score
	return newscore


def clean_intersection(fintfil,fclu,fint2):
	import pybedtools
	items=[]
	scores=[]
	for interval in fintfil :
		chrom,startp,endp,seq,score,strand=interval.chrom,interval.start,interval.end,interval.name,interval.score,interval.strand
		newseq=findseq(chrom,startp,endp,seq,strand,fclu)
		newscore=findscore(chrom,startp,endp,seq,strand,fint2)
		items.append((chrom,startp,endp,newseq,newscore,strand))
		scores.append(int(newscore))
	newfintfil=pybedtools.BedTool(items)
	return newfintfil,items,scores


def getthreshold(fback):
	scores=[]
	for interval in fback :
		scores.append(int(interval.score))
	threshold=(round(max(scores)/2))
	return threshold

def createhistogram(scores_hist,figurename,nbin,file1):
	import numpy
	import pylab
	n,bins,patches=pylab.hist(scores_hist,bins=nbin)
	pylab.title('scores histogram for '+file1)
	pylab.xlabel('background scores')
	pylab.savefig(figurename+'.png')

def DeleteDuplicatesList(x):
    import pybedtools
    inputlist= []
    items =[]
    for interval in x :	
	items.append((interval.chrom,interval.start,interval.end,interval.name,interval.score,interval.strand))
    for a in items:
        if a not in inputlist:
            inputlist.append(a)
    cleanbed=pybedtools.BedTool(inputlist)	
    return cleanbed,inputlist



def deleteandmerge(bedint):
	import pybedtools
# Remove short sequences (l <10)
	nline=0;
	total_lineas=len(bedint)
	outbed=[]
	while nline<total_lineas :
		line=bedint[nline]
		chrom=line[0]
		startp=int(line[1])
		endp=int(line[2])
		name=line[3]
		score=int(line[4])
		strand=line[5]
		startp,endp,score,nreads=check_overlap(bedint,chrom,startp,endp,score,strand,nline)
		nline=nline+nreads 
		if endp-startp > 10 :
			outbed.append((chrom,startp,endp,name,score,strand))				
		nline=nline+1;
	foutbed=pybedtools.BedTool(outbed)
	return foutbed



def check_overlap(lineas,chrom,startp,endp,score,strand,nline):
# check overlapping reads and find the one with the highest score
	nreads=0;
	flag=0;
	while ((flag==0) and ((nline+nreads+1)<len(lineas)) ) :
		line=lineas[nline+nreads+1]
		chromtemp=line[0]
		startptemp=int(line[1])
		endptemp=int(line[2])
		nametemp=line[3]
		scoretemp=int(line[4])
		strandtemp=line[5]
		if ((chromtemp==chrom) and (strandtemp==strand)) :	 
			if ((startp==startptemp) and (endp==endptemp)) :
				nreads=nreads+1;
				if(scoretemp>score) :
					startp,endp,score=startptemp,endptemp,scoretemp
			else :
				if ((startp>=startptemp) and (startp<=endptemp) and (endp>=endptemp)) :
					nreads=nreads+1;
					if(scoretemp>score) :
						startp,endp,score=startptemp,endptemp,scoretemp
				elif ((startp<=startptemp) and (startptemp<=endp) and (endp>=endptemp)) :	
					nreads=nreads+1;
					if(scoretemp>score) :
						startp,endp,score=startptemp,endptemp,scoretemp
				elif ((startp>=startptemp) and (endp<=endptemp) and (startp<=endptemp)) :
					nreads=nreads+1;
					if(scoretemp>score) :
						startp,endp,score=startptemp,endptemp,scoretemp
				elif ((startp<=startptemp) and (startptemp<=endp) and (endp<=endptemp)) :
					nreads=nreads+1;
					if(scoretemp>score) :
						startp,endp,score=startptemp,endptemp,scoretemp
				else :
					flag=1	
		else :
			flag=1
	return startp,endp,score,nreads

	

def process_inputfile(parameters):
	[alpha,histogram,filename,namebed,fbackground,fclusters,threshold]=getparameters(parameters)
	import pybedtools
	print "filebackground --- "+fbackground
	print "fileclusters --- "+fclusters
	fback=pybedtools.BedTool(fbackground)
	fclu =pybedtools.BedTool(fclusters)
# get score threshold (default)
	if threshold==0 :
		threshold=getthreshold(fback)	
# first file must have the sequence
	fint1=fclu.intersect(fback)
	fint2=fback.intersect(fclu)
# sort intersection
	fint1.sort()
	fint2.sort()
# remove duplicates
	fint1,fint1list=DeleteDuplicatesList(fint1)
	fint2,fint2list=DeleteDuplicatesList(fint2)
# check overlap 
	fintfil=deleteandmerge(fint1list)
# obtain smaller sequences
	newfintfil,items,scores=clean_intersection(fintfil,fclu,fint2)
	newfintfil.saveas(namebed)
# statistics from bed
	ci,pvalue=getstatistics(scores,alpha,threshold)
# Histogram from bed	
	if histogram == True :
# n bin set to the maximum score
		nbin=max(scores)
# hist title is the filename, but this can be modified
		histtitle=filename
		createhistogram(scores,filename,nbin,histtitle)
	return ci,threshold,alpha


def int_main(parameters):
	print "fileparameters --- "+parameters
	ci,threshold,alpha=process_inputfile(parameters)
	print "Result:"
	if ci[0]> threshold or ci[1]>threshold :
		print "Unusual presence of common background, ci > Threshold, (Threshold ="+str(threshold)+")"
		print "Confidence interval for the mean of the score (alpha="+str(alpha)+"), ci:"+str(ci) 
	else :
		print "Usual presence of common background, ci < Threshold, (Threshold ="+str(threshold)+")"
		print "Confidence interval for the mean of the score (alpha="+str(alpha)+"), ci:"+str(ci) 


#----------------------------------MAIN-------------------------------------------

 
import sys
int_main(sys.argv[1])

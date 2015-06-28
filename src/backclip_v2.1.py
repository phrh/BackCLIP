#-------------------------------------------------------------------------------------
#!/usr/bin/env python
# USAGE: python backclip.py filebackground fileclusters fileparameters
# Use this program to have a quantitative measure of the amount of background present in a determined dataset. The presence of background is considerable 
# INPUT:         fileparameters -> file with parameters such as:
#		 alpha: significance level(default value 0.01)
#		 threshold: to define if the amount of common background is significant (default (maximum score)/2) 
#		 histogram: give as output histogram of the intersection scores (default false)
#		 remove: remove sites in the intersection with score > threhsold from dataset. 
#		 	The output file name is clusters file name .backclip.bed (default false)
#		 filename: in case histogram is true, the name of the output file (default fileclusters.histbackground)	
# OUTPUT: confidence interval (significance level 99%)
# 	  histogram score distribution in the intersection file
#	  According to the user input parameters, the dataset without the sites in the intersection (score > threshold)		
#-------------------------------FUNCTIONS-LIST----------------------------------------

def getparameters(parameters):
	from datetime import date
	fparam=open(parameters,'r')
	lineas=fparam.readlines()	
#default parameters in case are deleted from parameters file
	alpha=0.01
	histogram=False
	smotif=False
	motif=''
	filename='histbackground '+str(date.today().ctime())
	namebed=filename+'intersection.bed'
	remove=False
	
	for line in lineas:
		paramname,value=line.split("=")
		if paramname.lower().strip() == "alpha" :
			alpha=float(value)
		elif paramname.lower().strip() == "histogram" :
			if value.lower().strip() == "true" :
				histogram=True
		elif paramname.lower().strip() == "filename" :		
			filename=value.strip()
		elif paramname.lower().strip() == "smotif" :
			smotif=True
		elif paramname.lower().strip() == "motif" :
			motif=value.strip()
		elif paramname.lower().strip() == "namebed" :
			namebed=value.strip()+'.intersection.bed'
		elif paramname.lower().strip() == "fbackground" :
			fbackground=value.strip()
		elif paramname.lower().strip() == "fclusters" :
			fclusters=value.strip()
		elif paramname.lower().strip() == "remove" :
			if value.lower().strip() == "true" :
				remove=True
	return alpha,histogram,filename,smotif,motif,namebed,fbackground,fclusters,remove


def getstatistics(scores,alpha,threshold):
	import scipy
	from scipy import stats	
	import math
	if len(scores) >= 30:
#Indicator 1 - Confidence interval for the proportion
		hscores=filterarray(scores,threshold)
		lenhscores=len(hscores)
		lenscores=len(scores)
		ascores=scipy.array(scores)	
		p=float(lenhscores)/float(lenscores)
		x=math.sqrt(p*(1-p)/float(lenscores))		
		if alpha==0.01 :
			factor=2.58
			ci=[p-factor*x,p+factor*x]
			ci=[x*100 for x in ci]
		elif alpha==0.05 :
			factor=1.96
			ci=[p-factor*x,p+factor*x]
			ci=[x*100 for x in ci]
		else :
			print "sample size must be >30"
			ci=[]
	return ci

def filterarray(array,filternumber):
	import numpy
	myarray=numpy.asarray(array)
	arrayfilteredindexes=myarray>filternumber
	arrayfiltered=myarray[arrayfilteredindexes].tolist()
	return arrayfiltered

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

def clean_intersection(fintfil,fclu,fint2,smotif):
	import pybedtools
	items=[]
	scores=[]
	for interval in fintfil :
		chrom,startp,endp,seq,score,strand=interval.chrom,interval.start,interval.end,interval.name,interval.score,interval.strand
		if smotif == True :
			newseq=findseq(chrom,startp,endp,seq,strand,fclu)
		else :
			newseq='NA'
		newscore=findscore(chrom,startp,endp,seq,strand,fint2)
		items.append((chrom,startp,endp,newseq,newscore,strand))
		scores.append(int(newscore))
	newfintfil=pybedtools.BedTool(items)
	return newfintfil,items,scores

def filterarray(array,filternumber):
	import numpy
	myarray=numpy.asarray(array)
	arrayfilteredindexes=myarray>filternumber
	arrayfiltered=myarray[arrayfilteredindexes].tolist()
	return arrayfiltered

# score based on the 99 percentile
def getthreshold(fback):
	import scipy 
	from scipy import stats	
	scores=[]
	ifilter=1
	for interval in fback :
		scores.append(int(interval.score))
	array2=filterarray(scores,ifilter)
	spscoresback=scipy.array(array2)
	percentile=stats.scoreatpercentile(array2,99)
	return percentile	


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


def countmotifs(fclu):
	iupacdict = {'M':'[AC]','R':'[AG]','W':'[AT]','S':'[CG]','Y':'[CT]','K':'[GT]','V':'[ACG]','H':'[ACT]','D':'[AGT]','B':'[CGT]','X':'[ACGT]','N':'[ACGT]'}    
        for key, value in iupacdict.iteritems():
                motif=motif.replace(key,value)	


def deleteandmerge(bedint):
	import pybedtools
# Remove short sequences (l <10()
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

#def process_inputfile(fbackground,fclusters,parameters):
def process_inputfile(parameters):
	[alpha,histogram,filename,smotif,motif,namebed,fbackground,fclusters,remove]=getparameters(parameters)
	import pybedtools
	print "filebackground --- "+fbackground
	print "fileclusters --- "+fclusters
	fback=pybedtools.BedTool(fbackground)
	fclu =pybedtools.BedTool(fclusters)
	totalclusters=len(fclu)
# get score threshold (default) - 99 percentile
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
	newfintfil,items,scores=clean_intersection(fintfil,fclu,fint2,smotif)
	newfintfil.saveas(namebed)
# Remvove sites above threshold from dataset 
	if remove == True : 
		newdataset=fclu.subtract(newfintfil)
		newdataset.saveas(fclusters+'.backclip.bed')
# Indicator 2
	rsize=float(float(len(scores))/float(totalclusters))*100
	ci=getstatistics(scores,alpha,threshold)
# Histogram from bed	
	if histogram == True :
# n bin set to the maximum score
		nbin=max(scores)
# hist title is the filename, but this can be modified
		histtitle=filename
		createhistogram(scores,filename,nbin,histtitle)
	return ci,threshold,alpha,rsize


def int_main(parameters):
	print "fileparameters --- "+parameters
	ci,threshold,alpha,rsize=process_inputfile(parameters)
	print "Result:"	
	if ci[0]> 50 or ci[1]>50 :
		print "Unusual presence of common background, ci proportion> Threshold, (Threshold =50%)"
		print "Confidence interval for the proportion of sites with a score higher than the threshold(alpha="+str(alpha)+"), ci:"+str(ci) 
	else :
		print "Confidence interval for the proportion of sites with a score higher than the threshold(alpha="+str(alpha)+"), ci:"+str(ci)
	if rsize>30 :
		print "Unusual number of sites in the intersection (with the common background) relative to the number of clusters"+str(rsize)
	else :
		print "number of sites in the intersection (with the common background) relative to the number of clusters"+str(rsize)


#----------------------------------MAIN-------------------------------------------

 
import sys
int_main(sys.argv[1])

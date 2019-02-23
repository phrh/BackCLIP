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
#
#
#Created by Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra, Santiago Herrera on 5 February 2015 
#Copyright (c) 2015 Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra, Santiago Herrera. All rights reserved.
#BackCLIP is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 2.
#--------------------------------------------------------------------------------------


from datetime import date
import scipy
from scipy import stats
import math
import numpy		
import pylab
import pybedtools


"""
This method takes the parameters defined in setup file and save them into execution variables
"""
def get_parameters(parameters):
	#access to the file
	fparam = open(parameters,'r')
	lines = fparam.readlines()	
	
	#default parameters in case are deleted from parameters file
	alpha = 0.01
	histogram = False
	smotif = False
	motif = ''
	filename = 'histbackground ' + str(date.today().ctime())
	namebed = filename + 'intersection.bed'
	remove = False
	
	#read line by line of the file, and assign each value depending of variable setup   
	for line in lines:
		paramname,value = line.split("=")
		
		if paramname.lower().strip() == "alpha":
			alpha = float(value)
		elif paramname.lower().strip() == "histogram":
			if value.lower().strip() == "true" :
				histogram = True
		elif paramname.lower().strip() == "filename":		
			filename = value.strip()
		elif paramname.lower().strip() == "smotif":
			smotif = True
		elif paramname.lower().strip() == "motif":
			motif = value.strip()
		elif paramname.lower().strip() == "namebed":
			namebed = value.strip()+'.intersection.bed'
		elif paramname.lower().strip() == "fbackground":
			fbackground = value.strip()
		elif paramname.lower().strip() == "fclusters":
			fclusters = value.strip()
		elif paramname.lower().strip() == "remove":
			if value.lower().strip() == "true":
				remove = True
				
	return alpha, histogram, filename, smotif, motif, namebed, fbackground, fclusters, remove



'''
'''
def get_statistics(scores, alpha, threshold):
	if len(scores) >= 30:
		#Indicator 1 - Confidence interval for the proportion
		hscores = filter_array(scores,threshold)
		lenhscores = len(hscores)
		lenscores = len(scores)
		ascores = scipy.array(scores)	
		p = float(lenhscores) / float(lenscores)
		x = math.sqrt(p * (1-p) / float(lenscores))		
	
		if alpha == 0.01:
			factor = 2.58
			ci = [p - factor * x, p + factor * x]
			ci=[x * 100 for x in ci]
		elif alpha == 0.05:
			factor = 1.96
			ci=[p - factor * x, p + factor * x]
			ci=[x * 100 for x in ci]
		else:
			print("sample size must be > 30")
			ci = []
	
	return ci



'''
'''
def filter_array(array, filter_number):
	my_array = numpy.asarray(array)
	array_filtered_indexes = my_array > filter_number
	array_filtered = my_array[ array_filtered_indexes ].tolist()
	
	return array_filtered



'''
'''
def find_seq(chrom, startp, endp, seq, strand, fclu):
	flag=0
	for interval in fclu:
		if seq in interval.name:
			cchrom, cstartp, cendp, cseq, cscore, cstrand = interval.chrom, interval.start, interval.end, interval.name, interval.score, interval.strand
			if (cstartp <= startp) and (endp <= cendp) and (cchrom == chrom): 
				new_seq = cseq[startp - cstartp : endp - cstartp + 1]
				flag = 1
							
	return new_seq



'''
'''
def find_score(chrom, startp, endp, seq, strand, fint2):
	new_score = 0
	for interval in fint2 :
		if (startp == interval.start) and (endp == interval.end) and (chrom == interval.chrom): 
			new_score = interval.score

	return new_score



'''
'''
def clean_intersection(fintfil, fclu, fint2, smotif):
	items = []
	scores = []
	for interval in fintfil :
		chrom, startp, endp, seq, score, strand = interval.chrom, interval.start, interval.end, interval.name, interval.score, interval.strand
		if smotif == True:
			new_seq = find_seq(chrom, startp, endp, seq, strand, fclu)
		else :
			new_seq = 'NA'
		
		new_score = find_score(chrom, startp, endp, seq, strand, fint2)
		items.append((chrom, startp, endp, new_seq, new_score, strand))
		scores.append(int(new_score))
		
	new_fintfil = pybedtools.BedTool(items)
	return new_fintfil, items, scores




'''
score based on the 99 percentile
'''
def get_threshold(fback):
	scores = []
	ifilter = 1
	for interval in fback :
		scores.append(int(interval.score))
	
	array2 = filter_array(scores, ifilter)
	spscoresback = scipy.array(array2)
	percentile = stats.scoreatpercentile(array2, 99)
	print("ads " + str(percentile))
	
	return percentile	



'''
'''
def create_histogram(scores_hist, figurename, nbin, file):
	n, bins, patches = pylab.hist(scores_hist, bins = nbin)
	pylab.title('Scores histogram for ' + file)
	pylab.xlabel('Background scores')
	pylab.savefig(figurename + '.png')



'''
'''
def delete_duplicates_list(x):
    inputlist= []
    items =[]
    for interval in x:	
        items.append((interval.chrom, interval.start, interval.end, interval.name, interval.score, interval.strand))
    
    for a in items:
        if a not in inputlist:
            inputlist.append(a)
    
    cleanbed = pybedtools.BedTool(inputlist)	
    return cleanbed, inputlist



'''
'''
def count_motifs(fclu):
	iupacdict = {'M':'[AC]','R':'[AG]','W':'[AT]','S':'[CG]','Y':'[CT]','K':'[GT]','V':'[ACG]','H':'[ACT]','D':'[AGT]','B':'[CGT]','X':'[ACGT]','N':'[ACGT]'}    
	for key, value in iupacdict.iteritems():
		motif=motif.replace(key,value)	



'''
'''
def delete_and_merge(bedint):
	# Remove short sequences (l <10()
	nline=0;
	total_lineas = len(bedint)
	outbed = []
	while nline<total_lineas :
		line = bedint[nline]
		chrom = line[0]
		startp = int(line[1])
		endp = int(line[2])
		name = line[3]
		score = int(line[4])
		strand = line[5]
		startp, endp, score, nreads =check_overlap(bedint, chrom, startp, endp, score, strand, nline)
		nline = nline + nreads 
		if (endp - startp) > 10 :
			outbed.append((chrom, startp, endp, name, score, strand))
							
		nline=nline+1;
		
	foutbed = pybedtools.BedTool(outbed)
	return foutbed



'''
'''
def check_overlap(lines, chrom, startp, endp, score, strand, nline):
	# check overlapping reads and find the one with the highest score
	nreads = 0
	flag = 0
	while (flag == 0) and ((nline + nreads + 1) < len(lines)):
		line = lines[nline + nreads + 1]
		chrom_temp = line[0]
		startp_temp = int(line[1])
		endp_temp = int(line[2])
		name_temp = line[3]
		score_temp = int(line[4])
		strand_temp = line[5]
		
		if (chrom_temp == chrom) and (strand_temp == strand):	 
			if (startp == startp_temp) and (endp == endp_temp):
				nreads += 1
				if score_temp > score:
					startp, endp, score = startp_temp, endp_temp, score_temp
			else:
				if (startp >= startp_temp) and (startp <= endp_temp) and (endp >= endp_temp):
					nreads += 1
					if score_temp > score:
						startp, endp, score = startp_temp, endp_temp, score_temp
				elif (startp <= startp_temp) and (startp_temp <= endp) and (endp >= endp_temp):	
					nreads += 1
					if score_temp > score:
						startp, endp, score = startp_temp, endp_temp, score_temp
				elif (startp >= startp_temp) and (endp <= endp_temp) and (startp <= endp_temp):
					nreads += 1
					if score_temp > score:
						startp, endp, score = startp_temp, endp_temp, score_temp
				elif (startp <= startp_temp) and (startp_temp <= endp) and (endp <= endp_temp) :
					nreads += 1
					if score_temp > score:
						startp, endp, score = startp_temp, endp_temp, score_temp
				else:
					flag = 1	
		else:
			flag=1
	return startp, endp, score, nreads



'''
This method takes input file and apply over it a cleaning process
'''
def process_inputfile(parameters):
	[alpha, histogram, filename, smotif, motif, namebed, fbackground, fclusters, remove] = get_parameters(parameters)
	print("file background --- %s" % fbackground)
	print("file clusters --- %s" % fclusters)
	fback = pybedtools.BedTool(fbackground)
	fclu = pybedtools.BedTool(fclusters)
	total_clusters = len(fclu)
	
	# get score threshold (default) - 99 percentile
	threshold = get_threshold(fback)
		
	# first file must have the sequence
	fint1 = fclu.intersect(fback)
	fint2 = fback.intersect(fclu)
	
	# sort intersection
	fint1.sort()
	fint2.sort()

	# remove duplicates
	fint1, fint1list = delete_duplicates_list(fint1)
	fint2, fint2list = delete_duplicates_list(fint2)
	
	# check overlap 
	fintfil = delete_and_merge(fint1list)

	# obtain smaller sequences
	newfintfil, items, scores= clean_intersection(fintfil, fclu, fint2, smotif)
	newfintfil.saveas(namebed)

	# Remvove sites above threshold from dataset 
	if remove: 
		new_dataset = fclu.subtract(newfintfil)
		new_dataset.saveas(fclusters + '.backclip.bed')
		
	# Indicator 2
	rsize = float(float(len(scores)) / float(total_clusters)) * 100
	ci = get_statistics(scores, alpha, threshold)

	# Histogram from bed	
	if histogram:
		# n bin set to the maximum score
		nbin = max(scores)
		# hist title is the filename, but this can be modified
		histtitle = filename
		create_histogram(scores, filename, nbin, histtitle)
	
	return ci, threshold, alpha, rsize


"""
This method is the initial called of execution
"""
def int_main(parameters):
	print("file parameters ---  %s" % parameters)
	
	ci,threshold,alpha,rsize = process_inputfile(parameters)
	print("Result:")	
	
	if ci[0] > 50 or ci[1] > 50:
		print("Unusual presence of common background, ci proportion> Threshold, (Threshold =50%)")
		print("Confidence interval for the proportion of sites with a score higher than the threshold(alpha = %s ), ci: %s" % (str(alpha), str(ci)) )
	else:
		print("Confidence interval for the proportion of sites with a score higher than the threshold(alpha = %s ), ci: %s" % (str(alpha), str(ci)) )
	if rsize > 30:
		print("Unusual number of sites in the intersection (with the common background) relative to the number of clusters %s" % str(rsize))
	else:
		print("Number of sites in the intersection (with the common background) relative to the number of clusters %d" % str(rsize))



#----------------------------------MAIN-------------------------------------------
import sys
int_main(sys.argv[1]) #calling main method
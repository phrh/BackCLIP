#-------------------------------------------------------------------------------------
#!/usr/bin/env python
# Usage: python delete_and_merge_minimum.py filebackground result minimum
#This file prints only one record per overlapping reads and assings the maximum score
#Please sort the background bed file 
#sort -k1,1 -k2n,2n GLOBALBACKGROUNDGROUPS.bed > GLOBALBACKGROUNDGROUPS.bed.sorted
#-------------------------------FUNCTIONS-LIST----------------------------------------
#def process_inputfile(file1,fout,minimum):
#def int_main(file1,fout):
#---------------------------------FUNCTIONS-------------------------------------------

def check_overlap(lineas,chrom,startp,endp,score,strand,nline):
# check overlapping reads and find the one with the highest score
	nreads=0;
	flag=0;
	while ((flag==0) and ((nline+nreads+1)<len(lineas)) ) :
		templine=lineas[nline+nreads+1]
		chromtemp, startptemp, endptemp,nametemp,scoretemp,strandtemp = templine.split("\t")
		startptemp, endptemp, scoretemp = int(startptemp), int(endptemp), int(scoretemp)
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
			


def process_inputfile(file1,fout,minimo):
	g = open(file1,'r')
	outfile = open(fout,'a+')
	lineas=g.readlines()
	nline=0;
	total_lineas=len(lineas)
	while nline<total_lineas :
		line=lineas[nline]
		chrom, startp, endp,name,score,strand = line.split("\t")
		startp, endp, score = int(startp), int(endp), int(score)
		startp,endp,score,nreads=check_overlap(lineas,chrom,startp,endp,score,strand,nline)
		#change line number - skipping overlapping reads
		nline=nline+nreads 
		if endp-startp > int(minimo) :
			outfile.write(chrom+'\t'+str(startp)+'\t'+str(endp)+'\tNA\t'+str(score)+'\t'+strand)				
		nline=nline+1;
	g.close()
	outfile.close()
#MAIN
def int_main(file1,fout,minimo):
	print "file1 - "+file1
	print "fout -"+fout
	process_inputfile(file1,fout,minimo)

	
#----------------------------------MAIN-------------------------------------------

import sys
int_main(sys.argv[1],sys.argv[2],sys.argv[3])



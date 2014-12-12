#-------------------------------------------------------------------------------------
#!/usr/bin/env python
# Usage: python count_motif_occurrences.py motif filebed outfile
#This program finds the number of motif occurrences in a bed format file (name is the sequence). However, the maximum number of motif occurrences in each cluster is 1..
#-------------------------------FUNCTIONS-LIST----------------------------------------
#def process_inputfile(motif,fin):
#def int_main(motif,fin):
#---------------------------------FUNCTIONS-------------------------------------------

def findocurrences(motifseqs,seq):
# find motif occurrences and count the number of occurrences
	import re
	i=[]
	for motifseq in motifseqs:
		occurrences=len(re.findall(motifseq,seq))
		i.append(occurrences)
	return i		

def find_paths(sequences,path,pospath):
    templates=sequences
    sequences=[]
    for j in range(len(templates)):
        temp=templates[j]
        for i in range(len(path)):
            seq=temp[:pospath]+path[i]+temp[pospath+1:]
            sequences.append(seq)
    return  sequences   

def obtain_paths(pattern):    
        #obtain all possible paths for a motif
        paths=[]
        npaths=[]
        pospaths=[]
        tempstring=[]
        len1=len(pattern)
        flag=0
        match=''
        count=0
        #Obtaining template(match), paths and positions
        for i in range(len1):
            if flag==1 :
                if pattern[i]==']':
                    paths.append(tempstring)
                    npaths.append(len(tempstring))
                    tempstring=[]
                    flag=0                
                else:
                    tempstring.append(pattern[i])
            elif pattern[i]=='[':
                 count=count+1
                 flag=1
                 match=match+'-'
                 pospaths.append(count-1)
            else:
                match=match+pattern[i]
                count=count+1
        #Finding paths and saving everything in the variable sequences
        sequences=[]
        sequences.append(match)    
        for j in range(len(paths)):
            sequences=find_paths(sequences,paths[j],pospaths[j])
        return sequences     



def reverse_complement(seq):    
	alt_map = {'ins':'0'}
#PAULA Prueba esto!!!!!
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A','[':']',']':'['} 
	for k,v in alt_map.iteritems():
        	seq = seq.replace(k,v)
	bases = list(seq) 
   	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)
	for k,v in alt_map.iteritems():
        	bases = bases.replace(v,k)
	return bases

def obtain_unique_match(motif):
        #Replacing IUPAC dictionary
	iupacdict ={'U':'T','M':'[AC]','R':'[AG]','W':'[AT]','S':'[CG]','Y':'[CT]','K':'[GT]','V':'[ACG]','H':'[ACT]','D':'[AGT]','B':'[CGT]','X':'[ACGT]','N':'[ACGT]'}    
	for key, value in iupacdict.iteritems():
		motif=motif.replace(key,value)
	#obtain paths from iupac
	seqs=obtain_paths(motif)
	#obtain reverse complement
	rmotif=reverse_complement(motif)
	rseqs=obtain_paths(rmotif)
	seqs.extend(rseqs)
	#obtain unique sequences
	sequences=list(set(seqs))
	return sequences
	


def process_inputfile(motif,fin,foutname):
	g = open(fin,'r')
	fout=open(foutname,'w')	
	lineas=g.readlines()
	clustersone=0
	clusterstwo=0
	clustersmore=0
	clustersonerev=0
	clusterstworev=0
	clustersmorerev=0
	motifseqs=obtain_unique_match(motif)
	for m in motifseqs :
		fout.write(m+'\t')
	fout.write('\n')
	for line in lineas:
		chrom, startp, endp,seq,score,strand = line.split("\t")
		moccurrences=findocurrences(motifseqs,seq)
		for occurrences in moccurrences :
			fout.write(str(occurrences)+'\t')
		fout.write('\n')
	g.close()
	fout.close()

#MAIN
def int_main(motif,fin,fout):
	print "motif - "+motif
	print "fin -"+fin
	print "fout -"+fout
	import re
	motif=re.sub('U','T',motif)
	process_inputfile(motif,fin,fout)

	
#----------------------------------MAIN-------------------------------------------

import sys
int_main(sys.argv[1],sys.argv[2],sys.argv[3])



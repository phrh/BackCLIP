# count_motif_occurrences - (extra) if you have information regarding the motif, this program counts
#				    the motif occurrences in the input bed
#
# Created by Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra, Santiago Herrera on 5 February 2015
# Copyright (c) 2015 Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra, Santiago Herrera. All rights reserved.
#
#
# count_motif_ocurrences is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 2.
#-------------------------------------------------------------------------------------
#!/usr/bin/env python
# Usage: python count_motif_occurrences.py motif filebed outfile
# This program finds the number of motif occurrences in a bed format file (name is the sequence). 
#-------------------------------FUNCTIONS-LIST----------------------------------------
#	process_inputfile(motif,fin)
#	int_main(motif,fin)
#	find_ocurrences(motifseqs,seq)
#	find_paths(sequences,path,pospath)
#	obtain_paths(pattern)
#	reverse_complement(seq)
#	obtain_unique_match(motif)
#	 
#---------------------------------FUNCTIONS-------------------------------------------
import re


'''
'''
def find_ocurrences(motif_seqs, seq):
	# find motif occurrences and count the number of occurrences
	i = []
	for motif_seq in motif_seqs:
		occurrences = len(re.findall(motif_seq, seq))
		i.append(occurrences)
	return i		



'''
'''
def find_paths(sequences, path, pospath):
    templates = sequences
    sequences = []
    for j in range(len(templates)):
        temp = templates[j]
        for i in range(len(path)):
            seq = temp[:pospath] + path[i] + temp[pospath + 1:]
            sequences.append(seq)
    return  sequences   



'''
'''
def obtain_paths(pattern):    
        #obtain all possible paths for a motif
        paths = []
        npaths = []
        pospaths = []
        temp_string = []
        len1 = len(pattern)
        flag = 0
        match = ''
        count = 0
        
        #Obtaining template(match), paths and positions
        for i in range(len1):
            if flag == 1:
                if pattern[i] == ']':
                    paths.append(temp_string)
                    npaths.append(len(temp_string))
                    temp_string = []
                    flag = 0                
                else:
                    temp_string.append(pattern[i])
            elif pattern[i] == '[':
                 count = count + 1
                 flag = 1
                 match = match + '-'
                 pospaths.append(count - 1)
            else:
                match = match + pattern[i]
                count = count+1
        #Finding paths and saving everything in the variable sequences
        sequences = []
        sequences.append(match)    
        for j in range(len(paths)):
            sequences = find_paths(sequences, paths[j], pospaths[j])
   
        return sequences     



'''
'''
def reverse_complement(seq):    
	alt_map = {'ins':'0'}
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','U':'A','[':']',']':'['}
	 
	for k,v in alt_map.iteritems():
        	seq = seq.replace(k,v)
	
	bases = list(seq) 
   	bases = reversed([complement.get(base,base) for base in bases])
	bases = ''.join(bases)
	
	for k,v in alt_map.iteritems():
        	bases = bases.replace(v,k)
	
	return bases



'''
'''
def obtain_unique_match(motif):
        #Replacing IUPAC dictionary
	iupacdict = {'U':'T','M':'[AC]','R':'[AG]','W':'[AT]','S':'[CG]','Y':'[CT]','K':'[GT]','V':'[ACG]','H':'[ACT]','D':'[AGT]','B':'[CGT]','X':'[ACGT]','N':'[ACGT]'}    

	for key, value in iupacdict.iteritems():
		motif = motif.replace(key,value)
	
	#obtain paths from iupac
	seqs=obtain_paths(motif)
	
	#obtain reverse complement
	rmotif = reverse_complement(motif)
	rseqs = obtain_paths(rmotif)
	seqs.extend(rseqs)
	
	#obtain unique sequences
	sequences = list(set(seqs))
	
	return sequences
	


'''
'''
def process_inputfile(motif, f_in, f_outname):
	g = open(f_in,'r')
	f_out = open(f_outname,'w')	
	lines = g.readlines()
	clusters_one = 0
	clusters_two = 0
	clusters_more = 0
	clusters_onerev = 0
	clusters_tworev = 0
	clusters_morerev = 0
	motif_seqs = obtain_unique_match(motif)
	
	for m in motif_seqs:
		f_out.write(m + '\t')
	f_out.write('\n')
	
	for line in lines:
		chrom, startp, endp, seq, score, strand = line.split("\t")
		m_occurrences = find_ocurrences(motif_seqs, seq)
		for occurrences in m_occurrences :
			f_out.write(str(occurrences) + '\t')
		f_out.write('\n')
	g.close()
	f_out.close()



'''
'''
def int_main(motif,f_in,f_out):
	print("motif - %s" % motif)
	print("fin -%s" % f_in)
	print("fout - %s" + f_out)
	motif = re.sub('U','T', motif)
	process_inputfile(motif, f_in, f_out)

	
#----------------------------------MAIN-------------------------------------------
import sys
int_main(sys.argv[1], sys.argv[2], sys.argv[3])
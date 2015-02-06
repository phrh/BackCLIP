BackCLIP
========

Tools to identify background presence in PAR-CLIP datasets


----------------
#### Requirements

- Python 2.7 and above
- Python  Packages (To install plese use: pip install -r /path/to/[requirements.txt](https://github.com/phrh/BackCLIP/blob/master/requirements.txt))

----------------
#### Install

Download python scritp (scr/[backclip_v1.0.py](https://github.com/phrh/BackCLIP/blob/master/src/backclip_v1.0.py))

----------------
#### Usage


- **backclip_v1.0.py**.  This python program identifies the presence of common background in a PAR-CLIP dataset
- shell script will search all the restriction sites from the input file (patternfilename) in every genome from the input file (genomefilename). As a result the script provides the following files:
 
	- ALL.aligned.txt, ALL.failed.txt, ALL.processed.txt,  ALL.suppressed.txt - each file with a table summarizing bowtie output(reads aligned, failed, processed and suppressed) for each genome.
	- ALL.count.txt - contains a table with the number of restriciton sites found in each genome
	- ALL.size.txt - contains a table with the size of each genome

	The input arguments are: 
	- genomefilename: name of file with table with two columns (1) species code and (2) link to whole genome fasta file 
	- 
	

#### License

Created by Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra on 5 February 2015
Copyright (c) 2015 Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra. All rights reserved.

BackCLIP is free software: you can redistribute it and/or modify  it under the terms of the GNU General Public License as published by # the Free Software Foundation, version 2.


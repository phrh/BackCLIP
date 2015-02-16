BackCLIP
========

The PAR-CLIP protocol derives a transcriptome wide set of binding sites for RNA Binding Proteins. However,   non-specific RNA background remains. We propose a tool, BackCLIP, to identify the presence of common RNA background in a PAR-CLIP dataset. 

----------------
### Content

* Supplementary data
* Requirements
* Install
* Usage
* License
* Extra

----------------
#### Supplementary data

* Details for Background Initial Set and Additional Results [_Supplementary_Information_BACKCLIP.doc_](https://github.com/phrh/BackCLIP/blob/master/SupplementaryData/Supplementary_%20Information_BACKCLIP.doc)
* CommonBackground: _BackgroundTraining_19datasets.bed_ (Common background set used in the paper),  _Background_49datsets.bed_ (Common Background built on 49 datasets)


----------------
#### Requirements

- Python 2.7 and above
- Python  Packages (To install plese (1) download  [requirements.txt](https://github.com/phrh/BackCLIP/blob/master/requirements.txt) and (2) use the command  pip install -r /path/to/requirements.txt )

----------------
#### Install

Download python scritp (scr/[backclip_v1.0.py](https://github.com/phrh/BackCLIP/blob/master/src/backclip_v1.0.py))

----------------
#### Usage

- **backclip_v1.0.py**.  This python program identifies the presence of common background in a PAR-CLIP dataset. As a result provides:
	- confidence interval for the mean distribution 
	- histogram score distribution in the intersection file
	- common background in a dataset (bed file)
	
- The input argument is: 
	- parameters file name  (see _src/[parameters_v1.0](https://github.com/phrh/BackCLIP/blob/master/src/parameters_v1.0)_): this file contains the following parameters :
		- alpha=significance level (default value 0.01)
		- threshold=to define if the amount of common background is significant (default (maximum score)/2)
		- histogram=give as output histogram of the intersection scores (default false)
		- fbackground=(.BED)file background and corresponding scores (see _src/example/[GLOBALBACKGROUNDGROUPS.bed.sorted.delete.min10](https://github.com/phrh/BackCLIP/blob/master/src/example/GLOBALBACKGROUNDGROUPS.bed.sorted.delete.min10)_. Nevertheless, the complete version of the background is _CommonBackground/[Background_49datsets.bed](https://github.com/phrh/BackCLIP/blob/master/CommonBackground/Background_49datsets.bed)_)
		- fclusters=(.BED)file with clusters detected from PAR-CLIP dataset (see _src/example/[QKI_SRR048972_Bowtie_Score_Cleanmin10_Prueba1.bed](https://github.com/phrh/BackCLIP/blob/master/src/example/QKI_SRR048972_Bowtie_Score_Cleanmin10_Prueba1.bed)_)
		- filename=in case histogram is true, the name of the output file (default fileclusters.histbackground)
		- namebed=name bed with set of common background found in fileclusters
		
		
	To run, just write on shell

	_python backclip_v1.0.py parametersfilename_

----------------
#### Extra

- ClusterDetection (CD_Bg.jar)
- _count_motif_occurrences.py_

----------------
#### License

Created by Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra, Santiago Herrera on 5 February 2015
Copyright (c) 2015 Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra, Santiago Herrera. All rights reserved.

BackCLIP is free software: you can redistribute it and/or modify  it under the terms of the GNU General Public License as published by # the Free Software Foundation, version 2.


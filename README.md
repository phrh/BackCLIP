BackCLIP
========

Tools to identify background presence in PAR-CLIP datasets


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
	
- The input arguments is: 
	- parameters file name: this file contains the following parameters:
		- alpha=significance level (default value 0.01)
		- threshold=to define if the amount of common background is significant (default (maximum score)/2)
		- histogram=give as output histogram of the intersection scores (default false)
		- fbackground=(.BED)file background and corresponding scores
		- fclusters=(.BED)file with clusters detected from PAR-CLIP dataset
		- filename=in case histogram is true, the name of the output file (default fileclusters.histbackground)
		- namebed=name bed with set of common background found in fileclusters
	
----------------

#### License

Created by Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra on 5 February 2015
Copyright (c) 2015 Paula H. Reyes-Herrera, Cesar A. Speck Hernandez, Carlos Sierra. All rights reserved.

BackCLIP is free software: you can redistribute it and/or modify  it under the terms of the GNU General Public License as published by # the Free Software Foundation, version 2.


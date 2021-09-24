####
# Using this script, you will be able to run Neighborhood Enrichment from Sandra's ChromHMM dataset (10 states) in your bed dataset
# This script should output png, svg, and txt files that display the enrichment state realtive to your input bed file
from __future__ import division

import pandas as pd
import numpy as np
import os
import sys
import subprocess
from subprocess import Popen

#Modify the path to your bed files

pathbed=' /oak/stanford/groups/oro/anncoll/venn/'

#Where you want to save. The last line is the new title name, not the target folder.

save=' /oak/stanford/groups/oro/anncoll/venn/ENCODE'

#Here is where you put the name of your bed file, without the .bed extension. Bed files have to be chr, start, stop, ID. Of course, the extension has to be .bed

beds = ['DMP_hg38_2kb']

cmds=[]

#can modify names (ns and a) depending on how you structured your cell table or on how many states you used.
#the "pos" path is the directory where you binarizebed output is stored
#add the -s option to change the x axis limits
for bed in beds:
	ns=['D7WT']
	for n in ns:
		ChromHMM='/oak/stanford/groups/oro/anncoll/NextSeqChIP/chromHMM/ChromHMM.jar'
		pos='/oak/stanford/groups/oro/anncoll/chromHMM/hg38/states/10/POSTERIOR'
		cell='java -Xmx256m -jar '+ ChromHMM+' NeighborhoodEnrichment -a D7WT -s 200 -posterior -colfields 0,1 -nostrand '+ pos+' '+ pathbed+bed+'.bed' +' ' +save+bed+'_'+n
		cmds.append(cell)

for cmd in cmds:
	process = Popen(cmd.split(), stdout=subprocess.PIPE)
	process.wait()
	print ('run command: %s'%(cmd))

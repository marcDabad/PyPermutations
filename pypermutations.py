#!/usr/bin/env python


"""

PyPermutations has been designed to perform permutations test in situations with big datasets and two groups. Finally it executes a FDR.
Initially has been designed to perform DE.

More information in github.io/marcDabad/PyPermutations

USAGE:
  pypermutations2.py -f <inputFile> [ --c1 <condition1> --c2 <condition2> -n <int> -H <int> -s <statistic> ] -o <outFile>
    
Inputs:
  -f, --file <inputFile>        Space separated table
    
Outputs:
  -o, --output <outFile>        Output table

Options:
  --c1	<condition1>            Condition that determine one of two groups [default: case]
  --c2	<condition2>            Conditions that determine other group [default: control]
  -n	<int>                   Maximum number of permutations [default: 1000000]
  -H	<int>                   Maximum hits to stop permutations [default: 5]
  -s, --stats	<statistic>     Statistic to compute pvalue median|perc25|perc75|ratio_median|ratio_perc25|ratio_perc75 [default: median]
    
"""

import numpy as np
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')
import pandas as pd
import sys
from docopt import docopt
args=docopt(__doc__)


#####FUNCTIONS##########################################################################################################

def median_diff(n):
    return abs(f[n][cd1].median()-f[n][cd2].median())

def percentil_diff(n, percen):
    return abs(np.percentile(f[n][cd1], percen)-np.percentile(f[n][cd2], percen))

def median_ratio(n):
    return f[n][cd1].median()/f[n][cd2].median()

def percentil_ratio(n,percen):
    return np.percentile(f[n][cd1], percen)/np.percentile(f[n][cd2], percen)

def compare_pv(gene):
    pv=pvalue(gene)
    ori=out_df[col_out[0]][gene]
    to_return=0
    if args['--stats']=="median" or args['--stats']=="perc75" or args['--stats']=="perc25":
	if pv>=ori:
	    to_return=1
    elif args['--stats']=="ratio_median" or args['--stats']=="ratio_perc75" or args['--stats']=="ratio_perc25":
	if ori>1:
	    if pv>=ori:
		to_return=1
	else:
	    if pv<=ori:
		to_return=1
    return to_return
    
def pvalue(gene):
    if args['--stats']=='median':
	return median_diff(gene)
    elif args['--stats']=="perc25":
	return percentil_diff(gene, 25)
    elif args['--stats']=="perc75":
	return percentil_diff(gene, 75)
    elif args['--stats']=='ratio_median':
	return median_ratio(gene)
    elif args['--stats']=="ratio_perc25":
	return percentil_ratio(gene, 25)
    elif args['--stats']=="ratio_perc75":
	return percentil_ratio(gene, 75)


def assign_columns(): 
    if args['--stats']=='median':
	return ["diff_median","N_Perm","median_pv_diff","median_padjust_diff"]
    elif args['--stats']=="perc25":
	return ["diff_LowQ","N_Perm","lowerq_pv_diff","lowerq_padjust_diff"]
    elif args['--stats']=="perc75":
	return["diff_UpQ","N_Perm","upperq_pv_diff","upper_padjust_diff"]
    elif args['--stats']=='ratio_median':
	return ["ratio_median","N_Perm","median_pv_ratio", "median_padjust_ratio"]
    elif args['--stats']=="ratio_perc25":
	return ["ratio_LowQ","N_Perm","lowerq_pv_ratio", "lower_padjust_ratio"]
    elif args['--stats']=="ratio_perc75":
	return ["ratio_UpQ","N_Perm","upperq_pv_ratio", "upper_pv_ratio", "upper_padjust_ratio"]
    else:
	return None
    
####MAIN#####################################################################################################

col_out=assign_columns()	## Two reasons: Check if stats is correct and create colnames
if col_out==None:
    sys.exit("You have no selected correct statistic")

f=pd.read_csv(args['--file'], sep=" ")

sys.stderr.write('Table read\n')

cd1= f['group']==args['--c1']  ##We have a boolean cond1=condition
cd2=[ not i for i in cd1 ]  ##We have a boolean cond1=condition

out=open(args['--output'],'w')
names=f.columns[4:]
out_df=pd.DataFrame(index=names, columns=col_out)
hits=dict()

out_df[col_out[0]]= [ pvalue(i)  for i in names ]

#print out_df[col_out[0]]
sys.stderr.write('Permutations begin...\n')
for i in names:
    hits[i]=0

for step in range(0,int(args["-n"])):
    if len(hits.keys())==0:
	break
    np.random.shuffle(cd1)
    cd2=[not i for i in cd1]
    for n in hits.keys():
	hits[n]+=compare_pv(n)
	if hits[n]==int(args["-H"]):
	    out_df[col_out[1]][n]=step+1
	    out_df[col_out[2]][n]=hits[n]/float(out_df[col_out[1]][n])
	    del hits[n]
	    

for n in hits.keys():
    out_df[col_out[1]][n]=int(args["-n"])
    out_df[col_out[2]][n]=(hits[n]+1)/float(out_df[col_out[1]][n])
    del hits[n]

sys.stderr.write('Done\n')
out_df[col_out[-1]] = pd.Series(stats.p_adjust(FloatVector(out_df[col_out[2]]), method = 'BH'), index=out_df.index)
out_df.to_csv(args['--output'], sep='\t', index_label="genes")


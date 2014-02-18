#!/usr/bin/env python


"""

PyPermutations has been designed to perform permutations test in situations with big datasets and two groups. Finally it executes a FDR.
Initially has been designed to perform DE.

More information in github.io/marcDabad/PyPermutations

USAGE:

    PyPermutations.py -f <inputFile> --c1 <condition1> --c2 <condition2> [ -n <int> -H <int> -s <statistic> ] -o <outputfile>
    
Inputs:
    -f, --file <inputFile>	Space separated table
    
Output:
    -o, --output <outFile>	Output table
Options:
    --c1	<condition1>	Condition that determine one of two groups
    --c2	<condition2>	Conditions that determine other group
    -n	<int>			Maximum number of permutations [default: 1000000]
    -H	<int>			Maximum hits to stop permutations
    -s, --stats	<statistic>	Statistic to compute pvalue median|perc25|perc75|ratio_median|ratio_perc25|ratio_perc75 [default: median]
    
"""

import numpy as np
import argparse
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')
import pandas as pd
import sys
from docopt import docopt
args=docopt(__doc__)

print args

f=pd.read_csv(args['--file'], sep=" ")
condition=f['group']
condShuffle=f['group'].copy()

def median_diff(array, labels, c1, c2):
    return abs(array[labels==c1].median()-array[labels==c2].median())

def percentil_diff(array, labels, c1, c2, percen):
    return abs(np.percentile(array[labels==c1], percen)-np.percentile(array[labels==c2], percen))

def median_ratio(array, labels, c1, c2):
    return array[labels==c1].median()/array[labels==c2].median()

def percentil_ratio(array, labels, c1, c2, percen):
    return np.percentile(array[labels==c1], percen)/np.percentile(array[labels==c2], percen)

  
    
def perm_med_diff(x , name, h, n , c1='control', c2='case'):
 
    ori=median_diff(x, condition, c1, c2)

    pv=''
    nperm=''
    per=0
    for p in range(n):
	np.random.shuffle(condShuffle)
	md=median_diff(x, condShuffle, c1, c2)
	if md>=ori:
		per+=1
		#sys.stderr.write(name+"\t"+str(per)+"\n")
	if per==h:
	    pv=h/float(p+1)
	    nperm=p+1
	    break
    
    if not pv:
	pv=(per+1)/float(n)
	nperm=n
    #sys.stderr.write(name+"\t"+str(per)+"\t")
    return (ori, pv, nperm)

def perc25_diff(x , name, h, n , c1='control', c2='case'):
 
    ori=percentil_diff(x, condition, c1, c2, 25)

    pv=''
    nperm=''
    per=0
    for p in range(n):
	np.random.shuffle(condShuffle)
	md=percentil_diff(x, condShuffle, c1, c2, 25)
	if md<=ori:
		per+=1

	if per==h:
	    pv=h/float(p+1)
	    nperm=p+1
	    
	    break
    if not pv:
	pv = (per+1)/float(n)
	nperm=n

    return (ori, pv, nperm)


def perc75_diff(x , name, h, n , c1='control', c2='case'):
 
    ori=percentil_diff(x, condition, c1, c2, 75)

    pv=''
    nperm=''
    per=0
    for p in range(n):
	np.random.shuffle(condShuffle)
	md=percentil_diff(x, condShuffle, c1, c2, 75)
	if md>=ori:
		per+=1

	if per==h:
	    pv=h/float(p+1)
	    nperm=p+1
	    
	    break
    if not pv:
	pv=per+1/float(n)
	nperm=n

    return (ori, pv, nperm)


def perm_med_ratio(x , name, h, n , c1='control', c2='case'):
 
    ori=median_ratio(x, condition, c1, c2)

    pv=''
    nperm=''
    per=0
    for p in range(n):
	np.random.shuffle(condShuffle)
	md=median_diff(x, condShuffle, c1, c2)
	if md>1:
	    if md>=ori:
		per+=1
	else:
	    if md<=ori:
		per+=1
	if per==h:
	    pv=h/float(p+1)
	    nperm=p+1
	    break
    if not pv:
	pv=per+1/float(n)
	nperm=n

    return (ori, pv, nperm)

def perm_perc25_ratio(x , name, h, n , c1='control', c2='case'):
 
    ori=percentil_ratio(x, condition, c1, c2, 25)

    pv=''
    nperm=''
    per=0
    for p in range(n):
	np.random.shuffle(condShuffle)
	md=percentil_ratio(x, condShuffle, c1, c2, 25)
	if md>1:
	    if md>=ori:
		per+=1
	else:
	    if md<=ori:
		per+=1

	if per==h:
	    pv=h/float(p+1)
	    nperm=p+1
	    
	    break
    if not pv:
	pv=per+1/float(n)
	nperm=n

    return (ori, pv, nperm)


def perm_perc75_ratio(x , name, h, n , c1='control', c2='case'):
 
    ori=percentil_ratio(x, condition, c1, c2, 75)

    pv=''
    nperm=''
    per=0
    for p in range(n):
	np.random.shuffle(condShuffle)
	md=percentil_ratio(x, condShuffle, c1, c2, 75)
	if md>1:
	    if md>=ori:
		per+=1
	else:
	    if md<=ori:
		per+=1

	if per==h:
	    pv=h/float(p+1)
	    nperm=p+1
	    
	    break
    if not pv:
	pv=per+1/float(n)
	nperm=n

    return (ori, pv, nperm)



out=open(args['--output'],'w')


if args['--stats']=='median':
    out.write("genes\tdiff_median\tN_Perm\tmedian_pv_diff\n")
    for i in f.columns[[ 'ENSG' in x for x in f.columns ]]:
	(medORIdiff, pmeddiff,nperm)=perm_med_diff(f[i], i, h=int(args['-H']), n=int(args['-n']), c1=args['--c1'], c2=args['--c2'])
	out.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n")
	#sys.stderr.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n") 

elif args['--stats']=="perc25":
    out.write("genes\tdiff_LowQ\tN_Perm\tlowerq_pv_diff\n")
    for i in f.columns[[ 'ENSG' in x for x in f.columns ]]:
	(medORIdiff, pmeddiff,nperm)=perc25_diff(f[i], i, h=int(args['-H']), n=int(args['-n']), c1=args['--c1'], c2=args['--c2'])
	out.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n")
	#sys.stderr.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n") 

elif args['--stats']=="perc75":
    out.write("genes\tdiff_UpQ\tN_Perm\tupperq_pv_diff\n")
    for i in f.columns[[ 'ENSG' in x for x in f.columns ]]:
	(medORIdiff, pmeddiff,nperm)=perc75_diff(f[i], i, h=int(args['-H']), n=int(args['-n']), c1=args['--c1'], c2=args['--c2'])
	out.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n")
	#sys.stderr.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n") 

elif args['--stats']=='ratio_median':
    out.write("genes\tratio_median\tN_Perm\tmedian_pv_ratio\n")
    for i in f.columns[[ 'ENSG' in x for x in f.columns ]]:
	(medORIdiff, pmeddiff,nperm)=perm_med_diff(f[i], i, h=int(args['-H']), n=int(args['-n']), c1=args['--c1'], c2=args['--c2'])
	out.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n")
	#sys.stderr.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n") 

elif args['--stats']=="ratio_perc25":
    out.write("genes\tratio_LowQ\tN_Perm\tlowerq_pv_ratio\n")
    for i in f.columns[[ 'ENSG' in x for x in f.columns ]]:
	(medORIdiff, pmeddiff,nperm)=perc25_diff(f[i], i, h=int(args['-H']), n=int(args['-n']), c1=args['--c1'], c2=args['--c2'])
	out.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n")
	#sys.stderr.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n") 

elif args['--stats']=="ratio_perc75":
    out.write("genes\tratio_UpQ\tN_Perm\tupperq_pv_ratio\n")
    for i in f.columns[[ 'ENSG' in x for x in f.columns ]]:
	(medORIdiff, pmeddiff,nperm)=perc75_diff(f[i], i, h=int(args['-H']), n=int(args['-n']), c1=args['--c1'], c2=args['--c2'])
	out.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n")
	#sys.stderr.write(i+"\t"+str(medORIdiff)+"\t"+str(nperm)+"\t"+str("%f" % pmeddiff)+"\n") 

else:
    print " No has seleccionat el statistic correcte!!!!!!"
    
    
  
out.close()


data=pd.read_csv(args['--output'], sep="\t")

if args['--stats']=='median':
    data['median_padjust_diff'] = pd.Series(stats.p_adjust(FloatVector(data['median_pv_diff']), method = 'BH'), index=data.index)
elif args['--stats']=='perc25':
    data['lowerq_padjust_diff'] = pd.Series(stats.p_adjust(FloatVector(data['lowerq_pv_diff']), method = 'BH'), index=data.index)
elif args['--stats']=='perc75':
    data['upperq_padjust_diff'] = pd.Series(stats.p_adjust(FloatVector(data['upperq_pv_diff']), method = 'BH'), index=data.index)
elif args['--stats']=='ratio_median':
    data['median_padjust_ratio'] = pd.Series(stats.p_adjust(FloatVector(data['median_pv_ratio']), method = 'BH'), index=data.index)
elif args['--stats']=='perc25':
    data['lowerq_padjust_ratio'] = pd.Series(stats.p_adjust(FloatVector(data['lowerq_pv_ratio']), method = 'BH'), index=data.index)
elif args['--stats']=='perc75':
    data['upperq_padjust_ratio'] = pd.Series(stats.p_adjust(FloatVector(data['upperq_pv_ratio']), method = 'BH'), index=data.index)

else:
    print "NO has seleccionat el statistic correcte!!!!!"
    
    
data.to_csv(args['--output'], sep='\t', index=False)

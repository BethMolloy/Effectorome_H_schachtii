import os
import json
import glob, os
import statistics
import pandas as pd

os.chdir('/dir/of/pdb/files')

threshold = 0.5 # adjust me for strictness of pdb quality. 0.5 means ptm limit is 0.5 and plddt limit is 50.


with open('ptmValues.csv','w') as outfile:
        for file in glob.glob("*000.json"):
                with open(file, 'r') as f:
                        loadJSON = json.load(f)
                        ptmValue = loadJSON['ptm']
                        plddt = statistics.mean(loadJSON['plddt'])
                        #print(plddt)
                        outfile.write(file+','+str(ptmValue)+','+str(plddt)+'\n')
'''
!!! ### warnig below removes any PDB files that are below treshold ### !!!
!!! ### warnig below removes any PDB files that are below treshold ### !!!
!!! ### warnig below removes any PDB files that are below treshold ### !!!
'''


df = pd.read_csv('ptmValues.csv',header=None)



largerThanTresholdDF = df[df[1] > threshold] 
lowerThanTresholdDF = df[(df[1] <= threshold) | (df[2] <= threshold*100)]

for i in range(0,len(lowerThanTresholdDF)):
        searchString = (lowerThanTresholdDF[0].iloc[i].split('_scores_')[0]+'*')
        for f in glob.glob(searchString):
                if ('pdb' in f) or ('json' in f):
                        #os.remove(f) #uncomment me to actually delete
                        print(f)

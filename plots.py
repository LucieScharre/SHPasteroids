#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 16:39:29 2021

@author: luciescharre
"""

import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import math
import seaborn as sns
sns.set()
sns.set_style('whitegrid')


# read in hits be tokens[0]
# read palermo for those
# can use counter to count number of hits per asteroid
# can reduce to individual names ESA_MPC_old_Sentry = list(dict.fromkeys(object_ids))
plt.rcParams["font.family"] = "serif"

###################################### Histogram of all Sentry
object_hits = []

sentry_hits = []
palermo_hits = []

Sentry_files = ['POSS-I Sentry bright.txt','POSS-II Sentry bright.txt',
                    'ESO Sentry bright.txt',
                     'UKST all Sentry.txt' ]

Sentry_files = [
                'POSS-II Sentry.txt',
                'ESO Sentry.txt',
                'UKST all Sentry.txt'  ]

#'POSS-I Sentry.txt'

for filename in Sentry_files:
    with open(filename, 'r') as f:
        for line in f:
            tokens = line.split(',')
            sentry_hits.append(str(tokens[0]).rstrip())
            object_hits.append(str(tokens[0]).rstrip())
        
 

object_ids = []
palermo = [] 


with open('sentry 2.txt', 'r') as f:
    for line in f:
        tokens = line.split('   ')
        object_ids.append(str(tokens[0]).rstrip())
        palermo.append(float(tokens[1]))           
        
for hit in sentry_hits:
    palermo_hits.append(palermo[object_ids.index(hit)])
    
bins = np.arange(-13,-1,0.5)
ylim = [0,45]
xlim =[-13,-1]

plt.hist(palermo_hits, bins = bins)
plt.xlim(xlim)
plt.ylim(ylim)
plt.title('All hits from the NASA Sentry List')
#plt.minorticks_on()
#plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
#plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
plt.xlabel('Palermo Rating')
plt.ylabel('Number of Hits')
plt.show()
     
#print(len(palermo_hits))



###################################### Histogram of all Sentry & ESA and just ESA
  
#'POSS-I ESA.txt',
ESA_files = ['UKST all ESA.txt','POSS-II ESA.txt']   

ESA_hits = []

for filename in ESA_files:
    with open(filename, 'r') as f:
        for line in f:
            tokens = line.split(',')
            ESA_hits.append(str(tokens[0]).rstrip())
            object_hits.append(str(tokens[0]).rstrip())
            
"""            
with open('ESA.txt', 'r') as f:
    for line in f:
        object_ids.append(str(line[0:4]).rstrip()+' '+str(line[4:9]).rstrip())
        palermo.append(float(line[76:84]))"""

object_ids = []
palermo = [] 
ESA_palermo_hits = []

with open('ESA.txt', 'r') as f:
    for line in f:

        object_ids.append(str(line[0:4]).rstrip()+' '+str(line[4:9]).rstrip())
        palermo.append(float(line[75:83]))

        
for hit in ESA_hits:
    palermo_hits.append(palermo[object_ids.index(hit)])  
    ESA_palermo_hits.append(palermo[object_ids.index(hit)])  
    

print(len(palermo_hits))
plt.hist(palermo_hits, bins =bins)
plt.xlim(xlim)
plt.ylim(ylim)
#plt.minorticks_on()
#plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
#plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
#plt.title('All hits from NASA Sentry & ESA list')
plt.xlabel('Palermo Rating')
plt.ylabel('Number of Hits')
plt.show()

plt.hist(ESA_palermo_hits, bins = bins)
plt.xlim(xlim)
plt.ylim(ylim)
#plt.minorticks_on()
#plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
#plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
plt.title('All hits from ESA List')
plt.xlabel('Palermo Rating')
plt.ylabel('Number of Hits')
plt.show()
        
###################################### Histogram of Sentry&ESA in UKST


object_hits = []

sentry_hits = []
palermo_hits = []

Sentry_files = ['UKST all Sentry.txt' ]

for filename in Sentry_files:
    with open(filename, 'r') as f:
        for line in f:
            tokens = line.split(',')
            sentry_hits.append(str(tokens[0]).rstrip())
            object_hits.append(str(tokens[0]).rstrip())
        
 

object_ids = []
palermo = [] 


with open('sentry 2.txt', 'r') as f:
    for line in f:
        tokens = line.split('   ')
        object_ids.append(str(tokens[0]).rstrip())
        palermo.append(float(tokens[1]))           
        
for hit in sentry_hits:
    palermo_hits.append(palermo[object_ids.index(hit)])
    
    

object_ids = []
palermo = [] 
ESA_palermo_hits = []

with open('ESA.txt', 'r') as f:
    for line in f:

        object_ids.append(str(line[0:4]).rstrip()+' '+str(line[4:9]).rstrip())
        palermo.append(float(line[75:83]))

        
for hit in ESA_hits:
    palermo_hits.append(palermo[object_ids.index(hit)])  
    ESA_palermo_hits.append(palermo[object_ids.index(hit)])  
    
  
plt.hist(palermo_hits, bins = bins)
plt.xlim(xlim)
plt.ylim(ylim)
#plt.minorticks_on()
#plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
#plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
#plt.title('Hits in UKST catalog from Sentry NASA & ESA List')
plt.xlabel('Palermo Rating')
plt.ylabel('Number of Hits')
plt.show()
     
print(len(palermo_hits))
  
 
###################################### Histogram of low uncertainty Sentry in UKST (none in ESA)

object_hits = []

sentry_hits = []
palermo_hits = []

Sentry_files = ['UKST all Sentry low del.txt' ]

for filename in Sentry_files:
    with open(filename, 'r') as f:
        for line in f:
            tokens = line.split(',')
            sentry_hits.append(str(tokens[0]).rstrip())
            object_hits.append(str(tokens[0]).rstrip())
        
 

object_ids = []
palermo = [] 


with open('sentry 2.txt', 'r') as f:
    for line in f:
        tokens = line.split('   ')
        object_ids.append(str(tokens[0]).rstrip())
        palermo.append(float(tokens[1]))           
        
for hit in sentry_hits:
    palermo_hits.append(palermo[object_ids.index(hit)])
    
  
plt.hist(palermo_hits, bins = bins)
plt.xlim(xlim)
plt.ylim(ylim)
#plt.minorticks_on()
#plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
#plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
#plt.title('Sentry: Number of hits in UKST catalog with low del vs Palermo')
plt.title('Hits in UKST catalog from Sentry NASA & ESA List with low error')
plt.xlabel('Palermo Rating')
plt.ylabel('Number of Hits')
plt.show()
     
print(len(palermo_hits))
    
################################ Scatter plot uncertainties    

# scatter plot in degrees? presented measurements in RA and DEC
# in degrees ra and dec both the same, maybe combine with histogram and plt versus palermo, can color code for catalog or list
# collect hits from txt documents and recompute the uncertainty, collect that in degrees


delta_files = ['UKST Survey deltas.txt','UKST all deltas.txt','POSS-II deltas.txt','ESO deltas.txt']

objects = []
delta = []

delta_UKST = []
delta_UKST_S = []
delta_ESO = []
delta_POSS_II = []

palermos_POSS_II = []
palermos_UKST = []
palermos_UKST_S = []
palermos_ESO = []

catalogue = []
palermos = []

object_ids = []
palermo = [] 


with open('sentry 2.txt', 'r') as f:
    for line in f:
        tokens = line.split('   ')
        object_ids.append(str(tokens[0]).rstrip())
        palermo.append(float(tokens[1]))      
        

with open('ESA.txt', 'r') as f:
    for line in f:
        if str(line[0:4]).rstrip()+' '+str(line[4:9]).rstrip() not in object_ids:
            object_ids.append(str(line[0:4]).rstrip()+' '+str(line[4:9]).rstrip())
            palermo.append(float(line[75:83]))

for filename in delta_files:
    with open(filename, 'r') as f:
        for line in f:
            
            tokens = line.split(',')
            sentry_hits.append(str(tokens[0]).rstrip())
            
            if str(tokens[2]).rstrip()==' UKST':
                delta_UKST.append(float(tokens[1]))
                palermos_UKST.append(palermo[object_ids.index(str(tokens[0]).rstrip())])
                
            if str(tokens[2]).rstrip()==' POSS-II'  :
                delta_POSS_II.append(float(tokens[1]))
                palermos_POSS_II.append(palermo[object_ids.index(str(tokens[0]).rstrip())])
            
            if str(tokens[2]).rstrip()==' ESO'  : 
                delta_ESO.append(float(tokens[1]))
                palermos_ESO.append(palermo[object_ids.index(str(tokens[0]).rstrip())])
                
            if str(tokens[2]).rstrip()==' UKST Survey'  : 
                delta_UKST_S.append(float(tokens[1]))
                palermos_UKST_S.append(palermo[object_ids.index(str(tokens[0]).rstrip())])
                
            delta.append(float(tokens[1]))
            catalogue.append(str(tokens[2]).rstrip())
            palermos.append(palermo[object_ids.index(str(tokens[0]).rstrip())])


plt.scatter(palermos,np.log(delta))

plt.hlines(math.log(3),xlim[0],xlim[1])
plt.hlines(math.log(360),xlim[0],xlim[1],color='r')
plt.xlim(xlim)
plt.title('Uncertainty of hits in Sentry NASA & ESA List')
plt.xlabel('Palermo Rating')
plt.ylabel('Uncertainty in log(Degrees)')
plt.minorticks_on()
plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
plt.show()


#plt.title('Uncertainty of hits in Sentry NASA & ESA List')
plt.scatter(palermos_UKST,np.log(delta_UKST),label ='UKST')
plt.scatter(palermos_UKST_S,np.log(delta_UKST_S),label ='UKST Survey', color='k')
plt.scatter(palermos_ESO,np.log(delta_ESO),label ='ESO')
plt.scatter(palermos_POSS_II,np.log(delta_POSS_II),label ='POSS-II')
plt.scatter(-2.79,np.log(0.0019633333333333334),label='1979 XB hit', color = 'r')

plt.hlines(math.log(3),xlim[0],xlim[1])
plt.hlines(math.log(360),xlim[0],xlim[1],color='r')
plt.xlim(xlim)
#plt.grid(b=True, which='major', color='#666666', linestyle='-')
#plt.minorticks_on()
#plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

plt.legend(loc='lower left')
plt.xlabel('Palermo Rating')
plt.ylabel('Uncertainty in log(Degrees)')
#plt.minorticks_on()
#plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
#plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
plt.show()
    

# make another scatter plot of top 10 most dangerous, either same graph with colour or just the ones
# highlight prime plates to check and name them
# maybe rerun to include magnitude and plate  
# distinguish between survey plates and non survey
# make sure to mention that those plates have been checked


objects = []
delta = []

delta_UKST = []
delta_UKST_S = []
delta_ESO = []
delta_POSS_II = []

palermos_POSS_II = []
palermos_UKST = []
palermos_UKST_S = []
palermos_ESO = []

catalogue = []
palermos = []


for filename in delta_files:
    with open(filename, 'r') as f:
        for line in f:
            
                #print(float(tokens[1]))
                tokens = line.split(',')
                if float(tokens[1])<3.0:
                    sentry_hits.append(str(tokens[0]).rstrip())
                    
                    if str(tokens[2]).rstrip()==' UKST':
                        delta_UKST.append(float(tokens[1]))
                        palermos_UKST.append(palermo[object_ids.index(str(tokens[0]).rstrip())])
                        
                    if str(tokens[2]).rstrip()==' POSS-II'  :
                        delta_POSS_II.append(float(tokens[1]))
                        palermos_POSS_II.append(palermo[object_ids.index(str(tokens[0]).rstrip())])
                    
                    if str(tokens[2]).rstrip()==' ESO'  : 
                        delta_ESO.append(float(tokens[1]))
                        palermos_ESO.append(palermo[object_ids.index(str(tokens[0]).rstrip())])
                        
                    if str(tokens[2]).rstrip()==' UKST Survey'  : 
                        delta_UKST_S.append(float(tokens[1]))
                        palermos_UKST_S.append(palermo[object_ids.index(str(tokens[0]).rstrip())])
                        
                    delta.append(float(tokens[1]))
                    catalogue.append(str(tokens[2]).rstrip())
                    palermos.append(palermo[object_ids.index(str(tokens[0]).rstrip())])
            
                
plt.title('Uncertainty of hits in Sentry NASA & ESA List')
plt.scatter(palermos_UKST,delta_UKST,label ='UKST')
plt.scatter(palermos_UKST_S,delta_UKST_S,label ='UKST Survey', color='k')
plt.scatter(palermos_ESO,delta_ESO,label ='ESO')
plt.scatter(palermos_POSS_II,delta_POSS_II,label ='POSS-II')
plt.scatter(-2.79,0.0019633333333333334,label='1979 XB hit', color = 'r')
plt.legend(loc='upper right')
plt.xlabel('Palermo Rating')
plt.ylabel('Uncertainty in Degrees')
#plt.grid(b=True, which='major', color='#666666', linestyle='-')
#plt.minorticks_on()
#plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
plt.show()
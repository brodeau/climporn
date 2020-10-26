#!/usr/bin/env python3
#
#     CLIMPORN
#
#    L. Brodeau, November 2020

import sys
import numpy as nmp
from matplotlib import cm

vcm  = cm.bone(range(256))


(nl,nc) = nmp.shape(vcm)

print('\n *** Shape of vcm =',nl,nc)

if nl != 256 or nc != 4:
    print("Problem of shape")
    sys.exit(0)

irgb = nmp.zeros((3,),dtype=int)

f = open('bone.ncmap', 'w')

for jl in range(nl):

    for jc in range(3):
        irgb[jc] = int(round(255.*vcm[jl,jc],0))

    #print(irgb[:])
        
    cline = str(irgb[0])+' '+ str(irgb[1])+' '+ str(irgb[2])+'\n'

    
    
    f.write(cline)

f.close()

#vcm
#print(vcm)

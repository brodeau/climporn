#!/usr/bin/env python3
#
#     CLIMPORN
#
#  Prepare 2D maps (monthly) that will later become a movie!
#  NEMO output and observations needed
#
#    L. Brodeau, January 2019
#
import sys
from os import path, getcwd, mkdir
#import argparse as ap
import numpy as nmp
from PIL import Image

#from netCDF4 import Dataset


#import clprn_colmap as bcm
#import clprn_tool as bt
import clprn_ncio as bnc

# ClimPorn:
import nemo_hboxes as nhb

# (0-255)     R   G   B  transparency
#rgb_f = [ 0 , 0 , 0 , 255 ] ; # solid black
rgb_f = [ 255, 237, 0 , 255 ] ; # yellow ON (255,237,0)
npf = 2 ; # number of poins for frame...



nratio  = 5

nrim_chld = 1
#nrim_chld = 3
#nrim_chld = 10

ip = 35 ; # location (Fortran indexing) of bottom left corner AGRIF_FixedGrids.in)
jp = 43 ; # of the child box into parent box in Agrif ! (as in 

cdate = '2016-12-22_00'

cf_chld = './figs/CURLOF/CURLOF_CALEDO10-_ALL_'+cdate+'_on2.png'
cf_prnt = './figs/CURLOF/CURLOF_TROPICO05_NST-_ALL_'+cdate+'_on2.png'

chld = Image.open(cf_chld)
prnt = Image.open(cf_prnt)

(nyc,nxc,nrgb_c) = nmp.shape(chld)
(nyp,nxp,nrgb_p) = nmp.shape(prnt)


print(' *** Image array shape for child  = ', nyc,nxc,nrgb_c)
print(' *** Image array shape for parent = ', nyp,nxp,nrgb_p)

#if len(vshape_pic) == 3:    
if nrgb_c != 4: print(' Problem #1 with your child image, not what we expected!') ;
if nrgb_p != 4: print(' Problem #1 with your parent image, not what we expected!') ;

if nxp%nratio!=0 or nyp%nratio!=0: print(' Problem #1 with your parent image, it shoud be a multiple of '+str(nratio)+'!') ;



xchld = nmp.array(chld)
xprnt = nmp.array(prnt)


print("xchld", nmp.shape(xchld))
print("xprnt", nmp.shape(xprnt))

print(xprnt[0,0,:])
print(xprnt[100,100,:])
#sys.exit(0)

ip = ip * nratio - 1
jp = jp * nratio - 1


# Child array we keep after croping:
XC = xchld[nrim_chld:-nrim_chld,nrim_chld:-nrim_chld]
(Ny,Nx,Nr) = nmp.shape(XC)

print(" shape XC =", Ny,Nx,Nr)


# Drawing black frame:
#    ex: blue [  0 137 183 255]
#for jj in [0,Ny-1]:
#    for i in range(4): XC[jj,:,i] = rgb_f[i]
#for ji in [0,Nx-1]:
#    for i in range(4): XC[:,ji,i] = rgb_f[i]

for i in range(4):
    XC[     0:npf,:,i] = rgb_f[i]
    XC[Ny-npf:Ny ,:,i] = rgb_f[i]

    XC[:,     0:npf,i] = rgb_f[i]
    XC[:,Nx-npf:Nx ,i] = rgb_f[i]


    
# Fitting child array into parent array:
#j1 = nyp-jp-nyc
#xprnt[j1:j1+nyc,ip:ip+nxc,:] = xchld[:,:,:]


j1 = nyp-jp-Ny - 2*nrim_chld
i1 = ip        + nrim_chld 
xprnt[j1:j1+Ny,i1:i1+Nx,:] = XC[:,:,:]




#ifield8 = xchld.astype(nmp.uint8)
#image_chld = Image.fromarray(nmp.flipud(ifield8))

image_chld = Image.fromarray(xchld[nrim_chld:-nrim_chld,nrim_chld:-nrim_chld])

image_all  = Image.fromarray(xprnt)





cf_chld = './child.png'
# Then save it:
image_chld.save(cf_chld)
print(' *** Image '+cf_chld+' saved!\n')



cf_all = './parent+child.png'
# Then save it:
image_all.save(cf_all)
print(' *** Image '+cf_all+' saved!\n')



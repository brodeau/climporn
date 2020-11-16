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


nrim_nest = 3

cf_nst = './figs/CURLOF/CURLOF_CALEDO10-_ALL_2016-12-22_00_on2.png'

nest = Image.open(cf_nst)


(ny,nx,nrgb) = nmp.shape(nest)


print(' *** Image array shape = ', ny,nx,nrgb)

#if len(vshape_pic) == 3:    
if nrgb != 4: print(' Problem #1 with your image, not what we expected!') ;

print(" *** shape of nest: ", (ny,nx))

xnest = nmp.array(nest)    


#ifield8 = xnest.astype(nmp.uint8)
#image_out = Image.fromarray(nmp.flipud(ifield8))

image_out = Image.fromarray(xnest[nrim_nest:-nrim_nest,nrim_nest:-nrim_nest])


cf_im = './boo2.png'
# Then save it:
image_out.save(cf_im)
print(' *** Image '+cf_im+' saved!\n')

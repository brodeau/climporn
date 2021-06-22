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
import glob
#import argparse as ap
import numpy as nmp
from PIL import Image
#
# ClimPorn:
from   climporn.utils import chck4f
import climporn.nemo_hboxes as nhb

#nghost = 3

# About files to use:
fpref_p = './figs/CURLOF/CURLOF_TROPICO05_NST-_ALL_'
fsuff_p = '_on2.png'

fpref_c = './figs/CURLOF/CURLOF_CALEDO10-_ALL_'
fsuff_c = '_on2.png'

fpref_o = './figs/montage_CURLOF_TROPICO05-CALEDO10_'

nratio  = 5

ncrop_chld = 3
#ncrop_chld = 1
#ncrop_chld = 3
#ncrop_chld = 10

ipc = 35 ; # location (Fortran indexing) of bottom left corner AGRIF_FixedGrids.in)
jpc = 43 ; # of the child box into parent box in Agrif ! (as in 

cdate = '2016-12-22_00'

# Final box to keep (mind that for images it's flipped upside down!):
j1b = 336 ; i1b = 0
j2b = 680 ; i2b = int(16./9.*float(j2b-j1b))

# COLOR FOR THE NESTING BOX LINE:
# (0-255)     R   G   B  transparency
#rgb_f = [ 0 , 0 , 0 , 255 ] ; # solid black
rgb_f = [ 255, 237, 0 , 255 ] ; # yellow ON (255,237,0)
npf = 2 ; # number of poins for frame...






list_im_p = glob.glob(fpref_p+"*"+fsuff_p)
list_im_c = glob.glob(fpref_c+"*"+fsuff_c)

nbf = len(list_im_c)
if len(list_im_p) != nbf:
    print(' Problem different number of images between parent and child!', len(list_im_p), nbf) ; sys.exit(0)
    



for j in range(nbf):

    cf_prnt = list_im_p[j]
    chck4f(cf_prnt)
    
    cdate = str.replace( cf_prnt, fpref_p, '')
    cdate = str.replace( cdate, fsuff_p, '')
    
    cf_chld =  fpref_c+cdate+fsuff_c ; # that's what we expect
    chck4f(cf_chld)

    cf_out = str.replace( cf_prnt, fpref_p, fpref_o)
    cf_out = str.replace( cf_out, fsuff_p, '.png')


    #print(cf_prnt,' ',cf_chld,' ', cdate, )
    #print(cf_out)    
    #sys.exit(0)




    chld = Image.open(cf_chld)
    prnt = Image.open(cf_prnt)

    (nyc,nxc,nrgb_c) = nmp.shape(chld)
    (nyp,nxp,nrgb_p) = nmp.shape(prnt)


    #print(' *** Image array shape for child  = ', nyc,nxc,nrgb_c)
    #print(' *** Image array shape for parent = ', nyp,nxp,nrgb_p)

    if nrgb_c != 4: print(' Problem #1 with your child image, not what we expected!') ; sys.exit(0)
    if nrgb_p != 4: print(' Problem #1 with your parent image, not what we expected!') ; sys.exit(0)

    if nxp%nratio!=0 or nyp%nratio!=0: print(' Problem #1 with your parent image, it shoud be a multiple of '+str(nratio)+'!') ; sys.exit(0)



    xchld = nmp.array(chld)
    xprnt = nmp.array(prnt)

    chld.close()
    prnt.close()

    print("xchld", nmp.shape(xchld))
    print("xprnt", nmp.shape(xprnt))
    #print(xprnt[0,0,:])
    #print(xprnt[100,100,:])


    #ip = (ipc - 1) * nratio + nghost
    #jp = (jpc - 1) * nratio + nghost
    ip = (ipc + 1) * nratio + 2  # I have no fucking clue why ????
    jp = (jpc + 1) * nratio + 2  # I have no fucking clue why ????
    #ip = ipc * nratio
    #jp = jpc * nratio


    # Child array we keep after croping:
    if ncrop_chld >= 1:
        XC = xchld[ncrop_chld:-ncrop_chld,ncrop_chld:-ncrop_chld]
    else:
        XC = xchld[:,:]
    (Ny,Nx,Nr) = nmp.shape(XC)

    #print(" shape XC =", Ny,Nx,Nr)


    # Drawing  frame:
    for i in range(4):
        XC[     0:npf,:,i] = rgb_f[i]
        XC[Ny-npf:Ny ,:,i] = rgb_f[i]        
        XC[:,     0:npf,i] = rgb_f[i]
        XC[:,Nx-npf:Nx ,i] = rgb_f[i]

    j1 = nyp-jp-Ny - ncrop_chld
    i1 = ip        + ncrop_chld 
    xprnt[j1:j1+Ny,i1:i1+Nx,:] = XC[:,:,:]

    image_out  = Image.fromarray(xprnt[j1b:j2b,i1b:i2b,:])

    #cf_chld = './child.png'
    ## Then save it:
    #image_chld.save(cf_chld)
    #print(' *** Image '+cf_chld+' saved!\n')




    # Then save it:
    image_out.save(cf_out)
    print(' *** Image '+cf_out+' saved!\n')

    
    
    del XC, xprnt, xchld

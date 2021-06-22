#!/usr/bin/env python3

#       B a r a K u d a
#
#       L. Brodeau, 2017]

import sys
import numpy as nmp
from PIL import Image
import os
from netCDF4 import Dataset
import datetime

#l_fake_coor = True
l_fake_coor = False

l_nemo_like = True

narg = len(sys.argv)
if not narg in [2, 3]:
    print('Usage: '+sys.argv[0]+' <image> (<field divider for field>)'); sys.exit(0)

cf_im = sys.argv[1]

idiv = 1
if narg == 3: idiv = int(sys.argv[2])

print(idiv)

cfname, cfext = os.path.splitext(cf_im)


#(nj,ni) = nmp.shape(nav_lon)

cf_nc = str.replace(os.path.basename(cf_im), cfext, '.nc')

# Opening Images:
print(' *** Opening image '+cf_im)
pic = Image.open(cf_im)




lcolor = False ; # if false it is a 

vshape_pic = nmp.shape(pic)

if len(vshape_pic) == 3:
    (ny,nx,nrgb) = vshape_pic
    if nrgb != 3: print(' Problem #1 with your image, not what we expected!') ; sys.exit(0)
    lcolor = True    ;  # RGB color picture => 3 2D array
    print("\n It's a RGB color picture!\n")
    
elif len(vshape_pic) == 2:
    lcolor = False   ;  # grey-scale picture (true black and white) => 1 2D array
    (ny,nx) = vshape_pic
    nrgb = 1
    print("\n It's a grey-scale B&W picture!\n")
else:
    print(' Problem #2 with your image, not what we expected!') ; sys.exit(0)




print(" *** shape of pic: ", (ny,nx))

xpic = nmp.array(pic)


if l_fake_coor:
    # Prepare coordinates if needed:
    vlon = nmp.zeros(nx) ; dx = 360./float(nx)
    for ji in range(nx): vlon[ji] = (float(ji) + 0.5)*dx
    
    vlat = nmp.zeros(ny) ; dy = 180./float(ny)
    for jj in range(ny): vlat[jj] = -90 + (float(jj) + 0.5)*dy
    #print(vlat[:])
    #sys.exit(0)


f_out = Dataset(cf_nc, 'w', format='NETCDF4')

# Dimensions:

cdim_x = 'longitude'
cdim_y = 'latitude'

if l_nemo_like:
    cdim_x = 'x'
    cdim_y = 'y'


#if l_fake_coor:
#    cdim_x = 'lon'
#    cdim_y = 'lat'


f_out.createDimension(cdim_x, nx)
f_out.createDimension(cdim_y, ny)

#if l_nemo_like: f_out.createDimension('t', None)

if l_fake_coor:
    id_lon  = f_out.createVariable('lon0','f4',(cdim_x,))
    id_lat  = f_out.createVariable('lat0','f4',(cdim_y,))
    id_lon[:] = vlon[:]
    id_lat[:] = vlat[:]



if lcolor:
    
    id_red  = f_out.createVariable('red','f4',(cdim_y,cdim_x,))
    id_red.long_name = 'Red (of RGB)'

    id_green  = f_out.createVariable('green','f4',(cdim_y,cdim_x,))
    id_green.long_name = 'Green (of RGB)'

    id_blue  = f_out.createVariable('blue','f4',(cdim_y,cdim_x,))
    id_blue.long_name = 'Blue (of RGB)'

    id_red[:,:]   = nmp.flipud(xpic[:,:,0])
    id_green[:,:] = nmp.flipud(xpic[:,:,1])
    id_blue[:,:]  = nmp.flipud(xpic[:,:,2])

else:
    #if l_nemo_like:
    #    id_bw  = f_out.createVariable('bw','i1',('t',cdim_y,cdim_x,))
    #    id_bw.long_name = 'Grey scale'
    #    #id_bw[0,:,:]   = nmp.flipud(xpic[:,:]) / idiv
    #    id_bw[0,:,:]   = 1 - (nmp.flipud(xpic[:,:]) + 1)/idiv
    #else:
    id_bw  = f_out.createVariable('bw','i1',(cdim_y,cdim_x,))
    id_bw.long_name = 'Grey scale'
    id_bw[:,:]   = 1 - (nmp.flipud(xpic[:,:]) + 1)/idiv


f_out.About  = 'Image '+cf_im+' converted to netcdf.'
f_out.Author = 'Generated with image_to_netcdf.py of BARAKUDA (https://github.com/brodeau/barakuda)'

f_out.close()



print(cf_nc+' created!!!')


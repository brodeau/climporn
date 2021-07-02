#!/usr/bin/env python3
#
#     CLIMPORN
#
# Converts a NetCDF 2D mask field into a black & white bitmap image
#  ==> idea is to edit this image with a powerful image-editing software like Gimp
#      and convert back the result to a new NetCDF field using `imageBW_to_NetCDF.py` !!!
#
#    L. Brodeau, 2021
#

import sys
import numpy as nmp
from PIL import Image
import os
from netCDF4 import Dataset

l_fake_coor = True
#l_fake_coor = False

narg = len(sys.argv)
if narg not in [4]:
    print('Usage: '+sys.argv[0]+' <netcdf_file.nc> <netcdf_variable> <image_extension (jpg,png,bmp,...)>'); sys.exit(0)
#if narg not in [4, 5]:
#    print('Usage: '+sys.argv[0]+' <netcdf_file.nc> <netcdf_variable> <image_extension (jpg,png,bmp,...)> (mutiple to field)'); sys.exit(0)

cf_nc = sys.argv[1]
cv_nc = sys.argv[2]
ciext = sys.argv[3]

imult = 255
#imult = 1
#if narg == 5: imult = int(sys.argv[4])


print(imult)
    
cfname, cncext = os.path.splitext(cf_nc)


cf_im = str.replace(os.path.basename(cf_nc), cncext, '.'+ciext)

print(' *** Will create image '+cf_im)






# Reading data array:
f_nc = Dataset(cf_nc)
Ndim = len(f_nc.variables[cv_nc].dimensions)
if   Ndim == 4:
    xfield = imult*f_nc.variables[cv_nc][0,0,:,:]
elif Ndim == 3:
    xfield = imult*f_nc.variables[cv_nc][0,:,:]
elif Ndim == 2:
    xfield = imult*f_nc.variables[cv_nc][:,:]
else:
    print(' ERROR (mk_zonal_average.py) => weird shape for your mask array!')
    sys.exit(0)
#xfield  = imult*f_nc.variables[cv_nc][:,:]
f_nc.close()



(ny,nx) = nmp.shape(xfield)

ifield = nmp.zeros((ny,nx), dtype=nmp.int16)

xfield[:,:] = nmp.round(xfield[:,:], 0)

ifield = xfield.astype(nmp.int16)

# Cleaning overshoots:
idx_too_small = nmp.where(ifield < 0)
ifield[idx_too_small] = 0
idx_too_large = nmp.where(ifield > 255)
ifield[idx_too_large] = 255

#print(ifield[:,22])

ifield8 = ifield.astype(nmp.uint8)


image = Image.fromarray(nmp.flipud(ifield8))

# Then save it:
image.save(cf_im)
print(' *** Image '+cf_im+' saved!\n')

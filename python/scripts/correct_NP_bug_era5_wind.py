#!/usr/bin/env python3
#
#     CLIMPORN
#
#  Generate a netCDF file containing a land-sea mask (water==1, land==0),
#  based on the specified field (out of value or `_FillValue` attribute)
#  of the specified input file
#
#    L. Brodeau, 2017
###############################################################################

import sys
import numpy as nmp
from os import path

from netCDF4 import Dataset

from climporn import utils as cpu

#l_fake_coor = True
#l_fake_coor = False

#l_use_fillval = True

#cv_coefr  = 'socoefr'
#cv_runoff = 'sorunoff'

#cv_lon = 'longitude'
#cv_lat = 'latitude'
cv_lon = 'lon'
cv_lat = 'lat'
cv_tim = 'time'

narg = len(sys.argv)
if narg not in [4]:
    print('Usage: '+sys.argv[0]+' <file.nc> <var> <u/v>')
    sys.exit(0)

cf_in = sys.argv[1]
cv_in = sys.argv[2]
cvect = sys.argv[3]
    
cf_out = 'NEW_'+cf_in
 
print(' *** Will create new file: '+cf_out)


f_in = Dataset(cf_in)
Ni   = f_in.dimensions[cv_lon].size
Nj   = f_in.dimensions[cv_lat].size
Nt   = f_in.dimensions[cv_tim].size

print("Shape of domain: Ni, Nj =", Ni, Nj)
print("Number of time records =", Nt)
#sys.exit(0)

vlon = nmp.zeros(Ni)
vlat = nmp.zeros(Nj)

vlon[:]   = f_in.variables[cv_lon][:]
vlat[:]   = f_in.variables[cv_lat][:]


xin = nmp.zeros((Nj,Ni))
###############################################################

f_out = Dataset(cf_out, 'w', format='NETCDF4')

# Dimensions:
f_out.createDimension(cv_tim, None)
f_out.createDimension(cv_lat, Nj)
f_out.createDimension(cv_lon, Ni)

id_tim  = f_out.createVariable(cv_tim, 'f4',(cv_tim,))
id_lat  = f_out.createVariable(cv_lat, 'f4',(cv_lat,))
id_lon  = f_out.createVariable(cv_lon, 'f4',(cv_lon,))

id_lat[:] = vlat[:]
id_lon[:] = vlon[:]

id_out      = f_out.createVariable(cv_in, 'f4',(cv_tim,cv_lat,cv_lon,), zlib=True, complevel=8)


# Need to find oposite longitude:
#  exp: 90 => 270
vidx_mirror = nmp.zeros(Ni,dtype=int)
for ji in range(Ni):
    rl = vlon[ji]%360.
    print('\n rl = ', rl)
    rt = rl + 180.
    if rt > 360.: rt = rt - 360.
    if not rt in [180.,360.]:
        print(' rt = ', rt)
        [[ii]] = nmp.argwhere(vlon == rt )
        print(" --- lon = ",rl, " ==> mirror =", vlon[ii], "   ", ii)
        vidx_mirror[ji] = ii
    else:
        vidx_mirror[ji] = -999
        
#print(xx)

sys.exit(0)

for jt in range(Nt):
    print(" *** doing record #", jt)
    xin[:,:] = f_in.variables[cv_in][jt,:,:]
    tin      = f_in.variables[cv_tim][jt]
    
    id_tim[jt]     = tin
    id_out[jt,:,:] = xin[:,:]


#id_tim.units = 'months'
#id_out.units = 'kg/m^2/s'


#f_out.About  = 'Land-sea mask built out of variable `'+cv_nc+'` of file `'+path.basename(cf_in)+'` !'
#f_out.Author = 'Generated with `'+path.basename(sys.argv[0])+'` of `climporn` (https://github.com/brodeau/climporn)'

f_out.close()



print(cf_out+' created!!!')

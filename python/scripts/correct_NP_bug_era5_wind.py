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

vlon360    = nmp.zeros(Ni)
vlon360[:] = vlon[:]
idx = nmp.where( vlon < 0. )
vlon360[idx] = vlon[idx] + 360.

#print(vlon[:])
#print(' ' )
#print(vlon360[:])
#sys.exit(0)

xin = nmp.zeros((Nj,Ni))
###############################################################

f_out = Dataset(cf_out, 'w', format='NETCDF4')

# Dimensions:
f_out.createDimension(cv_tim, None)
f_out.createDimension(cv_lat, Nj)
f_out.createDimension(cv_lon, Ni)

id_tim  = f_out.createVariable(cv_tim, 'f8',(cv_tim,))
id_lat  = f_out.createVariable(cv_lat, 'f8',(cv_lat,))
id_lon  = f_out.createVariable(cv_lon, 'f8',(cv_lon,))

id_lat[:] = vlat[:]
id_lon[:] = vlon[:]

id_out    = f_out.createVariable(cv_in, 'f4',(cv_tim,cv_lat,cv_lon,), zlib=True, complevel=9)


# Need to find oposite longitude:
#  exp: 90 => 270
#vidx_mirror = nmp.zeros(Ni,dtype=int)
#for ji in range(Ni):
#    rl = vlon360[ji]
#    if ji%20 == 0: print('\n rl = ', rl)
#    rt = rl + 180.
#    if rt > 360.: rt = rt - 360.
#    if rl != 180.:
#        [[ii]] = nmp.argwhere(vlon360 == rt )
#        if ji%20 == 0: print(' rt = ', rt, ' --- lon = ',rl, ' ==> mirror =', vlon360[ii], '   ', ii)#
#        vidx_mirror[ji] = ii
#    else:
#        vidx_mirror[ji] = -999
        

for jt in range(Nt):
    print(" *** doing record #", jt)
    xin[:,:] = f_in.variables[cv_in][jt,:,:]
    tin      = f_in.variables[cv_tim][jt]

    # Fixing North-Pole (jj==0)
    for ji in range(Ni):
        xin[0,ji] = (vlon[0]-vlon[1]) * ( 2.*(xin[1,ji]-xin[2,ji])/(vlon[1]-vlon[2]) - (xin[2,ji]-xin[3,ji])/(vlon[2]-vlon[3])) + xin[1,ji]
    
    id_tim[jt]     = tin
    id_out[jt,:,:] = xin[:,:]


id_tim.units = f_in.variables[cv_tim].units
id_lat.units = f_in.variables[cv_lat].units
id_lon.units = f_in.variables[cv_lon].units
id_out.units = f_in.variables[cv_in].units

#id_tim.long_name = f_in.variables[cv_tim].long_name
id_lat.long_name = f_in.variables[cv_lat].long_name
id_lon.long_name = f_in.variables[cv_lon].long_name
id_out.long_name = f_in.variables[cv_in].long_name

f_out.About  = 'Original ERA5 file with vector component fixed at the North-Pole raw (j=0) (Akima extrapolation)'
f_out.Author = 'Generated with `'+path.basename(sys.argv[0])+'` of `climporn` (https://github.com/brodeau/climporn)'

f_out.close()



print(cf_out+' created!!!')

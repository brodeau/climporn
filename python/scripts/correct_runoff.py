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

l_use_fillval = True

cv_coefr  = 'socoefr'
cv_runoff = 'sorunoff'

narg = len(sys.argv)
if narg not in [4]:
    print('Usage: '+sys.argv[0]+' <runoff.nc> <mesh_mask.nc> <highest_value_allowed>')
    sys.exit(0)

cf_ro_in = sys.argv[1]
cf_mm_in = sys.argv[2]
rmaxv = float(sys.argv[3])
    
cf_ro_out = 'NEW_'+cf_ro_in
 
print(' *** Will create new file: '+cf_ro_out)


# Read mask to enforce:
f_mm_in = Dataset(cf_mm_in)
xmsk   = f_mm_in.variables['tmask'][0,0,:,:]
xlon   = f_mm_in.variables['nav_lon'][:,:]
xlat   = f_mm_in.variables['nav_lat'][:,:]
f_mm_in.close()

f_ro_in = Dataset(cf_ro_in)
#xcoef   = f_ro_in.variables[cv_coefr] [:,:]
xrnof   = f_ro_in.variables[cv_runoff][:,:,:]
f_ro_in.close()




(Nt, Nj,Ni) = nmp.shape(xrnof)
print("Shape of domain: Ni, Nj =", Ni, Nj)

if nmp.shape(xmsk) != (Nj,Ni):
    print("ERROR: runoff and mask do not match in shape !!!")
    sys.exit(0)


# Cleaning
ids = nmp.where( xrnof > rmaxv )
xrnof[ids] = rmaxv

ids = nmp.where( xrnof < 0. )
xrnof[ids] = 0.

xcoef = nmp.zeros((Nj,Ni))

for jt in range(Nt):
    xrnof[jt,:,:] = xrnof[jt,:,:]*xmsk[:,:]
    iok = nmp.where( xrnof[jt,:,:] > 0. )
    xcoef[iok] = 0.5
    
xcoef[:,:]        =    xcoef[:,:]*xmsk[:,:]

#cpu.smoother(xcoef, xmsk*0+1, nb_smooth=10)

f_out = Dataset(cf_ro_out, 'w', format='NETCDF4')

# Dimensions:
cdim_t = 'time_counter'
cdim_y = 'y'
cdim_x = 'x'

f_out.createDimension(cdim_t, None)
f_out.createDimension(cdim_y, Nj)
f_out.createDimension(cdim_x, Ni)

id_tim  = f_out.createVariable('time_counter','f4',(cdim_t,))
id_lat  = f_out.createVariable('nav_lat',     'f4',(cdim_y,cdim_x,))
id_lon  = f_out.createVariable('nav_lon',     'f4',(cdim_y,cdim_x,))

id_lat[:,:] = xlat[:,:]
id_lon[:,:] = xlon[:,:]

id_co_out      = f_out.createVariable(cv_coefr,'f4',(cdim_y,cdim_x,), zlib=True, complevel=8)
id_co_out[:,:] = xcoef[:,:]

id_ro_out      = f_out.createVariable(cv_runoff,'f4',(cdim_t,cdim_y,cdim_x,), zlib=True, complevel=8)

for jt in range(Nt):
    id_tim[jt]        = float(jt+1)+0.5
    id_ro_out[jt,:,:] = xrnof[jt,:,:]


id_tim.units = 'months'
id_ro_out.units = 'kg/m^2/s'


#f_out.About  = 'Land-sea mask built out of variable `'+cv_nc+'` of file `'+path.basename(cf_ro_in)+'` !'
#f_out.Author = 'Generated with `'+path.basename(sys.argv[0])+'` of `climporn` (https://github.com/brodeau/climporn)'

f_out.close()



print(cf_ro_out+' created!!!')

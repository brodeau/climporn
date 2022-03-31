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

narg = len(sys.argv)
if narg not in [4]:
    print('Usage: '+sys.argv[0]+' <bathy_meter.nc> <force_mask.nc> <shallowest_value_allowed>')
    sys.exit(0)

cf_b_in = sys.argv[1]
cf_m_in = sys.argv[2]
rshallv = float(sys.argv[3])
    
cf_b_out = 'NEW_'+cf_b_in
 
print(' *** Will create mask '+cf_b_out)


# Read mask to enforce:
f_m_in = Dataset(cf_m_in)
xmsk   = f_m_in.variables['bw'][:,:]
f_m_in.close()


# Read mask to enforce:
f_b_in = Dataset(cf_b_in)
xlon   = f_b_in.variables['nav_lon'][:,:]
xlat   = f_b_in.variables['nav_lat'][:,:]
xbath  = f_b_in.variables['Bathymetry'][:,:]
f_b_in.close()




(Nj,Ni) = nmp.shape(xmsk)
print("Shape of domain: Ni, Nj =", Ni, Nj)

if nmp.shape(xbath) != (Nj,Ni):
    print("ERROR: bathy and mask do not match in shape !!!"); sys.exit(0)


ids = nmp.where( xbath < rshallv )

xbath[ids] = rshallv


xbath[:,:] = xbath[:,:]*xmsk[:,:]

cpu.smoother(xbath, xmsk*0+1, nb_smooth=10)


xbath[:,:] = xbath[:,:]*xmsk[:,:]

f_out = Dataset(cf_b_out, 'w', format='NETCDF4')

# Dimensions:
cdim_y = 'y'
cdim_x = 'x'

f_out.createDimension(cdim_y, Nj)
f_out.createDimension(cdim_x, Ni)

id_lon  = f_out.createVariable('nav_lon','f4',(cdim_y,cdim_x,))
id_lat  = f_out.createVariable('nav_lat','f4',(cdim_y,cdim_x,))
id_lon[:,:] = xlon[:,:]
id_lat[:,:] = xlat[:,:]


id_b_out  = f_out.createVariable('Bathymetry','f4',(cdim_y,cdim_x,), zlib=True, complevel=8)
id_b_out[:,:]   = xbath[:,:]

id_b_out.long_name = 'Bathymetry'


#f_out.About  = 'Land-sea mask built out of variable `'+cv_nc+'` of file `'+path.basename(cf_b_in)+'` !'
#f_out.Author = 'Generated with `'+path.basename(sys.argv[0])+'` of `climporn` (https://github.com/brodeau/climporn)'

f_out.close()



print(cf_b_out+' created!!!')

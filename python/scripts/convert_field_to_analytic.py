#!/usr/bin/env python3


# L. Brodeau, 2017

import sys
import os
import numpy as nmp
from netCDF4 import Dataset


if len(sys.argv) != 4:
    print('Usage: '+sys.argv[0]+' <file_in.nc> <var_in> <value>')
    sys.exit(0)
cf_in  = sys.argv[1]
cv_in  = sys.argv[2]
rval   = float(sys.argv[3])

cf_new = str.replace(cf_in, '.nc', '_NEW.nc')

print(rval, cf_new)

os.system('rm -f '+cf_new)
os.system('cp '+cf_in+' '+cf_new)

print('\n')





#(Nj,Ni) = nmp.shape(xbathy)

#xnew = nmp.zeros((Nj,Ni))

print('\n')


# Opening the Netcdf file:
f_new = Dataset(cf_new, 'r+')     # r+ => can read and write in the file... )
print('File ', cf_new, 'is open...\n')

list_dim_var = list( f_new.variables[cv_in].dimensions )
nb_dim = len(list_dim_var)
l_time_rec = ( ('time_counter' in list_dim_var) or ('time' in list_dim_var) )

#print(list_dim_var, l_time_rec, list_dim_var[0])
if nb_dim==3 and l_time_rec:
    Nt = f_new.dimensions[list_dim_var[0]].size
    Nj = f_new.dimensions[list_dim_var[1]].size
    Ni = f_new.dimensions[list_dim_var[2]].size    
else:
    print('FIXME: unknown combination of dimmensions!')
    sys.exit(0)

print(' *** Ni, Nj, Nt = ', Ni, Nj, Nt)


# Updating field `cv_in`:
if nb_dim==3 and l_time_rec:
    for jt in range(Nt):
        f_new.variables[cv_in][jt,:,:] = rval


f_new.close()

print(cf_new+' sucessfully created!')


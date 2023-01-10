#!/usr/bin/env python3
#
# Takes a file `file_in.nc`, makes a copy of it: `file_in_NEW.nc`, and in this new file
# Fill the rectangular region defined by <i1,j1,i2,j2> with extrapolated values from the surroundings
#
# L. Brodeau, 2023
#

import sys
import os
from re import split
import numpy as np
from netCDF4 import Dataset

import climporn as cp

ldebug = True

if not len(sys.argv) in [4,5]:
    print('Usage: '+sys.argv[0]+' <file_in.nc> <var_in> <i1,j1,i2,j2> (<ilev>)')
    sys.exit(0)
cf_in = sys.argv[1]
cv_in = sys.argv[2]
cbox  = sys.argv[3]

l_3d = False
if len(sys.argv)==5:
    l_3d = True
    jlev = int(sys.argv[4])

vv = np.array( split(',', cbox), dtype=int)

print('\n')

[i1,j1, i2,j2] = vv
print(' *** i1,j1, i2,j2 =', i1,j1, i2,j2)
if l_3d: print('     level =',jlev)

#######

cf_new = str.replace(cf_in, '.nc', '_NEW.nc')

#print(rval, cf_new)

os.system('rm -f '+cf_new)
os.system('cp '+cf_in+' '+cf_new)

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
    Nk = -1
    print(' *** Ni, Nj, Nt = ', Ni, Nj, Nt)
    #
if nb_dim==4 and l_time_rec:
    Nt = f_new.dimensions[list_dim_var[0]].size
    Nk = f_new.dimensions[list_dim_var[1]].size
    Nj = f_new.dimensions[list_dim_var[2]].size
    Ni = f_new.dimensions[list_dim_var[3]].size
    print(' *** Ni, Nj, Nk, Nt = ', Ni, Nj, Nk, Nt)
    #
else:
    print('FIXME: unknown combination of dimmensions!')
    sys.exit(0)





rfill_val = f_new.variables[cv_in]._FillValue

print(' rfill_val =', rfill_val)

for jt in range(Nt):

    XX = f_new.variables[cv_in][jt,jlev,j1-50:j2+50,i1-50:i2+50]
    #XX = f_new.variables[cv_in][jt,jlev,:,:]

    XW = XX.copy()


    if jt==0:
        mask_orig      = XX.copy()
        mask_orig[:,:] = 1.
        im = np.where(XX.data==rfill_val)
        mask_orig[im] = 0.
        
    jA = 50 ; jB = jA + j2 - j1
    iA = 50 ; iB = iA + i2 - i1


    mask = mask_orig.copy()        
    mask[jA:jB,iA:iB] = 0.

    
    

    ii = cp.drown(XW, mask, k_ew=-1, nb_max_inc=100, nb_smooth=10)

    im = np.where(mask_orig==0.)
    XW[im] = rfill_val

    if ldebug:
        cp.dump_2d_field('0field_zoom'+str(jt)+'.nc', XX)
        cp.dump_2d_field('0mask_zoom'+str(jt)+'.nc', mask)
        cp.dump_2d_field('0drown_zoom'+str(jt)+'.nc', XW)


    f_new.variables[cv_in][jt,jlev,j1-50:j2+50,i1-50:i2+50] = XW[:,:]
    

f_new.close()

print(cf_new+' sucessfully created!')


#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#     CLIMPORN
#
#  Plot a NEMO field directly on the NEMO grid (no projection)
#  => 1 image pixel per grid point 
#
#    L. Brodeau, May 2018
##################################################################

import sys
from os import path
import numpy as nmp

from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import climporn as cp

bathy_max = 5000. # m

color_top = 'w'
clr_yellow = '#ffed00'

# Defaults:
lknown = True
rfact_zoom = 1.
l_show_cb = False
l_show_nm = False
l_scientific_mode = False
l_show_ttl = False
vcb = [0.15, 0.96, 0.8, 0.02]
font_rat = 1.



l_show_msh = False
    
pt_sz_track = 30

l_read_lsm=False

l_use_xmsk=True


fig_type='png'

narg = len(sys.argv)
if not narg in [6,7,8]:
    print('Usage: '+sys.argv[0]+' <CONF> <file> <variable> <snapshot> <name_fig> (<LSM_file> <point>)')
    sys.exit(0)
CNEMO  = sys.argv[1]
cf_fld = sys.argv[2]
cv_in  = sys.argv[3]
jt = int(sys.argv[4])
cnfig  = sys.argv[5]

if narg >= 7 :
    l_read_lsm=True
    cf_lsm = sys.argv[6]

cpnt = 't'
if narg == 8 :
    cpnt = sys.argv[7]

if not cpnt in ['t','f','u','v']:
    print('ERROR: what to do with C-grid "'+cpnt+'" point!?')
    sys.exit(0)

l_bathy_var = [ 'Bathymetry', 'elevation' ]
    

#if not l_read_lsm and ( not cv_in in l_bathy_var):
#    print("It's only for bathymetric fields that you can skip providing the mesh_mask file!")
#    sys.exit(0)



dir_conf = path.dirname(cf_fld)
if dir_conf == '':  dir_conf = '.'
print('\n *** dir_conf =',dir_conf,'\n')


l_add_topo_land = False
rof_log = 150.
    
i2=0
j2=0

if CNEMO == 'NATL60':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.25 ; vcb = [0.6, 0.1, 0.39, 0.025] ; font_rat = 5.*rfact_zoom
    x_cnf = 160. ; y_cnf = 2300. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False
    
elif CNEMO == 'NANUK1':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 4. ; vcb = [0.5, 0.875, 0.49, 0.02] ; font_rat = 0.16*rfact_zoom
elif CNEMO == 'NANUK1h': # half proc [i/2]
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 8. ; vcb = [0.04, 0.06, 0.92, 0.015] ; font_rat = 0.1*rfact_zoom
    l_show_cb = True ; color_top = 'k'

elif CNEMO == 'eHUDSON4h': # half proc [i/2]
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 12. ; vcb = [0.04, 0.06, 0.92, 0.015] ; font_rat = 0.1*rfact_zoom
    l_show_cb = False ; color_top = 'k'

elif CNEMO == 'HUDSON12': # half proc [i/2]
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.04, 0.06, 0.92, 0.015] ; font_rat = 0.1*rfact_zoom
    l_show_cb = False ; color_top = 'k'

elif CNEMO == 'NANUK4':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.63, 0.95, 0.33, 0.015] ; font_rat = 0.1*rfact_zoom
    l_show_cb=True
    l_show_nm=True ; x_cnf=15*rfact_zoom ; y_cnf=250*rfact_zoom ; 
    
elif CNEMO == 'NANUK025':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.5, 0.875, 0.49, 0.02] ; font_rat = 0.5*rfact_zoom
    x_cnf = 30. ; y_cnf = 540. ; # where to put label of conf on Figure...

elif CNEMO == 'ROALD12':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.63, 0.95, 0.36, 0.02] ; font_rat = 1.*rfact_zoom
    x_cnf = 50. ; y_cnf = 1250. ; # where to put label of conf on Figure...
    
elif CNEMO in [ 'CREG025' ] :
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2.
    vcb = [0.6, 0.975, 0.38, 0.02] ; font_rat = 0.5*rfact_zoom
    x_cnf = 20. ; y_cnf = 560. ; # where to put label of conf on Figure...

elif CNEMO in [ 'CREG4' ] :
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2.
    vcb = [0.14, 0.05, 0.8, 0.02] ; font_rat = 0.5*rfact_zoom
    x_ttl = 210. ; y_ttl = 620. ; # where to put label of conf on Figure...
    l_show_nm = False ; l_show_msh = True
    l_scientific_mode = True ; l_show_ttl = True
    color_top = 'k'

elif CNEMO == 'eNATL4':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.25 ; vcb = [0.6, 0.11, 0.39, 0.025] ; font_rat = 0.5*rfact_zoom
    x_cnf = 20. ; y_cnf = 560. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False

elif CNEMO == 'eNATL36':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.3 ; vcb = [0.6, 0.1, 0.39, 0.025] ; font_rat = 2.5*rfact_zoom
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False

elif CNEMO == 'eNATL60':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.25 ; vcb = [0.6, 0.1, 0.39, 0.025] ; font_rat = 5.*rfact_zoom
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False

elif CNEMO == 'eNATL1':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 6.
    vcb = [0.62, 0.11, 0.35, 0.025] ; font_rat = 0.12*rfact_zoom
    x_cnf = 4. ; y_cnf = 120. ; # where to put label of conf on Figure...

elif CNEMO == 'SouthPac':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.5, 0.875, 0.49, 0.02] ; font_rat = 0.16*rfact_zoom
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False
    bathy_max = 8000. # m
    
elif CNEMO == 'SWEPAC2':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 12. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 1./rfact_zoom*10.
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False
    bathy_max = 6000. # m
    
elif CNEMO == 'TROPICO2':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 10. ; vcb = [0.35, 0.94, 0.6, 0.04] ; font_rat = 8./rfact_zoom
    x_cnf = 20. ; y_cnf = 3. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = True ; l_scientific_mode=False
    bathy_max = 6000. # m
    
elif CNEMO == 'TROPICO05':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 5. ; vcb = [0.35, 0.09, 0.4, 0.03] ; font_rat = 4./rfact_zoom
    x_cnf = 280. ; y_cnf = 135. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False ; l_scientific_mode=False
    bathy_max = 6000. # m
    
elif CNEMO == 'TROPICO12':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.35, 0.08, 0.4, 0.03] ; font_rat = 1.1/rfact_zoom
    x_cnf = 1400. ; y_cnf = 820. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False ; l_scientific_mode=False
    bathy_max = 6000. # m
    l_add_topo_land = True ; ftopo = 'z_ETOPO1_21601x10801-TROPICO12_ice.nc' ; rof_log = 300.

elif CNEMO == 'GEBCO':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.2, 0.06, 0.6, 0.03] ; font_rat = 1.5/rfact_zoom
    x_cnf = 900. ; y_cnf = 1350. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = True ; l_scientific_mode=False
    bathy_max = 5000. # m
    
elif CNEMO == 'Azores':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 1./rfact_zoom*10.
    x_cnf = 0. ; y_cnf = 0. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False
    bathy_max = 6000. # m

elif CNEMO == 'GulfS':
    i1 = 0 ; j1 = 420 ; i2 = 900 ; j2 = 980 ; rfact_zoom = 1. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3

elif CNEMO == 'Faroe':
    i1 = 0 ; j1 = 0 ; i2 = 421 ; j2 = 351 ; rfact_zoom = 1. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3

elif CNEMO == 'WestMed':
    i1 = 0 ; j1 = 0 ; i2 = 868 ; j2 = 796 ; rfact_zoom = 1. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3

elif CNEMO == 'SouthWestPac_G12':
    i1 = 0 ; j1 = 0 ; i2 = 601 ; j2 = 301 ; rfact_zoom = 1. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3

elif CNEMO == 'ORCA1':
    i1 = 0 ; j1 = 0 ; i2 = 362 ; j2 = 292 ; rfact_zoom = 2. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3

elif CNEMO == 'CALEDO10':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.2, 0.06, 0.6, 0.03] ; font_rat = 1.5/rfact_zoom
    x_cnf = 900. ; y_cnf = 1350. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False ; l_scientific_mode=False
    bathy_max = 6000. # m
    
elif CNEMO == 'CALEDO60':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.2, 0.06, 0.6, 0.03] ; font_rat = 1.5/rfact_zoom
    x_cnf = 900. ; y_cnf = 1350. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False ; l_scientific_mode=False
    bathy_max = 6000. # m
    
else:
    print('\n WARNING [nemo_imshow_2d_field.py]: "'+CNEMO+'" is an unknown config!\n     ==> falling back on default setup')
    lknown = False
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0



laplacian = False
l_log_field = False
l_pow_field = False
cextend='both'
l_hide_cb_ticks = False
tmin=0. ; tmax=1. ; df=0.01
cb_jump = 1


print(' cv_in = '+cv_in)

if cv_in in ['sosstsst','tos']:
    cfield = 'SST'
    tmin=0. ;  tmax=28.   ;  df = 1.
    cpal_fld = 'ncview_nrl'    
    cunit = r'SST ($^{\circ}$C)'


elif cv_in in l_bathy_var:   # 
    cfield = 'Bathymetry'       # 
    tmin=0. ;  tmax=bathy_max   ;  df = 100.
    #cpal_fld = 'ocean'
    #cpal_fld = 'Blues'
    #cpal_fld = 'PuBu'
    cpal_fld = 'on2_r'
    #cpal_fld = 'ncview_ssec'
    #cpal_fld = 'ncview_hotres'
    #cpal_fld = 'ncview_helix'
    cunit = r'Bathymetry (m)'
    l_pow_field = True
    pow_field = 1.3
    #l_log_field = False
    cextend='max'
    #l_hide_cb_ticks=True
    if cv_in == 'elevation':
        #cpal_fld = 'ncview_hotres' ; l_pow_field = False
        cpal_fld = 'ncview_nrl' ; l_pow_field = True
        #cpal_fld = 'ocean' ; l_pow_field = False
    
elif cv_in == 'sossheig':
    cfield = 'SSH'
    tmin=-0.3 ;  tmax=0.3   ;  df = 0.1
    cpal_fld = 'ncview_jaisnc'    
    cunit = r'SSH (m)'


elif cv_in == 'somxl010':
    cfield == 'MLD'
    tmin=50. ;  tmax=1500. ;  df = 50.
    cpal_fld = 'viridis_r'

elif cv_in == 'track':
    cfield = 'TRACK'
    cpal_fld = 'nipy_spectral'
    cunit = r'SST ($^{\circ}$C)'
    fig_type='svg'

elif cv_in in [ 'damage', 'damage-t', 'damage-f', 'dmg', 'dmgf', 'dmgt' ]:
    cfield = 'damage'
    tmin=0. ;  tmax=1.   ;  df = 0.1
    cpal_fld = 'inferno'
    #cpal_fld = 'cividis'
    cunit = r'damage'
    l_pow_field = True
    pow_field = 3.5

elif cv_in in [ 'siconc' ]:
    cfield = 'ice_frac'
    tmin=0. ;  tmax=1.   ;  df = 0.1
    #cpal_fld = 'viridis'
    cpal_fld = 'bone'
    cunit = r'sea-ice fraction'
    l_pow_field = True
    pow_field = 2.5
    
elif cv_in in [ 'zfU', 'zfV', 'ds11dx', 'ds22dy', 'ds12dx', 'ds12dy', 'zfUv', 'zfVu' ]:
    cfield = 'divS'
    tmin=-1. ;  tmax=1.   ;  df = 0.1
    cpal_fld = 'RdBu_r'
    cunit = r'[Pa]'

elif cv_in in [ 'u_oceT' ]:
    cfield = 'u_oce'
    tmin=-0.25 ;  tmax=0.25   ;  df = 0.05
    cpal_fld = 'RdBu_r'
    cunit = r'[m/s]'
elif cv_in in [ 'v_oceT' ]:
    cfield = 'v_oce'
    tmin=-0.25 ;  tmax=0.25   ;  df = 0.05
    cpal_fld = 'RdBu_r'
    cunit = r'[m/s]'

elif cv_in in [ 'taux_ai' ]:
    cfield = 'taux_ice'
    tmin=-0.5 ;  tmax=0.5   ;  df = 0.025
    cpal_fld = 'RdBu_r'
    cunit = r'[Pa]'
elif cv_in in [ 'tauy_ai' ]:
    cfield = 'tauy_ice'
    tmin=-0.5 ;  tmax=0.5   ;  df = 0.025
    cpal_fld = 'RdBu_r'
    cunit = r'[Pa]'
    
elif cv_in in [ 's11', 'sig1', 'sig11', 'sig1F' ]:
    cfield = 'sig11'
    tmin=-50000. ;  tmax=50000.   ;  df = 25000.
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_{11}$ (N)'

elif cv_in in [ 's22', 'sig2', 'sig22', 'sig2F' ]:
    cfield = 'sig22'
    tmin=-30000. ;  tmax=30000.   ;  df = 15000.
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_{22}$ (N)'

elif cv_in in [ 's12', 'sig12', 'sig12f', 'sig12T' ] or cv_in[0:5]=='sig12' :
    cfield = 'sig12'
    tmin=-30000. ;  tmax=30000.   ;  df = 15000.
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_{12}$ (N)'

elif cv_in in [ 'Uice', 'Ut', 'u_ice', 'uVice' ]:
    cfield = 'Uice'
    tmin=-0.3 ;  tmax=0.3   ;  df = 0.2
    cpal_fld = 'RdBu_r'
    cunit = r'$u_{ice}$ (m/s)'

elif cv_in in [ 'Vice', 'Vt', 'v_ice', 'vUice' ]:
    cfield = 'Vice'
    tmin=-0.3 ;  tmax=0.3   ;  df = 0.2
    cpal_fld = 'RdBu_r'
    cunit = r'$v_{ice}$ (m/s)'

elif cv_in in [ 'zsN', 'MC', 'sigI', 'sigII' ]:
    cfield = 'sigma'
    tmin=-30000. ;  tmax=30000.   ;  df = 15000.
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_N$ (N)'

elif cv_in == 'zsS':
    cfield = 'sigS'
    tmin=-30000. ;  tmax=30000.   ;  df = 15000.
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_S$ (N)'
elif cv_in == 'zsS2':
    cfield = 'sigS'
    tmin=-0.5E9 ;  tmax=-tmin   ;  df = 1.E8
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_S$ (N)'

elif cv_in in ['e11','e22','e12','div','ten','shr']:
    cfield = 'eps'
    tmin=-1.e-5 ;  tmax=-tmin   ;  df = 5.E-6
    cpal_fld = 'RdBu_r'
    cunit = r'$\epsilon$ (m/s^2)'

elif cv_in in ['shear2']:
    cfield = 'shear2'
    tmin=-1.e-10 ;  tmax=-tmin   ;  df = 5.E-9
    cpal_fld = 'RdBu_r'
    cunit = r'$\tau$^2'

elif cv_in in ['elasticity','elasticity-t','elasticity-f']:
    cfield = 'E'
    tmin=1.E8 ;  tmax=6.E8   ;  df = 1.E8
    cpal_fld = 'viridis'
    cunit = r'Elasticity'

elif cv_in in ['lambda-t','lambda-f']:
    cfield = 'lambda'
    tmin=1.E6 ;  tmax=1.E7   ;  df = 2.E6
    cpal_fld = 'viridis'
    cunit = r'$\lambda$'

elif cv_in in ['Pmax-t','Pmax-f']:
    cfield = 'Pmax'
    tmin=0. ;  tmax=1.E5   ;  df = 2.E4
    cpal_fld = 'ncview_tofino'
    cunit = r'$P_{max}$'

elif cv_in in ['Ptilde-t','Ptilde-f']:
    cfield = 'Ptilde'
    tmin=0. ;  tmax=1.   ;  df = 0.1
    cpal_fld = 'ncview_tofino'
    cunit = r'$\tilde{P}$'

elif cv_in in ['mult-t','mult-f']:
    cfield = 'mult'
    tmin=0. ;  tmax=1.   ;  df = 0.1
    cpal_fld = 'ncview_tofino'
    cunit = r'multiplicator'

elif cv_in in ['dmg-inc-t','dmg-inc-f']:
    cfield = 'd_inc'
    tmin=-1.e-3 ;  tmax=-tmin   ;  df = 2.e-4
    cpal_fld = 'ncview_tofino'
    cunit = r'multiplicator'

elif cv_in in ['dcrit-t','dcrit-f']:
    cfield = 'dcrit'
    tmin=-5 ;  tmax=5.   ;  df = 1.
    #cpal_fld = 'ncview_tofino'
    cpal_fld = 'RdBu_r'
    cunit = r'$d_{crit}$'

elif cv_in in ['Tdc-t','Tdc-f']:
    cfield = 'Tdc'
    tmin=0. ;  tmax=300.   ;  df = 50.
    cpal_fld = 'ncview_tofino'
    cunit = r'$Td_{c}$'

elif cv_in in ['mult-dmg-t','mult-dmg-f']:
    cfield = 'mult-dmg'
    tmin=0. ;  tmax=0.2   ;  df = 0.05
    cpal_fld = 'ncview_tofino'
    cunit = r'$MultDmg$'


    
else:
    print('ERROR: variable '+cv_in+' is not known yet...'); sys.exit(0)




# Time record stuff...
cp.chck4f(cf_fld)
id_fld = Dataset(cf_fld)
list_var = id_fld.variables.keys()
if 'time_counter' in list_var:    
    vtime = id_fld.variables['time_counter'][:]
    Nt = len(vtime)
    print('\n There is a "time_counter" in file '+cf_fld+' !')
    print('   => '+str(Nt)+' snapshots!')
    if jt <= 0 or jt > Nt: print(' PROBLEM: your time record does not exist!', jt) ; sys.exit(0)
else:
    print('\nWARNING: there is NO "time_counter" in file '+cf_fld+' !')
    print('   ==> setting Nt = 0 !\n')
    Nt = 0
id_fld.close()

ibath=1

if l_read_lsm:
    cp.chck4f(cf_lsm)
    cnmsk = cpnt+'mask'
    print('\n *** Reading "'+cnmsk+'" in meshmask file...')
    id_lsm = Dataset(cf_lsm)
    nb_dim = len(id_lsm.variables[cnmsk].dimensions)
    Ni = id_lsm.dimensions['x'].size
    Nj = id_lsm.dimensions['y'].size
    if i2 == 0: i2 = Ni
    if j2 == 0: j2 = Nj
    if nb_dim == 4: XMSK  = id_lsm.variables[cnmsk][0,0,j1:j2,i1:i2] ; # t, y, x
    if nb_dim == 3: XMSK  = id_lsm.variables[cnmsk][0,  j1:j2,i1:i2] ; # t, y, x
    if nb_dim == 2: XMSK  = id_lsm.variables[cnmsk][    j1:j2,i1:i2] ; # t, y, x
    if l_show_msh:
        Xlon = id_lsm.variables['glamu'][0,j1:j2,i1:i2]
        Xlat = id_lsm.variables['gphiv'][0,j1:j2,i1:i2]
    id_lsm.close()
    print('      done.')

elif cv_in in l_bathy_var:
    cp.chck4f(cf_fld)
    print('\n *** Will build mask from "Bathymetry"...')
    id_fld = Dataset(cf_fld)
    list_dim = list(id_fld.dimensions.keys()) ; #lolopy3
    nb_dim = len(id_fld.variables[cv_in].dimensions)
    for jd in range(nb_dim):
        if list_dim[jd] in ['x','lon','longitude','glamt']: cdim_x = list_dim[jd]
        if list_dim[jd] in ['y','lat','latitude' ,'gphit']: cdim_y = list_dim[jd]
    print(' *** x, y dims =>', cdim_x, cdim_y)
    #
    Ni = id_fld.dimensions[cdim_x].size
    Nj = id_fld.dimensions[cdim_y].size
    if i2 == 0: i2 = Ni
    if j2 == 0: j2 = Nj
    if nb_dim == 3: XBATH = id_fld.variables[cv_in][0,  j1:j2,i1:i2] ; # t, y, x
    if nb_dim == 2: XBATH = id_fld.variables[cv_in][    j1:j2,i1:i2] ; # t, y, x
    if l_show_msh:
        Xlon = id_fld.variables['nav_lon'][j1:j2,i1:i2]
        Xlat = id_fld.variables['nav_lat'][j1:j2,i1:i2]
    id_fld.close()

    # Does bathymetry come as positive or negative in file???
    bt_max = nmp.max(XBATH[:,:])
    bt_min = nmp.min(XBATH[:,:])
    print('min, max bathy =>', bt_max, bt_min)

    if abs(bt_max) > 12000. or abs(bt_min) > 12000.:
        print('PROBLEM #1: we dont know what do do with min and max of bathymetry!'); sys.exit(0)
    if abs(bt_min) > abs(bt_max):
        print(' *** Bathymetry seems negative!')
        ibath = -1
        XBATH[:,:] = XBATH[:,:]*ibath
    
    XMSK = nmp.zeros((Nj,Ni))
    idx_oce = nmp.where(XBATH > 0.2)
    XMSK[idx_oce] = 1.
    del XBATH
    print('      done.')

elif cv_in == 'track':
    cp.chck4f(cf_fld)
    id_fld = Dataset(cf_fld)
    xtmp = id_fld.variables[cv_in][:,:]
    (Nj,Ni) = xtmp.shape
    if i2 == 0: i2 = Ni
    if j2 == 0: j2 = Nj
    XMSK = nmp.zeros((Nj,Ni), dtype=nmp.int) ; XMSK[:,:] = 1
    XMSK[nmp.where(xtmp==-100.)] = 0
    del xtmp
    
else:
    #print('PROBLEM #2'); sys.exit(0)
    print('WARNING: no mask provided!!!')
    l_use_xmsk = False
    cp.chck4f(cf_fld)
    id_fld = Dataset(cf_fld)
    xtmp = id_fld.variables[cv_in][:,:]
    (Nj,Ni) = xtmp.shape
    if i2 == 0: i2 = Ni
    if j2 == 0: j2 = Nj


print('\n The shape of the domain is Ni, Nj =', Ni, Nj)
#sys.exit(0)

# Show topo ?
if l_add_topo_land and l_use_xmsk:
    cf_topo_land = dir_conf+'/'+ftopo
    print('\n We are going to show topography:\n'+'  ==> '+cf_topo_land)

    cp.chck4f(cf_topo_land)
    id_top = Dataset(cf_topo_land)
    print(' *** Reading "z" into:\n'+cf_topo_land)
    xtopo = id_top.variables['z'][0,j1:j2,i1:i2]
    id_top.close()
    if nmp.shape(xtopo) != (Nj,Ni):
        print('ERROR: topo and mask do not agree in shape!'); sys.exit(0)
    xtopo = xtopo*(1. - XMSK)
    xtopo[nmp.where( XMSK  > 0.01)] = nmp.nan
    #xtopo[nmp.where( xtopo < 0.01)] = nmp.nan ; # checked the *_filled stuff instead...
    print('')






# Stuff for size of figure respecting pixels...
print('  *** we are going to show: i1,i2,j1,j2 =>', i1,i2,j1,j2, '\n')
nx_res = i2-i1
ny_res = j2-j1
yx_ratio = float(ny_res)/float(nx_res)
if not lknown:
    rfact_zoom = round(1000./float(ny_res),1)
nxr = int(rfact_zoom*nx_res) ; # widt image (in pixels)
nyr = int(rfact_zoom*ny_res) ; # height image (in pixels)
rDPI = 400
rh  = float(nxr)/float(rDPI) ; # width of figure as for figure...

print('\n *** width and height of image to create:', nxr, nyr, '\n')




if l_use_xmsk:
    pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)
    idx_oce = nmp.where(XMSK[:,:] > 0.5)

#font_rat
#params = { 'font.family':'Ubuntu',
params = { 'font.family':'Open Sans',
           'font.weight':    'normal',
           'font.size':       int(12.*font_rat),
           'legend.fontsize': int(22.*font_rat),
           'xtick.labelsize': int(18.*font_rat),
           'ytick.labelsize': int(18.*font_rat),
           'axes.labelsize':  int(15.*font_rat) }
mpl.rcParams.update(params)
cfont_clb  = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*font_rat), 'color':color_top }
cfont_date = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(12.*font_rat), 'color':'w' }
cfont_mail = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*font_rat), 'color':'0.8'}
cfont_cnfn = { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(35.*font_rat), 'color':'w' }
cfont_axis  = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*font_rat), 'color':color_top }
cfont_ttl = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(25.*font_rat), 'color':color_top }


# Colormaps for fields:
pal_fld = cp.chose_colmap(cpal_fld)
if l_log_field:
    norm_fld = colors.LogNorm(  vmin = tmin, vmax = tmax, clip = False)
if l_pow_field:
    norm_fld = colors.PowerNorm(gamma=pow_field, vmin = tmin, vmax = tmax, clip = False)
else:
    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

pal_lsm = cp.chose_colmap('land_dark')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

pal_filled = cp.chose_colmap('gray_r')
norm_filled = colors.Normalize(vmin = 0., vmax = 0.1, clip = False)

cfig = cnfig+'.'+fig_type

#if Nt == 0:
#    cfig = cv_in+'_'+CNEMO+'_'+cpal_fld+'.'+fig_type    
#else:
#    cfig = 'snapshot_'+str(jt)+'_'+cv_in+'_'+CNEMO+'_'+cpal_fld+'.'+fig_type    

#rextra_height = 1.
#if l_scientific_mode: rextra_height = 1.12
#fig = plt.figure(num = 1, figsize=(rh,rh*yx_ratio*rextra_height), dpi=None, facecolor='w', edgecolor='0.5')


fsize = ( rh, rh*yx_ratio )

fig = plt.figure(num = 1, figsize=fsize, dpi=rDPI, facecolor='k', edgecolor='k')

if l_scientific_mode:
    ax  = plt.axes([0.09, 0.09, 0.9, 0.9], facecolor = 'r')
else:
    ax  = plt.axes([0., 0., 1., 1.],     facecolor = '0.4')

vc_fld = nmp.arange(tmin, tmax + df, df)


print('\n *** Opening file '+cf_fld)
id_fld = Dataset(cf_fld)
if Nt > 0:
    print('    => Reading record #'+str(jt)+' of '+cv_in+' in '+cf_fld)
    XFLD  = id_fld.variables[cv_in][jt-1,j1:j2,i1:i2] ; # t, y, x
else:
    print('    => Reading 2D field '+cv_in+' in '+cf_fld+' (no time records...)')
    XFLD  = id_fld.variables[cv_in][j1:j2,i1:i2] ; # t, y, x

id_fld.close()
print('          Done!\n')

if l_use_xmsk:
    if XMSK.shape != XFLD.shape:
        print('\n PROBLEM: field and mask do not agree in shape!')
        print(XMSK.shape , XFLD.shape)
        sys.exit(0)
    print('  *** Shape of field and mask => ', nmp.shape(XFLD))



if cv_in in l_bathy_var and ibath==-1: XFLD = ibath*XFLD

l_add_true_filled = False

if cfield == 'Bathymetry':
    (idy_nan,idx_nan) = nmp.where( nmp.isnan(XFLD) )
    #
    # LSM with different masking for true lsm and filled lsm...
    cf_mask_lbc = dir_conf+'/lsm_LBC_'+CNEMO+'.nc'
    if path.exists(cf_mask_lbc):
        print('\n *** '+cf_mask_lbc+' found !!!')
        l_add_true_filled = True 
        id_filled = Dataset(cf_mask_lbc)
        xtmp = id_filled.variables['lsm'][j1:j2,i1:i2]
        id_filled.close()
        pfilled = nmp.ma.masked_where(xtmp[:,:] != -1., xtmp[:,:]*0.+40.)
        if l_add_topo_land and l_use_xmsk:
            xtopo[nmp.where(xtmp[:,:] < -0.9)] = nmp.nan
            XMSK[:,:] = 1
            XMSK[nmp.where(xtmp[:,:]==0)] = 0 ; # updated mask
        del xtmp

        print('  => done filling "pfilled" !\n')

    
if cv_in == 'track':

    XFLD[nmp.where(nmp.isnan(XFLD))] = -1000
    indx = nmp.where( XFLD > 0 )
    (idy,idx) = indx
    
    tmin=nmp.amin(XFLD[indx]) ;  tmax=nmp.amax(XFLD[indx])
    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

    cf = plt.scatter(idx, idy, c=XFLD[indx], cmap=pal_fld, norm=norm_fld, alpha=0.5, marker='.', s=pt_sz_track )


else:
    #cf = plt.imshow(XFLD[:,:], cmap=pal_fld, norm=norm_fld, interpolation='none')
    cf = plt.pcolormesh(XFLD[:,:], cmap=pal_fld, norm=norm_fld )
    if cfield == 'Bathymetry':
        #lulu aspect='auto'111 , interpolation='nearest'
        if len(idy_nan) > 0:
            idd = nmp.where(idy_nan==1); idy_nan[idd] = int(10./rfact_zoom)/2  # just so the boundary line is not too thin on plot...
            plt.scatter(idx_nan, idy_nan, color=clr_yellow, marker='s', s=int(10./rfact_zoom))

            
if l_show_msh:
    ccx = plt.contour(Xlon[:,:], 60, colors='k', linewidths=0.5)
    ccy = plt.contour(Xlat[:,:], 30, colors='k', linewidths=0.5)





del XFLD


if l_use_xmsk:
    if l_add_topo_land:
        #print('Ploting topography over continents...')
        #cp.dump_2d_field( 'xtopo.nc', xtopo ); #, xlon=[], xlat=[], name='field', unit='', long_name='', mask=[] )
        #cp.dump_2d_field( 'xmsk.nc', XMSK ); #, xlon=[], xlat=[], name='field', unit='', long_name='', mask=[] )
        xtopo = nmp.log10(xtopo+rof_log)
        pal_topo = cp.chose_colmap('gray_r')
        norm_topo = colors.Normalize(vmin = nmp.log10(-100. + rof_log), vmax = nmp.log10(6000. + rof_log), clip = False)
        cm = plt.imshow(xtopo, cmap=pal_topo, norm=norm_topo, interpolation='none')
        plt.contour(XMSK, [0.9], colors='k', linewidths=0.5)
        #
    else:
        #cm = plt.imshow(pmsk, cmap=pal_lsm, norm=norm_lsm, interpolation='none')
        cm = plt.pcolormesh( pmsk, cmap=pal_lsm, norm=norm_lsm )
        del XMSK



if cfield == 'Bathymetry' and l_add_true_filled:
    # Ocean that has been filled turns black:
    cfl = plt.imshow(pfilled, cmap=pal_filled, norm=norm_filled, interpolation='none' ) #, interpolation='none')


if l_use_xmsk: del pmsk
if l_add_true_filled: del pfilled


plt.axis([ 0, Ni, 0, Nj])

if l_scientific_mode:
    plt.xlabel('i-points', **cfont_axis)
    plt.ylabel('j-points', **cfont_axis)

if l_show_cb:

    ax2 = plt.axes(vcb)

    if l_pow_field or l_log_field:
        clb = mpl.colorbar.ColorbarBase(ax=ax2,               cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend=cextend)
    else:
        clb = mpl.colorbar.ColorbarBase(ax=ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend=cextend)

    if cb_jump > 1:
        cb_labs = [] ; cpt = 0
        for rr in vc_fld:
            if cpt % cb_jump == 0:
                if df >= 1.: cb_labs.append(str(int(rr)))
                if df <  1.: cb_labs.append(str(rr))
            else:
                cb_labs.append(' ')
            cpt = cpt + 1
        clb.ax.set_xticklabels(cb_labs)
    zlw = 0.3 ; #linewidth for everything in colorbar
    clb.set_label(cunit, **cfont_clb)
    clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color    
    clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor
    clb.outline.set_linewidth(zlw)
    clb.ax.set_facecolor(color_top) ;
    if l_hide_cb_ticks:
        clb.ax.tick_params(axis=u'both', which=u'both',length=0) ; # remove ticks!
    else:
        clb.ax.tick_params(color=color_top, length=1, width=zlw) ; # good-looking setup...
    plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_top) ; # set colorbar ticklabels with right color...
        
    del cf
    
if l_show_nm:  ax.annotate(CNEMO, xy=(1, 4), xytext=(x_cnf, y_cnf), **cfont_cnfn)

if l_show_ttl: ax.annotate(CNEMO, xy=(1, 4), xytext=(x_ttl, y_ttl), **cfont_ttl)



#plt.savefig(cfig, dpi=rDPI, orientation='portrait', facecolor='b', transparent=True)
plt.savefig(cfig, dpi=rDPI, orientation='portrait') #, transparent=True)
print(cfig+' created!\n')
plt.close(1)



if l_use_xmsk: del cm

del fig, ax



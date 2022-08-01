#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#     CLIMPORN
#
#  Prepare 2D maps (monthly) that will later become a movie!
#  NEMO output and observations needed
#
#    L. Brodeau, August 2022
##################################################################

import sys
from os import path
import numpy as nmp

from re import split

from netCDF4 import Dataset

import csv

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

fig_type='png'

narg = len(sys.argv)
if not narg in [6]:
    print('Usage: '+sys.argv[0]+' <CONF> <file_run_csv> <file_run,var> <name_fig> <LSM_file>')
    sys.exit(0)
#
# `file_run` => just to plot a field
    
CCONF  = sys.argv[1]
cf_csv = sys.argv[2]
vv = split(',',sys.argv[3])
cf_run = vv[0]
cv_run = vv[1]
cnfig  = sys.argv[4]
#
cf_lsm = sys.argv[5]

#cpnt = 't'
#if narg == 8 :
#    cpnt = sys.argv[7]

#if not cpnt in ['t','f','u','v']:
#    print('ERROR: what to do with C-grid "'+cpnt+'" point!?')
#    sys.exit(0)


    
print('\n Field to show in background: "'+cv_run+'" of file "'+cf_run+'" !\n')


# Testing CSV file:
ft = open( cf_csv, newline='' )
truc = csv.reader(ft, delimiter=',')

jl=0
jrec=0
for row in truc:
    jl = jl + 1
    itraj = int(row[0])      ; # ID of trajectory as an integer
    ctraj = '%4.4i'%(itraj)  ; #       "      as a 4-char long string
    if itraj == 1:
        jrec = jrec + 1
        print('\n record #',jrec,':')
        # This is a new time record as we are dealing with first trajectory (again)
    
    print( 'Line:', jl, ' => traj #',itraj, ' ROW => ', row )
    #print( row[0] , '\n')
Nrec_traj = jrec

# Trajectories to follow:
follow = [ 1, 2, 3 ] ;# fixme
NbTraj = len(follow)

ITRID = nmp.zeros( NbTraj, Nrec_traj, dtype=int ) ; #trajectory ID (integer)
COORX = nmp.zeros( NbTraj, Nrec_traj )
COORY = nmp.zeros( NbTraj, Nrec_traj )


jl=0
jrec=0
for row in truc:
    jl = jl + 1
    itraj = int(row[0])      ; # ID of trajectory as an integer
    if itraj == 1:
        jrec = jrec + 1
        # This is a new time record as we are dealing with first trajectory (again)
        jtraj=0
        
    if itraj in follow:
        ITRID[jtraj,jrec-1] = itraj
        COORX[jtraj,jrec-1] = float(row[1])
        COORY[jtraj,jrec-1] = float(row[2])

        jtraj = jtraj+1   # itteration of 1 trajectory for this particular record

sys.exit(0)

dir_conf = path.dirname(cf_csv)
if dir_conf == '':  dir_conf = '.'
print('\n *** dir_conf =',dir_conf,'\n')

i2=0
j2=0

if CCONF == 'NATL60':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.25 ; vcb = [0.6, 0.1, 0.39, 0.025] ; font_rat = 5.*rfact_zoom
    x_cnf = 160. ; y_cnf = 2300. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False
    
elif CCONF == 'NANUK1':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 4. ; vcb = [0.5, 0.875, 0.49, 0.02] ; font_rat = 0.16*rfact_zoom
elif CCONF == 'NANUK1h': # half proc [i/2]
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 8. ; vcb = [0.04, 0.06, 0.92, 0.015] ; font_rat = 0.1*rfact_zoom
    l_show_cb = True ; color_top = 'k'

elif CCONF == 'eHUDSON4h': # half proc [i/2]
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 12. ; vcb = [0.04, 0.06, 0.92, 0.015] ; font_rat = 0.1*rfact_zoom
    l_show_cb = False ; color_top = 'k'

elif CCONF == 'HUDSON12': # half proc [i/2]
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.04, 0.06, 0.92, 0.015] ; font_rat = 0.1*rfact_zoom
    l_show_cb = False ; color_top = 'k'

elif CCONF == 'NANUK4':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.2, 0.085, 0.6, 0.02] ; font_rat = 0.1*rfact_zoom
    l_show_cb=True
    
elif CCONF == 'NANUK025':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.5, 0.875, 0.49, 0.02] ; font_rat = 0.5*rfact_zoom
    x_cnf = 30. ; y_cnf = 540. ; # where to put label of conf on Figure...

elif CCONF == 'ROALD12':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.63, 0.95, 0.36, 0.02] ; font_rat = 1.*rfact_zoom
    x_cnf = 50. ; y_cnf = 1250. ; # where to put label of conf on Figure...
    
elif CCONF in [ 'CREG025' ] :
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2.
    vcb = [0.6, 0.975, 0.38, 0.02] ; font_rat = 0.5*rfact_zoom
    x_cnf = 20. ; y_cnf = 560. ; # where to put label of conf on Figure...

elif CCONF in [ 'CREG4' ] :
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2.
    vcb = [0.14, 0.05, 0.8, 0.02] ; font_rat = 0.5*rfact_zoom
    x_ttl = 210. ; y_ttl = 620. ; # where to put label of conf on Figure...
    l_show_nm = False ; l_show_msh = True
    l_scientific_mode = True ; l_show_ttl = True
    color_top = 'k'

elif CCONF == 'eNATL4':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.25 ; vcb = [0.6, 0.11, 0.39, 0.025] ; font_rat = 0.5*rfact_zoom
    x_cnf = 20. ; y_cnf = 560. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False

elif CCONF == 'eNATL36':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.3 ; vcb = [0.6, 0.1, 0.39, 0.025] ; font_rat = 2.5*rfact_zoom
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False

elif CCONF == 'eNATL60':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.25 ; vcb = [0.6, 0.1, 0.39, 0.025] ; font_rat = 5.*rfact_zoom
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False

elif CCONF == 'eNATL1':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 6.
    vcb = [0.62, 0.11, 0.35, 0.025] ; font_rat = 0.12*rfact_zoom
    x_cnf = 4. ; y_cnf = 120. ; # where to put label of conf on Figure...

elif CCONF == 'SouthPac':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.5, 0.875, 0.49, 0.02] ; font_rat = 0.16*rfact_zoom
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False
    bathy_max = 8000. # m
    
elif CCONF == 'SWEPAC2':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 12. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 1./rfact_zoom*10.
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False
    bathy_max = 6000. # m
    
elif CCONF == 'TROPICO2':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 10. ; vcb = [0.35, 0.94, 0.6, 0.04] ; font_rat = 8./rfact_zoom
    x_cnf = 20. ; y_cnf = 3. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = True ; l_scientific_mode=False
    bathy_max = 6000. # m
    
elif CCONF == 'TROPICO05':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 5. ; vcb = [0.35, 0.09, 0.4, 0.03] ; font_rat = 4./rfact_zoom
    x_cnf = 280. ; y_cnf = 135. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False ; l_scientific_mode=False
    bathy_max = 6000. # m
    
elif CCONF == 'TROPICO12':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.35, 0.08, 0.4, 0.03] ; font_rat = 1.1/rfact_zoom
    x_cnf = 1400. ; y_cnf = 820. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False ; l_scientific_mode=False
    bathy_max = 6000. # m

elif CCONF == 'GEBCO':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.2, 0.06, 0.6, 0.03] ; font_rat = 1.5/rfact_zoom
    x_cnf = 900. ; y_cnf = 1350. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = True ; l_scientific_mode=False
    bathy_max = 5000. # m
    
elif CCONF == 'Azores':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 1./rfact_zoom*10.
    x_cnf = 0. ; y_cnf = 0. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False
    bathy_max = 6000. # m

elif CCONF == 'GulfS':
    i1 = 0 ; j1 = 420 ; i2 = 900 ; j2 = 980 ; rfact_zoom = 1. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3

elif CCONF == 'Faroe':
    i1 = 0 ; j1 = 0 ; i2 = 421 ; j2 = 351 ; rfact_zoom = 1. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3

elif CCONF == 'WestMed':
    i1 = 0 ; j1 = 0 ; i2 = 868 ; j2 = 796 ; rfact_zoom = 1. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3

elif CCONF == 'SouthWestPac_G12':
    i1 = 0 ; j1 = 0 ; i2 = 601 ; j2 = 301 ; rfact_zoom = 1. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3

elif CCONF == 'ORCA1':
    i1 = 0 ; j1 = 0 ; i2 = 362 ; j2 = 292 ; rfact_zoom = 2. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3

elif CCONF == 'CALEDO10':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.2, 0.06, 0.6, 0.03] ; font_rat = 1.5/rfact_zoom
    x_cnf = 900. ; y_cnf = 1350. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False ; l_scientific_mode=False
    bathy_max = 6000. # m
    
elif CCONF == 'CALEDO60':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.2, 0.06, 0.6, 0.03] ; font_rat = 1.5/rfact_zoom
    x_cnf = 900. ; y_cnf = 1350. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False ; l_scientific_mode=False
    bathy_max = 6000. # m
    
else:
    print('\n WARNING [nemo_imshow_2d_field.py]: "'+CCONF+'" is an unknown config!\n     ==> falling back on default setup')
    lknown = False
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0



laplacian = False
l_log_field = False
l_pow_field = False
cextend='both'
l_hide_cb_ticks = False
tmin=0. ; tmax=1. ; df=0.01
cb_jump = 1


print(' cv_run = '+cv_run)

if cv_run in ['sosstsst','tos','votemper']:
    cfield = 'SST'
    tmin=-2. ;  tmax=32.   ;  df = 1.
    cpal_fld = 'ncview_nrl'    
    cunit = r'SST ($^{\circ}$C)'
    
elif cv_run == 'sossheig':
    cfield = 'SSH'
    tmin=-0.3 ;  tmax=0.3   ;  df = 0.1
    cpal_fld = 'ncview_jaisnc'    
    cunit = r'SSH (m)'


elif cv_run == 'somxl010':
    cfield == 'MLD'
    tmin=50. ;  tmax=1500. ;  df = 50.
    cpal_fld = 'viridis_r'

elif cv_run == 'track':
    cfield = 'TRACK'
    cpal_fld = 'nipy_spectral'
    cunit = r'SST ($^{\circ}$C)'
    fig_type='svg'

elif cv_run in [ 'damage', 'damage-t', 'damage-f', 'dmg', 'dmgf', 'dmgt' ]:
    cfield = 'damage'
    tmin=0. ;  tmax=1.   ;  df = 0.1
    cpal_fld = 'inferno'
    #cpal_fld = 'cividis'
    cunit = r'damage'
    l_pow_field = True
    pow_field = 3.5

elif cv_run in [ 'siconc' ]:
    cfield = 'ice_frac'
    tmin=0. ;  tmax=1.   ;  df = 0.1
    #cpal_fld = 'viridis'
    cpal_fld = 'bone'
    cunit = r'sea-ice fraction'
    l_pow_field = True
    pow_field = 2.5
    
elif cv_run in [ 'zfU', 'zfV', 'ds11dx', 'ds22dy', 'ds12dx', 'ds12dy', 'zfUv', 'zfVu' ]:
    cfield = 'divS'
    tmin=-1. ;  tmax=1.   ;  df = 0.1
    cpal_fld = 'RdBu_r'
    cunit = r'[Pa]'

elif cv_run in [ 'u_oceT' ]:
    cfield = 'u_oce'
    tmin=-0.25 ;  tmax=0.25   ;  df = 0.05
    cpal_fld = 'RdBu_r'
    cunit = r'[m/s]'
elif cv_run in [ 'v_oceT' ]:
    cfield = 'v_oce'
    tmin=-0.25 ;  tmax=0.25   ;  df = 0.05
    cpal_fld = 'RdBu_r'
    cunit = r'[m/s]'

elif cv_run in [ 'taux_ai' ]:
    cfield = 'taux_ice'
    tmin=-0.5 ;  tmax=0.5   ;  df = 0.025
    cpal_fld = 'RdBu_r'
    cunit = r'[Pa]'
elif cv_run in [ 'tauy_ai' ]:
    cfield = 'tauy_ice'
    tmin=-0.5 ;  tmax=0.5   ;  df = 0.025
    cpal_fld = 'RdBu_r'
    cunit = r'[Pa]'
    
elif cv_run in [ 's11', 'sig1', 'sig11', 'sig1F' ]:
    cfield = 'sig11'
    tmin=-50000. ;  tmax=50000.   ;  df = 25000.
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_{11}$ (N)'

elif cv_run in [ 's22', 'sig2', 'sig22', 'sig2F' ]:
    cfield = 'sig22'
    tmin=-30000. ;  tmax=30000.   ;  df = 15000.
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_{22}$ (N)'

elif cv_run in [ 's12', 'sig12', 'sig12f', 'sig12T' ] or cv_run[0:5]=='sig12' :
    cfield = 'sig12'
    tmin=-30000. ;  tmax=30000.   ;  df = 15000.
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_{12}$ (N)'

elif cv_run in [ 'Uice', 'Ut', 'u_ice', 'uVice' ]:
    cfield = 'Uice'
    tmin=-0.3 ;  tmax=0.3   ;  df = 0.2
    cpal_fld = 'RdBu_r'
    cunit = r'$u_{ice}$ (m/s)'

elif cv_run in [ 'Vice', 'Vt', 'v_ice', 'vUice' ]:
    cfield = 'Vice'
    tmin=-0.3 ;  tmax=0.3   ;  df = 0.2
    cpal_fld = 'RdBu_r'
    cunit = r'$v_{ice}$ (m/s)'

elif cv_run in [ 'zsN', 'MC', 'sigI', 'sigII' ]:
    cfield = 'sigma'
    tmin=-30000. ;  tmax=30000.   ;  df = 15000.
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_N$ (N)'

elif cv_run == 'zsS':
    cfield = 'sigS'
    tmin=-30000. ;  tmax=30000.   ;  df = 15000.
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_S$ (N)'
elif cv_run == 'zsS2':
    cfield = 'sigS'
    tmin=-0.5E9 ;  tmax=-tmin   ;  df = 1.E8
    cpal_fld = 'RdBu_r'
    cunit = r'$\sigma_S$ (N)'

elif cv_run in ['e11','e22','e12','div','ten','shr']:
    cfield = 'eps'
    tmin=-1.e-5 ;  tmax=-tmin   ;  df = 5.E-6
    cpal_fld = 'RdBu_r'
    cunit = r'$\epsilon$ (m/s^2)'

elif cv_run in ['shear2']:
    cfield = 'shear2'
    tmin=-1.e-10 ;  tmax=-tmin   ;  df = 5.E-9
    cpal_fld = 'RdBu_r'
    cunit = r'$\tau$^2'

elif cv_run in ['elasticity','elasticity-t','elasticity-f']:
    cfield = 'E'
    tmin=1.E8 ;  tmax=6.E8   ;  df = 1.E8
    cpal_fld = 'viridis'
    cunit = r'Elasticity'

elif cv_run in ['lambda-t','lambda-f']:
    cfield = 'lambda'
    tmin=1.E6 ;  tmax=1.E7   ;  df = 2.E6
    cpal_fld = 'viridis'
    cunit = r'$\lambda$'

elif cv_run in ['Pmax-t','Pmax-f']:
    cfield = 'Pmax'
    tmin=0. ;  tmax=1.E5   ;  df = 2.E4
    cpal_fld = 'ncview_tofino'
    cunit = r'$P_{max}$'

elif cv_run in ['Ptilde-t','Ptilde-f']:
    cfield = 'Ptilde'
    tmin=0. ;  tmax=1.   ;  df = 0.1
    cpal_fld = 'ncview_tofino'
    cunit = r'$\tilde{P}$'

elif cv_run in ['mult-t','mult-f']:
    cfield = 'mult'
    tmin=0. ;  tmax=1.   ;  df = 0.1
    cpal_fld = 'ncview_tofino'
    cunit = r'multiplicator'

elif cv_run in ['dmg-inc-t','dmg-inc-f']:
    cfield = 'd_inc'
    tmin=-1.e-3 ;  tmax=-tmin   ;  df = 2.e-4
    cpal_fld = 'ncview_tofino'
    cunit = r'multiplicator'

elif cv_run in ['dcrit-t','dcrit-f']:
    cfield = 'dcrit'
    tmin=-5 ;  tmax=5.   ;  df = 1.
    #cpal_fld = 'ncview_tofino'
    cpal_fld = 'RdBu_r'
    cunit = r'$d_{crit}$'

elif cv_run in ['Tdc-t','Tdc-f']:
    cfield = 'Tdc'
    tmin=0. ;  tmax=300.   ;  df = 50.
    cpal_fld = 'ncview_tofino'
    cunit = r'$Td_{c}$'

elif cv_run in ['mult-dmg-t','mult-dmg-f']:
    cfield = 'mult-dmg'
    tmin=0. ;  tmax=0.2   ;  df = 0.05
    cpal_fld = 'ncview_tofino'
    cunit = r'$MultDmg$'


    
else:
    print('ERROR: variable '+cv_run+' is not known yet...'); sys.exit(0)




# Time record stuff...
cp.chck4f(cf_run)
id_fld = Dataset(cf_run)
list_var = id_fld.variables.keys()
if 'time_counter' in list_var:    
    vtime = id_fld.variables['time_counter'][:]
    Nt = len(vtime)
    print('\n There is a "time_counter" in file '+cf_run+' !')
    print('   => '+str(Nt)+' snapshots!')
else:
    print('\nWARNING: there is NO "time_counter" in file '+cf_run+' !')
    print('   ==> setting Nt = 0 !\n')
    Nt = 0
id_fld.close()


ibath=1

cp.chck4f(cf_lsm)
cnmsk = 'tmask'
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

print('\n The shape of the domain is Ni, Nj =', Ni, Nj)


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

pmsk    = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)
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

fsize  = ( rh, rh*yx_ratio )
vc_fld = nmp.arange(tmin, tmax + df, df)





print('\n *** Opening file '+cf_run)
id_fld = Dataset(cf_run)


# Loop over time records:

for jt in range(Nt):

    ct   = '%4.4i'%(jt+1)
    cfig = cnfig+'_'+ct+'.'+fig_type

    fig = plt.figure(num = 1, figsize=fsize, dpi=rDPI, facecolor='k', edgecolor='k')

    if l_scientific_mode:
        ax  = plt.axes([0.09, 0.09, 0.9, 0.9], facecolor = 'r')
    else:
        ax  = plt.axes([0., 0., 1., 1.],     facecolor = '0.4')



    print('    => Reading record #'+str(jt)+' of '+cv_run+' in '+cf_run)
    XFLD  = id_fld.variables[cv_run][jt-1,j1:j2,i1:i2] ; # t, y, x
    print('          Done!\n')


    if XMSK[:,:].shape != XFLD.shape:
        print('\n PROBLEM: field and mask do not agree in shape!')
        print(XMSK.shape , XFLD.shape)
        sys.exit(0)
    print('  *** Shape of field and mask => ', nmp.shape(XFLD))
    
    
    
    l_add_true_filled = False
    
    if cfield == 'Bathymetry':
        (idy_nan,idx_nan) = nmp.where( nmp.isnan(XFLD) )
        #
        # LSM with different masking for true lsm and filled lsm...
        cf_mask_lbc = dir_conf+'/lsm_LBC_'+CCONF+'.nc'
        if path.exists(cf_mask_lbc):
            print('\n *** '+cf_mask_lbc+' found !!!')
            l_add_true_filled = True 
            id_filled = Dataset(cf_mask_lbc)
            xtmp = id_filled.variables['lsm'][j1:j2,i1:i2]
            id_filled.close()
            pfilled = nmp.ma.masked_where(xtmp[:,:] != -1., xtmp[:,:]*0.+40.)
            del xtmp
    
            print('  => done filling "pfilled" !\n')
    
        
    if cv_run == 'track':
        
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
    
    
    #cm = plt.imshow(pmsk, cmap=pal_lsm, norm=norm_lsm, interpolation='none')
    cm = plt.pcolormesh( pmsk, cmap=pal_lsm, norm=norm_lsm )
    
    if cfield == 'Bathymetry' and l_add_true_filled:
        # Ocean that has been filled turns black:
        cfl = plt.imshow(pfilled, cmap=pal_filled, norm=norm_filled, interpolation='none' ) #, interpolation='none')
    

    if l_add_true_filled: del pfilled
    
    
    plt.axis([ 0, Ni, 0, Nj])
    
    if l_scientific_mode:
        plt.xlabel('i-points', **cfont_axis)
        plt.ylabel('j-points', **cfont_axis)
    
    if l_show_cb:
    
        ax2 = plt.axes(vcb)
    
        if l_pow_field or l_log_field:
            clb = mpl.colorbar.ColorbarBase(ax=ax2,               cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend='neither')
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
        clb.set_label(cunit, **cfont_clb)
        clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color    
        #clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor
        if l_hide_cb_ticks: clb.ax.tick_params(axis=u'both', which=u'both',length=0) ; # remove ticks!
        plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_top) ; # set colorbar ticklabels
            
        del cf
        
    if l_show_nm:  ax.annotate(CCONF, xy=(1, 4), xytext=(x_cnf, y_cnf), **cfont_cnfn)
    
    if l_show_ttl: ax.annotate(CCONF, xy=(1, 4), xytext=(x_ttl, y_ttl), **cfont_ttl)
    
    
    
    #plt.savefig(cfig, dpi=rDPI, orientation='portrait', facecolor='b', transparent=True)
    plt.savefig(cfig, dpi=rDPI, orientation='portrait') #, transparent=True)
    print(cfig+' created!\n')
    plt.close(1)



    del cm, fig, ax
# END OF LOOP !!!

id_fld.close()

#!/usr/bin/env python

#       B a r a K u d a
#
#  Prepare 2D maps (monthly) that will later become a movie!
#  NEMO output and observations needed
#
#    L. Brodeau, January 2019

import sys
from os import path, getcwd
from string import replace
import argparse as ap
import numpy as nmp

from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.image as image
import matplotlib.cbook as cbook

import warnings
warnings.filterwarnings("ignore")

from calendar import isleap
import datetime

from re import split

import barakuda_colmap as bcm

import barakuda_tool as bt
import barakuda_ncio as bnc


cwd = getcwd()

l_smooth = False

vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

fig_type='png'
dpi = 110
color_top = 'white'
#color_top = 'k'

jt0 = 0


i2=0
j2=0
l_get_name_of_run = False
l_show_lsm = True
l_do_ice  = False
l_show_cb = True
l_log_field = False
l_pow_field = False
l_annotate_name = True
l_show_clock = True


cdir_logos = cwd+'/logos'
l_add_logo = True
cf_logo_on  = cdir_logos+'/ocean-next_trans_white_281x191.png'
l_add_logo_ige = True
cf_logo_ige = cdir_logos+'/IGE_blanc_notext.png'
l_add_logo_prace = True
cf_logo_prace = cdir_logos+'/PRACE_blanc.png'

rof_log = 150.



l_save_nc = False ; # save the field we built in a netcdf file !!!

l_apply_lap   = False
l_apply_hgrad = False

cb_extend = 'both' ;#colorbar extrema


################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--fin' , required=True,                help='specify the NEMO netCDF file to read from...')
requiredNamed.add_argument('-w', '--what', required=True, default="SST", help='specify the field/diagnostic to plot (ex: SST)')

parser.add_argument('-C', '--conf', default="eNATL60",      help='specify NEMO config (ex: eNATL60)')
parser.add_argument('-b', '--box' , default="ALL",          help='specify extraction box name (ex: ALL)')
parser.add_argument('-m', '--fmm' , default="mesh_mask.nc", help='specify the NEMO mesh_mask file (ex: mesh_mask.nc)')
parser.add_argument('-s', '--sd0' , default="20090101",     help='specify initial date as <YYYYMMDD>')
parser.add_argument('-l', '--lev' , type=int, default=0,    help='specify the level to use if 3D field (default: 0 => 2D)')
parser.add_argument('-z', '--zld' ,                         help='specify the topography netCDF file to use (field="z")')

args = parser.parse_args()

CNEMO = args.conf
CBOX  = args.box
CWHAT = args.what
cf_in = args.fin
cf_mm = args.fmm
csd0  = args.sd0
jk    = args.lev
cf_topo_land = args.zld

print ''
print ' *** CNEMO = ', CNEMO
print ' *** CBOX  = ', CBOX
print ' *** CWHAT = ', CWHAT
print ' *** cf_in = ', cf_in
print ' *** cf_mm = ', cf_mm
print ' *** csd0 = ', csd0
print ' ***   jk  = ', jk
l_add_topo_land = False
if args.zld != None:
    print ' *** cf_topo_land = ', cf_topo_land
    l_add_topo_land = True
l3d = False
if jk > 0: l3d=True
###############################################################################################################################################







x_logo  = 50 ; y_logo  = 50

cv_msk = 'tmask'


if   CWHAT == 'MLD':
    cv_in = 'somxl010'
    tmin=0. ;  tmax=1800.  ;  df = 50. ; cpal_fld = 'ncview_hotres' ;     cb_jump = 4
    cunit = r'MLD (m)'
    l_show_cb = True

elif CWHAT == 'SST':
    cv_in = 'sosstsst' ; #in ['sosstsst','tos']:    
    tmin=-2 ;  tmax=32.   ;  df = 1. ; cpal_fld = 'ncview_nrl' ;     cb_jump = 2
    #tmin=0. ;  tmax=32.   ;  df = 2. ; cpal_fld = 'viridis'
    #tmin=4. ;  tmax=20.   ;  df = 1. ; cpal_fld = 'PuBu'
    cunit = r'SST ($^{\circ}$C)'
    l_show_cb = True

elif CWHAT == 'SSS':
    cv_in = 'sosaline' ; #in ['sosstsst','tos']:    
    tmin=20. ;  tmax=40.   ;  df = 2. ; cpal_fld = 'ncview_ssec' ;     cb_jump = 2
    #tmin=0. ;  tmax=32.   ;  df = 2. ; cpal_fld = 'viridis'
    #tmin=4. ;  tmax=20.   ;  df = 1. ; cpal_fld = 'PuBu'
    cunit = r'Sea surface salinity'
    l_show_cb = True
    #l_save_nc = True

elif CWHAT == 'S':
    cv_in = 'vosaline' ; #in ['sosstsst','tos']:    
    #tmin=32. ;  tmax=38.   ;  df = 1. ; cpal_fld = 'ncview_ssec' ;     cb_jump = 1
    tmin=33.5 ;  tmax=36.5   ;  df = 0.5 ; cpal_fld = 'ncview_helix2' ; cb_jump = 1
    cunit = r'Salinity'
    l_show_cb = True
    l_save_nc = False

elif CWHAT == 'GRAD_SST':
    cv_in = 'sosstsst'
    l_apply_hgrad = True
    l_save_nc = True
    l_smooth = True ; nb_smooth  = 5
    tmin=0. ;  tmax=0.001 ;  df = 0.0001 ; cpal_fld = 'ncview_hotres' ; cb_jump = 1
    #tmin=0. ;  tmax=32.   ;  df = 2. ; cpal_fld = 'viridis'
    #tmin=4. ;  tmax=20.   ;  df = 1. ; cpal_fld = 'PuBu'
    cunit = r'$\left|\vec{\nabla}SST\right|$ (K/m)'
    l_show_cb = True

elif CWHAT == 'SSU':
    cv_in = 'sozocrtx'
    tmin=-0.8 ;  tmax=-tmin  ;  df = 0.1 ; cpal_fld = 'bone' ; cb_jump = 1
    cunit = r'Zonal component of current speed [m/s]'
    cv_msk = 'umask' ; l_show_cb = True

elif CWHAT == 'SSV':
    cv_in = 'somecrty'
    tmin=-1. ;  tmax=-tmin  ;  df = 0.2 ; cpal_fld = 'bone' ; cb_jump = 1
    cunit = r'Meridional component of current speed [m/s]'
    cv_msk = 'vmask' ; l_show_cb = True


elif CWHAT == 'SSH':
    cv_in = 'sossheig'    
    #cpal_fld = 'ncview_jaisnc'
    #cpal_fld = 'PuBu'
    #cpal_fld = 'RdBu'
    #cpal_fld = 'BrBG'
    #cpal_fld = 'bone_r' ; tmin=-0.9 ;  tmax=-tmin   ;  df = 0.05 ; l_apply_lap = True ; l_pow_field = True ; pow_field = 2.
    cpal_fld = 'RdBu_r' ; tmin=-3. ;  tmax=-tmin   ;  df = 0.5 ;
    #cpal_fld = 'RdBu_r' ; tmin=-1.5 ;  tmax=-tmin   ;  df = 0.5 ; 
    l_show_cb = True ; cb_jump = 1 ; x_logo = 2190 ; y_logo  = 1230
    cunit = r'SSH (m)'

elif CWHAT == 'LAP_SSH':
    cv_in = 'sossheig'
    l_apply_lap = True
    #cpal_fld = 'on3' ; tmin=-1.2 ;  tmax=2.3   ;  df = 0.05 ; l_apply_lap = True
    cpal_fld = 'on2' ; tmin=-1.2 ;  tmax=1.2   ;  df = 0.05 ; 
    #cpal_fld = 'coolwarm' ; tmin=-1. ;  tmax=1.   ;  df = 0.05 ; l_apply_lap = True
    #cpal_fld = 'gray_r' ; tmin=-0.3 ;  tmax=0.3   ;  df = 0.05 ; l_apply_lap = True

    

elif CWHAT == 'Amplitude':
    cv_in = 'r'    
    cpal_fld = 'RdBu_r' ; tmin=-0.5 ;  tmax=-tmin   ;  df = 0.1 ; cb_jump = 1
    cunit = r'Amplitude (m)'

elif CWHAT == 'Phase':
    cv_in = 'phi'    
    cpal_fld = 'RdBu_r' ; tmin=-30. ;  tmax=-tmin   ;  df = 5. ; cb_jump = 1
    #
    cunit = r'Phase (deg.)'
    
else:
    print 'ERROR: we do not know variable "'+str(cv_in)+'" !'
    sys.exit(0)


    
if CNEMO == 'eNATL60':

    # Defaults:
    Ni0 = 8354
    Nj0 = 4729
    l_show_clock = True
    x_clock = 1600 ; y_clock = 200 ; x_logo = 2200 ; y_logo  = 50
    cdt = '1h'
    l_get_name_of_run = True

    # Boxes:
    if   CBOX == 'ALL':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1440./float(Nj0) ; vcb=[0.61, 0.1, 0.36, 0.018]  ; font_rat=8.*rfact_zoom
        x_clock = 1600 ; y_clock = 200 ; x_logo=2200 ; y_logo=1200

    elif CBOX == 'EUROPA':
        i2=6400 ; j2=4000 ; i1=i2-2*1920; j1=j2-2*1080; rfact_zoom=0.5   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16 ; l_annotate_name=False
        if CWHAT == 'SST': tmin=0. ;  tmax=25.   ;  df = 1. ; cpal_fld = 'ncview_nrl'
        
    elif CBOX == 'EUROPAs':
        i2=6500 ; j2=3600 ; i1=i2-1920; j1=j2-1080; rfact_zoom=1.  ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16 ; l_annotate_name=False
        if CWHAT == 'SST': tmin=6. ;  tmax=18.   ;  df = 1. ; cpal_fld = 'ncview_nrl'
        
    elif   CBOX == 'ALLFR':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.39, 0.018]  ; font_rat=8.*rfact_zoom
        x_clock = 4000 ; y_clock = 200 ; x_logo = 6000 ; y_logo  = 50; l_show_clock=False ; l_annotate_name=False; l_add_logo=False

    elif   CBOX == 'SALL':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1080./float(Nj0) ; vcb=[0.59, 0.1, 0.39, 0.018]  ; font_rat=8.*rfact_zoom
        x_clock = 1600 ; y_clock = 200 ; x_logo = 2200 ; y_logo  = 50

    elif CBOX == 'FullMed':
        i1=5520; j1=1525; i2=i1+2560 ; j2=j1+1440 ; rfact_zoom=1. ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        l_annotate_name=False

    elif CBOX == 'Med+BS':
        i1=5400; j1=1530; i2=Ni0 ; j2=3310 ; rfact_zoom=1440./float(j2-j1)   ; vcb=[0.025, 0.06, 0.4, 0.02] ; font_rat=3.*rfact_zoom
        l_annotate_name=False ; l_add_logo_ige=False ; l_add_logo_prace=False
        x_clock = 100 ; y_clock = 170 ; x_logo = 2090 ; y_logo  = 10
        if CWHAT == 'SST': tmin=7. ;  tmax=25.   ;  df = 1. ; cpal_fld = 'ncview_nrl'

    elif CBOX == 'LabSea':
        # 1818,3600 -> 3630,4722
        #i2=3770 ; j2=Nj0 ; i1=i2-1920; j1=j2-1200;  rfact_zoom=1200./float(j2-j1)   ; vcb=[0.015, 0.07, 0.28, 0.02] ; font_rat=3.*rfact_zoom
        i2=3770 ; j2=Nj0 ; i1=i2-1920; j1=j2-1080;  rfact_zoom=1080./float(j2-j1)   ; vcb=[0.014, 0.09, 0.25, 0.02] ; font_rat=3.*rfact_zoom
        l_annotate_name=False ; l_add_logo_ige=False ; l_add_logo_prace=False
        x_clock = 40 ; y_clock = 170 ; x_logo = 2090 ; y_logo  = 10
        if CWHAT == 'SST': tmin=-2. ;  tmax=15.;  df = 1. ; cpal_fld = 'ncview_nrl'
        if CWHAT == 'SSS': tmin=28. ; tmax=35.5 ;  df = 1. ; cb_jump = 1 ; cpal_fld = 'gist_stern_r'
        l_save_nc = False

    elif CBOX == 'BlackSea':
        i1=Ni0-1920; j1=3330-1080; i2=Ni0 ; j2=3330 ; rfact_zoom=1.   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16 ; l_annotate_name=False

    elif CBOX == 'Brittany':
        i1=5400; j1=2850; i2=5700 ; j2=3100 ; rfact_zoom=4.   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=1.*rfact_zoom
        x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16 ; l_annotate_name=False
        if CWHAT == 'SST': tmin=7. ;  tmax=13.   ;  df = 1. ; cpal_fld = 'ncview_nrl'

    elif CBOX == 'Portrait':
        i1=2760; j1=1000; i2=4870; j2=4000 ; rfact_zoom=1.     ; vcb=[0.59, 0.1, 0.38, 0.018]  ; font_rat=1.*rfact_zoom
        l_annotate_name=False; l_show_clock=False

    elif CBOX == 'EATL':
        i1=3100; j1=2160; i2=i1+1920; j2=j1+1200 ; rfact_zoom=1. ; vcb=[0.3, 0.06, 0.38, 0.018] ; font_rat = 2. ; l_annotate_name=False
        x_clock = 1570 ; y_clock = 1150 ; x_logo = 1620 ; y_logo = 16

    elif CBOX == 'EATL2':
        i1=2740; j1=1600; i2=i1+2560; j2=j1+1440 ; rfact_zoom=1. ; vcb=[0.3, 0.06, 0.38, 0.018] ; font_rat = 2.5 ; l_annotate_name=False
        x_clock = 2140 ; y_clock = 1390 ; x_logo = 2220 ; y_logo = 16
        
    elif CBOX == 'GrlIcl':
        i1=3370; j1=3941; i2=5062; j2=Nj0 ; rfact_zoom=1. ; vcb=[0.3, 0.1, 0.38, 0.018] ; font_rat = 2. ; l_annotate_name=False
        x_clock = 1350 ; y_clock = 750 ; x_logo = 1400 ; y_logo = 16 ; l_save_nc=True
        if CWHAT == 'SST': tmax = 12.

    elif CBOX == 'Azores':
        # 785 x 1190 => comparison of two => 1600 x 1200 (5px for frame)
        i2=4410; j2=2240 ; i1=i2-785; j1=j2-1190; rfact_zoom=1. ; vcb=[0.05, 0.05, 0.9, 0.018] ; font_rat = 2. ; l_annotate_name=False
        l_add_logo=False; l_add_logo_prace=False; l_add_logo_ige=False
        x_clock = 400 ; y_clock = 120 ; l_save_nc=False
        if CWHAT == 'SST': tmin = 15. ; tmax = 25. ; df=0.5 ; cb_jump = 2

    elif CBOX == 'Band':
        i1=5100-1920; j1=2200; i2=5100; j2=j1+1080 ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.38, 0.018] ; font_rat = 2.           ; l_annotate_name=False
        l_show_clock = False ; l_add_logo = False ; #x_clock = 1420 ; y_clock = 1030 ; x_logo = 1500 ; y_logo = 16

    elif CBOX == 'Balear':
        i1=5750; j1=1880; i2=6470; j2=2600 ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.38, 0.018] ; font_rat = 2. ; l_annotate_name=False
        x_clock = 1420 ; y_clock = 1030 ; x_logo = 1500 ; y_logo = 16

    else:
        print ' ERROR: unknow box "'+CBOX+'" for config "'+CNEMO+'" !!!'
        sys.exit(0)



elif CNEMO == 'NATL60':
    Ni0 = 5422
    Nj0 = 3454
    #l_pow_field = True ; pow_field = 1.5
    l_show_clock = False
    cdt = '1h'
    #CBOX = 'zoom1' ; i1 = 1800 ; j1 = 950 ; i2 = i1+1920 ; j2 = j1+1080 ; rfact_zoom = 1. ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom ; l_show_lsm = False
    #CBOX = 'zoom1' ; i1 = 1800 ; j1 = 950 ; i2 = i1+2560 ; j2 = j1+1440 ; rfact_zoom = 1. ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom
    CBOX = 'ALL' ; i1=0 ; j1=0 ; i2=Ni0 ; j2=Nj0 ; rfact_zoom = 0.4 ; vcb=[0.59, 0.1, 0.38, 0.018] ; font_rat = 4.*rfact_zoom
    x_clock = 350 ; y_clock = 7 ; # where to put the date


elif CNEMO == 'NANUK025':
    l_do_ice = True
    cdt = '3h'; CBOX = 'ALL' ; i1 = 0 ; j1 = 0 ; i2 = 492 ; j2 = 614 ; rfact_zoom = 2. ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom
    x_clock = 350 ; y_clock = 7 ; # where to put the date




elif CNEMO == 'eNATL4':
    # Defaults:
    Ni0 = 559
    Nj0 = 313
    l_annotate_name = False
    l_show_clock = False
    l_add_logo = False
    x_clock = 1600 ; y_clock = 200 ; x_logo = 2200 ; y_logo  = 50
    cdt = '1h'

    # Boxes:
    if   CBOX == 'ALL':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=3. ; vcb=[0.59, 0.1, 0.39, 0.018]  ; font_rat=0.7*rfact_zoom
        l_show_cb = True
    else:
        print ' ERROR: unknow box "'+CBOX+'" for config "'+CNEMO+'" !!!'
        sys.exit(0)

else:
    print 'ERROR: we do not know NEMO config "'+str(CNEMO)+'" !'
    sys.exit(0)


CRUN = ''
if l_get_name_of_run:
    # Name of RUN:
    vv = split('-|_', path.basename(cf_in))
    if vv[0] != CNEMO:
        print 'ERROR: your file name is not consistent with "'+CNEMO+'" !!! ('+vv[0]+')' ; sys.exit(0)
    CRUN = vv[1]
    print '\n Run is called: "'+CRUN+'" !\n'


    

print '\n================================================================'
print '\n rfact_zoom = ', rfact_zoom
print ' font_rat = ', font_rat, '\n'


print ' i1,i2,j1,j2 =>', i1,i2,j1,j2

nx_res = i2-i1
ny_res = j2-j1

print ' *** nx_res, ny_res =', nx_res, ny_res

yx_ratio = float(nx_res+1)/float(ny_res+1)

rnxr = rfact_zoom*nx_res ; # widt image (in pixels)
rnyr = rfact_zoom*ny_res ; # height image (in pixels)

# Target resolution for figure:
rh_fig = round(rnyr/float(dpi),3) ; # width of figure
rw_fig = round(rh_fig*yx_ratio      ,3) ; # height of figure
rh_img = rh_fig*float(dpi)
rw_img = rw_fig*float(dpi)
while rw_img < round(rnxr,0):
    rw_fig = rw_fig + 0.01/float(dpi)
    rw_img = rw_fig*float(dpi)
while rh_img < round(rnyr,0):
    rh_fig = rh_fig + 0.01/float(dpi)
    rh_img = rh_fig*float(dpi)
    print ' *** size figure =>', rw_fig, rh_fig, '\n'
    print ' *** Forecasted dimension of image =>', rw_img, rh_img

print '\n================================================================\n\n\n'



cyr0=csd0[0:4]
cmn0=csd0[4:6]
cdd0=csd0[6:8]


l_3d_field = False



# Ice:
if l_do_ice:
    cv_ice  = 'siconc'
    cf_ice = replace(cf_in, 'grid_T', 'icemod')
    rmin_ice = 0.5
    cpal_ice = 'ncview_bw'
    vcont_ice = nmp.arange(rmin_ice, 1.05, 0.05)





if l_do_ice: bt.chck4f(cf_ice)

bt.chck4f(cf_mm)

l_notime=False
bt.chck4f(cf_in)
id_fld = Dataset(cf_in)
list_var = id_fld.variables.keys()
if 'time_counter' in list_var:
    vtime = id_fld.variables['time_counter'][:]
elif 'time' in list_var:
    vtime = id_fld.variables['time'][:]
else:
    l_notime=True
    print 'Did not find a time variable! Assuming no time and Nt=1'
id_fld.close()

Nt = 1
if not l_notime: Nt = len(vtime)




if l_show_lsm or l_apply_lap or l_apply_hgrad or l_add_topo_land:
    bt.chck4f(cf_mm)
    id_lsm = Dataset(cf_mm)
    nb_dim = len(id_lsm.variables[cv_msk].dimensions)
    print ' The mesh_mask has '+str(nb_dim)+' dimmensions!'
    print ' *** Reading '+cv_msk+' !'
    if nb_dim==4: XMSK = id_lsm.variables[cv_msk][0,jk,j1:j2,i1:i2]
    if nb_dim==3: XMSK = id_lsm.variables[cv_msk][jk,j1:j2,i1:i2]
    if nb_dim==2: XMSK = id_lsm.variables[cv_msk][j1:j2,i1:i2]
    (nj,ni) = nmp.shape(XMSK)

    if l_apply_lap:
        print ' *** Reading e1t and e2t !'
        XE1T2 = id_lsm.variables['e1t'][0,j1:j2,i1:i2]
        XE2T2 = id_lsm.variables['e2t'][0,j1:j2,i1:i2]
        XE1T2 = XE1T2*XE1T2
        XE2T2 = XE2T2*XE2T2
    if l_apply_hgrad:
        print ' *** Reading e1u and e2v !'
        XE1U = id_lsm.variables['e1u'][0,j1:j2,i1:i2]
        XE2V = id_lsm.variables['e2v'][0,j1:j2,i1:i2]
        print ' *** Reading umask and vmask !'
        if nb_dim==4:
            UMSK = id_lsm.variables['umask'][0,jk,j1:j2,i1:i2]
            VMSK = id_lsm.variables['vmask'][0,jk,j1:j2,i1:i2]
        if nb_dim==3:
            UMSK = id_lsm.variables['umask'][jk,j1:j2,i1:i2]
            VMSK = id_lsm.variables['vmask'][jk,j1:j2,i1:i2]
        if nb_dim==2:
            UMSK = id_lsm.variables['umask'][j1:j2,i1:i2]
            VMSK = id_lsm.variables['vmask'][j1:j2,i1:i2]
    if l_save_nc:
        Xlon = id_lsm.variables['glamt'][0,j1:j2,i1:i2]
        Xlat = id_lsm.variables['gphit'][0,j1:j2,i1:i2]
        
    if l3d:
        vdepth = id_lsm.variables['gdept_1d'][0,:]
        zdepth = vdepth[jk]
        cdepth = str(round(zdepth,1))+'m'
    id_lsm.close()

    print 'Shape Arrays => ni,nj =', ni,nj

    print 'Done!\n'



idx_land = nmp.where( XMSK < 0.01)

if l_add_topo_land:
    bt.chck4f(cf_topo_land)
    id_top = Dataset(cf_topo_land)
    print ' *** Reading "z" into:\n'+cf_topo_land
    xtopo = id_top.variables['z'][0,j1:j2,i1:i2]
    id_top.close()
    if nmp.shape(xtopo) != (nj,ni):
        print 'ERROR: topo and mask do not agree in shape!'; sys.exit(0)
    xtopo = xtopo*(1. - XMSK)
    bnc.dump_2d_field('topo_'+CBOX+'.nc', xtopo, name='z')    
    xtopo[nmp.where( XMSK > 0.01)] = nmp.nan


params = { 'font.family':'Helvetica Neue',
           'font.weight':    'normal',
           'font.size':       int(9.*font_rat),
           'legend.fontsize': int(9.*font_rat),
           'xtick.labelsize': int(9.*font_rat),
           'ytick.labelsize': int(9.*font_rat),
           'axes.labelsize':  int(9.*font_rat) }
mpl.rcParams.update(params)
cfont_clb  =  { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(7.*font_rat), 'color':color_top}
cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(9.*font_rat), 'color':'w' }
cfont_mail =  { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*font_rat), 'color':'0.8'}
cfont_titl =  { 'fontname':'Helvetica Neue', 'fontweight':'light', 'fontsize':int(30.*font_rat), 'color':'w' }


# Colormaps for fields:
pal_fld = bcm.chose_colmap(cpal_fld)
if   l_log_field:
    norm_fld = colors.LogNorm(                   vmin=tmin, vmax=tmax, clip=False)
elif l_pow_field:
    norm_fld = colors.PowerNorm(gamma=pow_field, vmin=tmin, vmax=tmax, clip=False)
else:
    norm_fld = colors.Normalize(                 vmin=tmin, vmax=tmax, clip=False)


if l_show_lsm or l_add_topo_land:
    if l_add_topo_land:
        xtopo = nmp.log10(xtopo+rof_log)
        pal_lsm = bcm.chose_colmap('gray_r')
        norm_lsm = colors.Normalize(vmin = nmp.log10(-100.+rof_log), vmax = nmp.log10(4000.+rof_log), clip = False)        
    else:
        pal_lsm = bcm.chose_colmap('land_dark')
        norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

if l_do_ice:
    pal_ice = bcm.chose_colmap(cpal_ice)
    norm_ice = colors.Normalize(vmin = rmin_ice, vmax = 1, clip = False)



if cdt == '3h':
    dt = 3
elif cdt == '1h':
    dt = 1
else:
    print 'ERROR: unknown dt!'


if l3d: cunit = cunit+' at '+cdepth
    
ntpd = 24/dt

vm = vmn
if isleap(int(cyr0)): vm = vml
#print ' year is ', vm, nmp.sum(vm)

jd = int(cdd0) - 1
jm = int(cmn0)

for jt in range(jt0,Nt):

    jh = (jt*dt)%24
    jdc = (jt*dt)/24 + 1

    if jt%ntpd == 0: jd = jd + 1

    if jd == vm[jm-1]+1 and (jt)%ntpd == 0 :
        jd = 1
        jm = jm + 1

    ch = '%2.2i'%(jh)
    #cdc= '%3.3i'%(jdc)
    cd = '%3.3i'%(jd)
    cm = '%2.2i'%(jm)

    #print '\n\n *** jt, ch, cd, cm =>', jt, ch, cd, cm


    ct = str(datetime.datetime.strptime(cyr0+'-'+cm+'-'+cd+' '+ch, '%Y-%m-%j %H'))
    ct=ct[:5]+cm+ct[7:] #lolo bug !!! need to do that to get the month and not "01"
    print ' ct = ', ct
    cday  = ct[:10]   ; print ' *** cday  :', cday
    chour = ct[11:13] ; print ' *** chour :', chour


    if l3d:
        cfig = 'figs/'+cv_in+'_'+CNEMO+'-'+CRUN+'_lev'+str(jk)+'_'+CBOX+'_'+cday+'_'+chour+'_'+cpal_fld+'.'+fig_type
    else:
        cfig = 'figs/'+cv_in+'_'+CNEMO+'-'+CRUN+'_'+CBOX+'_'+cday+'_'+chour+'_'+cpal_fld+'.'+fig_type

    ###### FIGURE ##############

    fig = plt.figure(num = 1, figsize=(rw_fig, rh_fig), dpi=None, facecolor='w', edgecolor='0.5')

    ax  = plt.axes([0., 0., 1., 1.], axisbg = '0.5')

    vc_fld = nmp.arange(tmin, tmax + df, df)


    
    print "Reading record #"+str(jt)+" of "+cv_in+" in "+cf_in
    if l3d: print '            => at level #'+str(jk)+' ('+cdepth+')!'
    id_fld = Dataset(cf_in)
    if l_notime:
        if l3d:
            Xplot  = id_fld.variables[cv_in][jk,j1:j2,i1:i2]
        else:
            Xplot  = id_fld.variables[cv_in][j1:j2,i1:i2]
    else:
        if l3d:
            Xplot  = id_fld.variables[cv_in][jt,jk,j1:j2,i1:i2] ; # t, y, x        
        else:
            Xplot  = id_fld.variables[cv_in][jt,j1:j2,i1:i2] ; # t, y, x        
    id_fld.close()
    print "Done!\n"

    if l_apply_lap:
        print ' *** Computing Laplacian of "'+cv_in+'"!'
        lx = nmp.zeros((nj,ni))
        ly = nmp.zeros((nj,ni))
        lx[:,1:ni-1] = 1.E9*(Xplot[:,2:ni] -2.*Xplot[:,1:ni-1] + Xplot[:,0:ni-2])/XE1T2[:,1:ni-1]
        ly[1:nj-1,:] = 1.E9*(Xplot[2:nj,:] -2.*Xplot[1:nj-1,:] + Xplot[0:nj-2,:])/XE2T2[1:nj-1,:]
        Xplot[:,:] = lx[:,:] + ly[:,:]
        del lx, ly

    if l_apply_hgrad:
        print ' *** Computing gradient of "'+cv_in+'"!'
        lx = nmp.zeros((nj,ni))
        ly = nmp.zeros((nj,ni))

        if l_smooth: bt.smoother(Xplot, XMSK, nb_smooth=nb_smooth)
        
        # Zonal gradient on T-points:
        lx[:,1:ni-1] = (Xplot[:,2:ni] - Xplot[:,0:ni-2]) / (XE1U[:,1:ni-1] + XE1U[:,0:ni-2]) * UMSK[:,1:ni-1] * UMSK[:,0:ni-2]
        lx[:,:] = XMSK[:,:]*lx[:,:]
        #bnc.dump_2d_field('dsst_dx_gridT.nc', lx, xlon=Xlon, xlat=Xlat, name='dsst_dx')
        # Meridional gradient on T-points:
        ly[1:nj-1,:] = (Xplot[2:nj,:] - Xplot[0:nj-2,:]) / (XE2V[1:nj-1,:] + XE2V[0:nj-2,:]) * VMSK[1:nj-1,:] * VMSK[0:nj-2,:]
        ly[:,:] = XMSK[:,:]*ly[:,:]
        #bnc.dump_2d_field('dsst_dy_gridT.nc', ly, xlon=Xlon, xlat=Xlat, name='dsst_dy')
        Xplot[:,:] = 0.0
        # Modulus of vector gradient:        
        Xplot[:,:] = nmp.sqrt(  lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] )
        #bnc.dump_2d_field('mod_grad_sst.nc', Xplot, xlon=Xlon, xlat=Xlat, name='dsst')
        del lx, ly


    print ''
    if not l_show_lsm and jt == jt0: ( nj , ni ) = nmp.shape(Xplot)
    print '  *** dimension of array => ', ni, nj, nmp.shape(Xplot)

    if l_save_nc:
        if l3d:
            cf_out = 'nc/'+CWHAT+'_NEMO_'+CNEMO+'_lev'+str(jk)+'_'+CBOX+'_'+cday+'_'+chour+'_'+cpal_fld+'.nc'
        else:
            cf_out = 'nc/'+CWHAT+'_NEMO_'+CNEMO+'_'+CBOX+'_'+cday+'_'+chour+'_'+cpal_fld+'.nc'
        print ' Saving in '+cf_out
        bnc.dump_2d_field(cf_out, Xplot, xlon=Xlon, xlat=Xlat, name=CWHAT)
        print ''


    

    print "Ploting"

    plt.axis([ 0, ni, 0, nj])

    Xplot[idx_land] = nmp.nan
    
    cf = plt.imshow(Xplot[:,:], cmap = pal_fld, norm = norm_fld, interpolation='none')
    del Xplot

    # Ice
    if not CWHAT == 'MLD' and l_do_ice:
        print "Reading record #"+str(jt)+" of "+cv_ice+" in "+cf_ice
        id_ice = Dataset(cf_ice)
        XICE  = id_ice.variables[cv_ice][jt,:,:] ; # t, y, x
        id_ice.close()
        print "Done!"

        #XM[:,:] = XMSK[:,:]
        #bt.drown(XICE, XM, k_ew=2, nb_max_inc=10, nb_smooth=10)
        #ci = plt.contourf(XICE[:,:], vcont_ice, cmap = pal_ice, norm = norm_ice) #

        pice = nmp.ma.masked_where(XICE < rmin_ice, XICE)
        ci = plt.imshow(pice, cmap = pal_ice, norm = norm_ice, interpolation='none') ; del pice, ci
        del XICE

    #LOLO: rm ???
    if l_show_lsm or l_add_topo_land:
        if l_add_topo_land:
            clsm = plt.imshow(nmp.ma.masked_where(XMSK>0.0001, xtopo), cmap = pal_lsm, norm = norm_lsm, interpolation='none')
            plt.contour(XMSK, [0.9], colors='k', linewidths=0.5)
        else:
            clsm = plt.imshow(nmp.ma.masked_where(XMSK>0.0001, XMSK), cmap = pal_lsm, norm = norm_lsm, interpolation='none')

    if l_show_cb:
        ax2 = plt.axes(vcb)
        clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend=cb_extend)
        cb_labs = []
        #if cb_jump > 1:
        cpt = 0
        for rr in vc_fld:
            if cpt % cb_jump == 0:
                if df >= 1.: cb_labs.append(str(int(rr)))
                if df <  1.: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./df)))+1) ))
            else:
                cb_labs.append(' ')
            cpt = cpt + 1
        #else:
        #    for rr in vc_fld: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./df)))+1) ))

        clb.ax.set_xticklabels(cb_labs, **cfont_clb)
        clb.set_label(cunit, **cfont_clb)
        clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color
        clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor
        clb.ax.tick_params(which = 'minor', length = 2, color = color_top )
        clb.ax.tick_params(which = 'major', length = 4, color = color_top )


    if l_show_clock:
        xl = float(x_clock)/rfact_zoom
        yl = float(y_clock)/rfact_zoom
        ax.annotate('Date: '+cday+' '+chour+':00', xy=(1, 4), xytext=(xl,yl), **cfont_clock)

    #ax.annotate('laurent.brodeau@ocean-next.fr', xy=(1, 4), xytext=(xl+150, 20), **cfont_mail)

    if l_annotate_name:
        xl = rnxr/20./rfact_zoom
        yl = rnyr/1.33/rfact_zoom
        ax.annotate(CNEMO, xy=(1, 4), xytext=(xl, yl), **cfont_titl)

    if l_add_logo:
        datafile = cbook.get_sample_data(cf_logo_on, asfileobj=False)
        im = image.imread(datafile)
        #im[:, :, -1] = 0.5  # set the alpha channel
        fig.figimage(im, x_logo, y_logo, zorder=9)
        del datafile, im
        #
    if l_add_logo_ige:
        datafile = cbook.get_sample_data(cf_logo_ige, asfileobj=False)
        im = image.imread(datafile)
        fig.figimage(im, x_logo+144, y_logo-150., zorder=9)
        del datafile, im
        #
    if l_add_logo_prace:
        datafile = cbook.get_sample_data(cf_logo_prace, asfileobj=False)
        im = image.imread(datafile)
        fig.figimage(im, x_logo-77, y_logo-140., zorder=9)
        del datafile, im

    plt.savefig(cfig, dpi=dpi, orientation='portrait', facecolor='k')
    print cfig+' created!\n'
    plt.close(1)

    if l_show_lsm: del clsm
    del cf, fig, ax
    if l_show_cb: del clb

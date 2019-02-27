#!/usr/bin/env python

#     CLIMPORN
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

vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

fig_type='png'
dpi = 110
color_top = 'white'
#color_top = 'k'

cv_out = 'unknown'

jt0 = 0


jk=0
i2=0
j2=0
l_get_name_of_run = False
l_show_lsm = True
l_show_cb = True
l_log_field = False
l_pow_field = False
l_annotate_name = True
l_show_clock = True
l_show_exp = False
l_fill_holes_black = False


cdir_logos = cwd+'/logos'
l_add_logo = True
cf_logo_on  = cdir_logos+'/ocean-next_trans_white_281x191.png'
l_add_logo_ige = True
cf_logo_ige = cdir_logos+'/IGE_blanc_notext.png'
l_add_logo_prace = True
cf_logo_prace = cdir_logos+'/PRACE_blanc.png'

rof_log = 150.
rof_dpt = 0.

l_3d_field = False

l_save_nc = False ; # save the field we built in a netcdf file !!!

romega = 2.*nmp.pi/86400.0

cb_extend = 'both' ;#colorbar extrema


################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-u', '--fiu' , required=True, help='specify the NEMO netCDF U file to read from...')
requiredNamed.add_argument('-v', '--fiv' , required=True, help='specify the NEMO netCDF V file to read from...')
requiredNamed.add_argument('-x', '--fldx', required=True, help='specify the name of the NEMO U field in U file')
requiredNamed.add_argument('-y', '--fldy', required=True, help='specify the name of the NEMO V field in V file')
requiredNamed.add_argument('-w', '--what', required=True, help='specify the field/diagnostic to plot (ex: CSPEED,CURLOF,ect.)')

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
cfx_in = args.fiu
cfy_in = args.fiv
cvx_in = args.fldx
cvy_in = args.fldy
cf_mm = args.fmm
csd0  = args.sd0
jk    = args.lev
cf_topo_land = args.zld

print ''
print ' *** CNEMO = ', CNEMO
print ' *** CBOX  = ', CBOX
print ' *** CWHAT = ', CWHAT
print ' *** cfx_in = ', cfx_in
print ' *** cvx_in = ', cvx_in
print ' *** cfx_in = ', cfy_in
print ' *** cvx_in = ', cvy_in
print ' *** cf_mm = ', cf_mm
print ' *** csd0 = ', csd0
print ' ***   jk  = ', jk
l_add_topo_land = False
if args.zld != None:
    print ' *** cf_topo_land = ', cf_topo_land
    l_add_topo_land = True
l3d = False
if jk > 0:
    l3d=True
else:
    jk=0
###############################################################################################################################################







x_logo  = 50 ; y_logo  = 50



l_do_crl=False
l_do_cof=False
l_do_cspd=False
l_do_tke=False
l_do_eke=False

if   CWHAT == 'CURL':
    l_do_crl = True  ; # do curl (relative-vorticity) !!!
    
elif CWHAT == 'CURLOF':
    l_do_cof = True  ; # do curl/f
    cpal_fld = 'on2' ; tmin=-1.2 ;  tmax=-tmin ;  df = 0.1 ; cb_jump = 2
    cunit = r'$\zeta/f$'
    
elif CWHAT == 'CURLOF_1000':
    l_do_cof = True  ; # do curl/f
    cpal_fld = 'on2' ; tmin=-0.6 ;  tmax=-tmin ;  df = 0.1 ; cb_jump = 1
    cunit = r'$\zeta/f$ at 1000 m'
    
elif CWHAT == 'CSPEED':
    l_do_cspd = True  ; # do current speed
    tmin=0. ; tmax=1.8 ; df = 0.2 ; cpal_fld = 'on3' ; # Poster full-res!!!
    cunit = 'Surface current speed [m/s]' ; cb_jump = 1
    cb_extend = 'max' ;#colorbar extrema

elif CWHAT == 'CSPEED_1000':
    l_do_cspd = True  ; # do current speed
    tmin=0. ; tmax=0.6 ; df = 0.1 ; cpal_fld = 'on3' ; # Poster full-res!!!
    cunit = 'Current speed at 1000 m [m/s]' ; cb_jump = 1
    cb_extend = 'max' ;#colorbar extrema

elif CWHAT == 'TKE':
    l_do_tke = True ; l_log_field = False
    tmin=0. ; tmax=3. ; df = 0.25 ; cpal_fld = 'ncview_hotres'
    #tmin=0. ; tmax=3. ; df = 0.25 ; cpal_fld = 'ncview_helix2'
    cunit = r'Turbulent Kinetic Energy [$m^2/s^2$]' ; cb_jump = 1
    cb_extend = 'max' ;#colorbar extrema

elif CWHAT == 'EKE':
    l_do_eke = True
    tmin=0. ; tmax=1. ; df = 0.1 ; cpal_fld = 'ncview_hotres' ; # Poster full-res!!!
    cunit = r'E Kinetic Energy [$m^2/s^2$]' ; cb_jump = 1
    cb_extend = 'max' ;#colorbar extrema
    
else:
    print ' ERROR: unknow diagnostic "'+CWHAT+'" !!!'
    sys.exit(0)
    
cv_out = CWHAT


    
if cvx_in=='vozocrtx' and cvy_in=='vomecrty': l_3d_field = True


elif cvx_in=='sozocrtx' and cvy_in=='somecrty' and l_do_cspd:
    # Current speed !
    if CBOX == 'Balear': tmin=0. ;  tmax=1.2 ;  df = 0.2 ; cpal_fld = 'ncview_hotres'
    #if CBOX == 'Balear': tmin=0. ;  tmax=1. ;  df = 0.2
    #if CBOX == 'Med+BS': tmin=0. ;  tmax=4. ;  df = 0.5 ; cpal_fld = 'ncview_hotres'
    if CBOX == 'Med+BS': tmin=0. ;  tmax=1.5 ;  df = 0.25 ; cpal_fld = 'on3'
        
elif cvx_in=='vozocrtx' and cvy_in=='vomecrty' and l_do_crl:
    l_3d_field = True
    #cpal_fld = 'on2' ; tmin=-0.025 ;  tmax=0.025 ;  df = 0.05
    #cpal_fld = 'ncview_bw' ; tmin=-0.025 ;  tmax=0.025 ;  df = 0.05
    #cpal_fld = 'gray' ; tmin=-0.025 ;  tmax=0.025 ;  df = 0.05
    cpal_fld = 'bone' ; tmin=-0.025 ;  tmax=0.025 ;  df = 0.05
    cunit = '';  cb_jump = 1
    l_show_clock = False
    l_add_logo   = False

#else:
#    print 'ERROR: we do not know cvx_in and cvy_in! ("'+cvx_in+'", "'+cvy_in+'")'
#    sys.exit(0)


if l3d and not l_3d_field:
    print 'ERROR: you cannot extract a level is the field is not 3D!!!'
    sys.exit(0)



if CNEMO == 'eNATL60':

    # Defaults:
    Ni0 = 8354
    Nj0 = 4729
    l_show_clock = True
    x_clock = 1600 ; y_clock = 200 ; x_logo = 2200 ; y_logo  = 50
    x_exp = 40 ; y_exp = 980
    cdt = '1h'
    l_get_name_of_run = True

    # Boxes:
    if   CBOX == 'ALL':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1440./float(Nj0) ; vcb=[0.61, 0.1, 0.36, 0.018]  ; font_rat=8.*rfact_zoom
        x_clock = 1600 ; y_clock = 200 ; x_logo=2200 ; y_logo=1200
        l_annotate_name=False ; l_show_exp = True ; x_exp = 1750 ; y_exp = 300
        l_fill_holes_black=True

    elif CBOX == 'EUROPA':
        i2=6400 ; j2=4000 ; i1=i2-2*1920; j1=j2-2*1080; rfact_zoom=0.5   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16 ; l_annotate_name=False
        
    elif CBOX == 'EUROPAs':
        i2=6500 ; j2=3600 ; i1=i2-1920; j1=j2-1080; rfact_zoom=1.  ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16 ; l_annotate_name=False
         
    elif   CBOX == 'ALLFR':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1. ; vcb=[0.6, 0.1, 0.38, 0.018]  ; font_rat=8.*rfact_zoom
        l_show_cb=True; l_add_logo=False ; l_add_logo_ige=False ; l_add_logo_prace=False ; l_show_clock=False ; l_annotate_name=False
        l_fill_holes_black=True

    elif   CBOX == 'SALL':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1080./float(Nj0) ; vcb=[0.59, 0.1, 0.39, 0.018]  ; font_rat=8.*rfact_zoom
        x_clock = 1600 ; y_clock = 200 ; x_logo = 2200 ; y_logo  = 50

    elif CBOX == 'Med':
        i2=8040; i1=i2-2560 ; j1=1525 ; j2=j1+1440 ; rfact_zoom=1. ; vcb=[0.025, 0.06, 0.4, 0.02] ; font_rat=2.*rfact_zoom
        l_annotate_name=False ; l_show_exp = True ; x_exp = 100 ; y_exp = 240
        x_clock = 200 ; y_clock = 170 ; x_logo = 1650 ; y_logo  = 1200
        if CWHAT == 'CURLOF': tmin=-1. ;  tmax=-tmin ; df = 0.1 ; cb_jump = 2
        if CWHAT == 'CSPEED': tmin= 0. ;  tmax=1.3   ; df = 0.1

    elif CBOX == 'Meddies':
        i2=5800; j1=1400; i1=i2-2560 ; j2=j1+1440 ; rfact_zoom=1. ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        l_annotate_name=False


    elif CBOX == 'Med+BS':
        #
        i2=Ni0 ; i1=5400; j1=1530;  j2=3310 ; rfact_zoom=1440./float(j2-j1)   ; vcb=[0.025, 0.06, 0.4, 0.02] ; font_rat=3.*rfact_zoom
        l_annotate_name=False ; l_add_logo_ige=False ; l_add_logo_prace=False
        x_clock = 100 ; y_clock = 170 ; x_logo = 2090 ; y_logo  = 10

    elif CBOX == 'BlackSea':
        i1=Ni0-1920; j1=3330-1080; i2=Ni0 ; j2=3330 ; rfact_zoom=1.   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16 ; l_annotate_name=False

    elif CBOX == 'Brittany':
        i1=5400; j1=2850; i2=5700 ; j2=3100 ; rfact_zoom=4.   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=1.*rfact_zoom
        x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16 ; l_annotate_name=False

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
        x_clock = 1350 ; y_clock = 750 ; x_logo = 1400 ; y_logo = 16

    elif CBOX == 'AzoresP':
        # Azores Portrait
        # 785 x 1190 => comparison of two => 1600 x 1200 (5px for frame), image is: 780x1190
        ## use: CMD="montage -geometry 780x1190+10+5 -background black <img1> <img2> <img_montage>"
        # => montage is then 1600x1200
        i2=4410; j2=2240 ; i1=i2-780; j1=j2-1190; rfact_zoom=1. ; vcb=[0.05, 0.05, 0.9, 0.018] ; font_rat = 2. ; l_annotate_name=False
        l_add_logo=False; l_add_logo_prace=False; l_add_logo_ige=False
        x_clock = 400 ; y_clock = 120
        if CWHAT=='CURLOF': tmin=-0.8 ;  tmax=-tmin ;  df = 0.1 ; cb_jump = 2

    elif CBOX == 'AzoresL':
        # Azores Landscape: (1920x1080)
        i2=5200; j2=2240 ; i1=i2-1920; j1=j2-1080; rfact_zoom=1. ; vcb=[0.57, 0.08, 0.4, 0.018] ; font_rat = 2. ; l_annotate_name=False
        l_add_logo=False; l_add_logo_prace=False; l_add_logo_ige=False
        x_clock = 1400 ; y_clock = 120 ; l_show_exp = True ; x_exp = 40 ; y_exp = 1030

    elif CBOX == 'AzoresS':
        # Azores small square:
        l_show_cb = True ; l_show_clock = False
        i2=4300; j2=2140 ; i1=i2-360; j1=j2-360; rfact_zoom=1. ; vcb=[0.05, 0.12, 0.9, 0.018] ; font_rat = 1.5 ; l_annotate_name=False
        l_add_logo=False; l_add_logo_prace=False; l_add_logo_ige=False
        x_clock = 400 ; y_clock = 120
        if CWHAT == 'SST': tmin = 15. ; tmax = 25. ; df=0.5 ; cb_jump = 2
        if CWHAT == 'W_1000': color_top = 'k'
        if CWHAT == 'CSPEED': tmin=0. ;  tmax=1.2   ;  df = 0.2 ; cb_jump = 1
        l_save_nc = True

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
    vv = split('-|_', path.basename(cfx_in))
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







bt.chck4f(cf_mm)
bt.chck4f(cfx_in)
bt.chck4f(cfy_in)

id_fx = Dataset(cfx_in)
vtime = id_fx.variables['time_counter'][:]
id_fx.close()

Nt = len(vtime)

if l_show_lsm or l_do_crl or l_do_cof or l_do_cspd or l_do_tke or l_do_eke:
    print "\nReading record metrics in "+cf_mm
    id_lsm = Dataset(cf_mm)
    nb_dim = len(id_lsm.variables['tmask'].dimensions)
    print ' The mesh_mask has '+str(nb_dim)+' dimmensions!'
    if l_show_lsm:
        cv_msk = 'tmask'
        if l_do_crl or l_do_cof:  cv_msk = 'fmask'
        print '\n Reading mask as "'+cv_msk+'" in:'; print cv_msk,'\n'
        if nb_dim==4: XMSK = id_lsm.variables[cv_msk][0,jk,j1:j2,i1:i2]
        if nb_dim==3: XMSK = id_lsm.variables[cv_msk][0,j1:j2,i1:i2]
        if nb_dim==2: XMSK = id_lsm.variables[cv_msk][j1:j2,i1:i2]
        (nj,ni) = nmp.shape(XMSK)
    if l_do_crl or l_do_cof:
        # e2v, e1u, e1f, e2f
        e2v = id_lsm.variables['e2v'][0,j1:j2,i1:i2]
        e1u = id_lsm.variables['e1u'][0,j1:j2,i1:i2]
        e1f = id_lsm.variables['e1f'][0,j1:j2,i1:i2]
        e2f = id_lsm.variables['e2f'][0,j1:j2,i1:i2]
    if l_do_cof:
        ## Coriolis Parameter:
        ff  = id_lsm.variables['gphif'][0,j1:j2,i1:i2]
        ff[:,:] = 2.*romega*nmp.sin(ff[:,:]*nmp.pi/180.0)
        (nj,ni) = nmp.shape(XMSK)
    if l_save_nc:
        Xlon = id_lsm.variables['glamt'][0,j1:j2,i1:i2]
        Xlat = id_lsm.variables['gphit'][0,j1:j2,i1:i2]
        
    if l3d:
        vdepth = id_lsm.variables['gdept_1d'][0,:]
        zdepth = vdepth[jk]
        if zdepth>990. and zdepth<1010.: rof_log = 1000.
        rof_dpt = zdepth
        cdepth = str(round(zdepth,1))+'m'
    id_lsm.close()

    print 'Shape Arrays => ni,nj =', ni,nj

    print 'Done!\n'



idx_land = nmp.where( XMSK <  0.5 )
#idx_ocea = nmp.where( XMSK >= 0.5 )

#(vjj_land, vji_land) = idx_land
#print idx_land
#print vjj_land
#print vji_land
#print len(vjj_land), len(vji_land)
#sys.exit(0)


XLSM = nmp.zeros((nj,ni)) ; # will define real continents not NEMO mask...

if l_add_topo_land:
    bt.chck4f(cf_topo_land)
    id_top = Dataset(cf_topo_land)
    print ' *** Reading "z" into:\n'+cf_topo_land
    xtopo = id_top.variables['z'][0,j1:j2,i1:i2]
    id_top.close()
    if nmp.shape(xtopo) != (nj,ni):
        print 'ERROR: topo and mask do not agree in shape!'; sys.exit(0)
    xtopo = xtopo*(1. - XMSK)
    #bnc.dump_2d_field('topo_'+CBOX+'.nc', xtopo, name='z')    
    if l3d: xtopo = xtopo + rof_dpt
    xtopo[nmp.where( XMSK > 0.01)] = nmp.nan
    if l_fill_holes_black and not l3d:
        XLSM[nmp.where( xtopo < 0.0)] = 1
        xtopo[nmp.where( xtopo < 0.0)] = nmp.nan

XLSM[nmp.where( XMSK > 0.5)] = 1

        
params = { 'font.family':'Helvetica Neue',
           'font.weight':    'normal',
           'font.size':       int(9.*font_rat),
           'legend.fontsize': int(9.*font_rat),
           'xtick.labelsize': int(9.*font_rat),
           'ytick.labelsize': int(9.*font_rat),
           'axes.labelsize':  int(9.*font_rat) }
mpl.rcParams.update(params)
cfont_clb  =  { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(7.*font_rat), 'color':color_top}
cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(9.*font_rat), 'color':color_top }
cfont_exp= { 'fontname':'Open Sans'  , 'fontweight':'light', 'fontsize':int(9.*font_rat), 'color':color_top }
cfont_mail =  { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*font_rat), 'color':'0.8'}
cfont_titl =  { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(30.*font_rat), 'color':color_top }


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
        #norm_lsm = colors.Normalize(vmin = nmp.log10(min(-100.+rof_dpt/3.,0.) + rof_log), vmax = nmp.log10(4000.+rof_dpt + rof_log), clip = False)
        norm_lsm = colors.Normalize(vmin = nmp.log10(-100. + rof_log), vmax = nmp.log10(4000.+rof_dpt + rof_log), clip = False)
    else:
        pal_lsm = bcm.chose_colmap('land_dark')
        norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)


if cdt == '3h':
    dt = 3
elif cdt == '1h':
    dt = 1
else:
    print 'ERROR: unknown dt!'










if l_do_eke:
    Umean = nmp.zeros((nj,ni))
    Vmean = nmp.zeros((nj,ni))
    # Need to go for the mean first!!!
    print '\n Goint for '+str(Nt)+' snaphots to compute the mean first!!!'
    for jt in range(jt0,Nt):
        print "Reading record #"+str(jt)+"..."

        id_fx = Dataset(cfx_in)
        if not l_3d_field:
            XFLD  = id_fx.variables[cvx_in][jt,j1:j2,i1:i2] ; # t, y, x
        else:
            XFLD  = id_fx.variables[cvx_in][jt,jk,j1:j2,i1:i2] ; # t, y, x
        id_fx.close()
        
        id_fy = Dataset(cfy_in)
        if not l_3d_field:
            YFLD  = id_fy.variables[cvy_in][jt,j1:j2,i1:i2] ; # t, y, x
        else:
            YFLD  = id_fy.variables[cvy_in][jt,jk,j1:j2,i1:i2] ; # t, y, x
        id_fy.close()
        
        print "Done!"
        # On T-points:
        Umean[:,2:ni] = Umean[:,2:ni] + 0.5*( XFLD[:,1:ni-1] + XFLD[:,2:ni] )/float(Nt)
        Vmean[2:nj,:] = Vmean[2:nj,:] + 0.5*( YFLD[1:nj-1,:] + YFLD[2:nj,:] )/float(Nt)
    print ' time averaging done...\n\n'



    
ntpd = 24/dt

vm = vmn
if isleap(int(cyr0)): vm = vml
#print ' year is ', vm, nmp.sum(vm)

jd = int(cdd0) - 1
jm = int(cmn0)








Xplot = nmp.zeros((nj,ni))


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
        cfig = 'figs/'+cv_out+'_'+CNEMO+'-'+CRUN+'_lev'+str(jk)+'_'+CBOX+'_'+cday+'_'+chour+'_'+cpal_fld+'.'+fig_type
    else:
        cfig = 'figs/'+cv_out+'_'+CNEMO+'-'+CRUN+'_'+CBOX+'_'+cday+'_'+chour+'_'+cpal_fld+'.'+fig_type

    ###### FIGURE ##############

    fig = plt.figure(num = 1, figsize=(rw_fig, rh_fig), dpi=None, facecolor='w', edgecolor='0.5')

    ax  = plt.axes([0., 0., 1., 1.], axisbg = 'k')

    vc_fld = nmp.arange(tmin, tmax + df, df)


    print "Reading record #"+str(jt)+" of "+cvx_in+" in "+cfx_in
    if l3d: print '            => at level #'+str(jk)+' ('+cdepth+')!'
    id_fx = Dataset(cfx_in)
    if not l_3d_field:
        XFLD  = id_fx.variables[cvx_in][jt,j1:j2,i1:i2] ; # t, y, x
    else:
        print 'j1:j2 =', j1,j2
        print 'i1:i2 =', i1,i2
        XFLD  = id_fx.variables[cvx_in][jt,jk,j1:j2,i1:i2] ; # t, y, x
    id_fx.close()
    print "Done!"

    print "Reading record #"+str(jt)+" of "+cvy_in+" in "+cfy_in
    if l3d: print '            => at level #'+str(jk)+' ('+cdepth+')!'
    id_fy = Dataset(cfy_in)
    if not l_3d_field:
        YFLD  = id_fy.variables[cvy_in][jt,j1:j2,i1:i2] ; # t, y, x
    else:
        YFLD  = id_fy.variables[cvy_in][jt,jk,j1:j2,i1:i2] ; # t, y, x
    id_fy.close()
    print "Done!"


    lx = nmp.zeros((nj,ni))
    ly = nmp.zeros((nj,ni))

    if l_do_crl or l_do_cof:
        print '\nComputing curl...'
        lx[:,1:ni-1] =   e2v[:,2:ni]*YFLD[:,2:ni] - e2v[:,1:ni-1]*YFLD[:,1:ni-1]
        ly[1:nj-1,:] = - e1u[2:nj,:]*XFLD[2:nj,:] + e1u[1:nj-1,:]*XFLD[1:nj-1,:]
        if l_do_cof: Xplot[:,:] = ( lx[:,:] + ly[:,:] )*XMSK[:,:] / ( e1f[:,:]*e2f[:,:]*ff[:,:] ) # Relative Vorticity...
        if l_do_crl: Xplot[:,:] = ( lx[:,:] + ly[:,:] )*XMSK[:,:] / ( e1f[:,:]*e2f[:,:] ) * 1000. # Curl...

    if l_do_cspd:
        print '\nComputing current speed at T-points ...'
        lx[:,2:ni] = 0.5*( XFLD[:,1:ni-1] + XFLD[:,2:ni] )
        ly[2:nj,:] = 0.5*( YFLD[1:nj-1,:] + YFLD[2:nj,:] )
        Xplot[:,:] = nmp.sqrt( lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] ) * XMSK[:,:]

    if l_do_tke:
        print '\nComputing TKE at T-points ...'
        lx[:,2:ni] = 0.5*( XFLD[:,1:ni-1] + XFLD[:,2:ni] )
        ly[2:nj,:] = 0.5*( YFLD[1:nj-1,:] + YFLD[2:nj,:] )
        Xplot[:,:] = ( lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] ) * XMSK[:,:]


    if l_do_eke:
        print '\nComputing TKE at T-points ...'
        lx[:,2:ni] = 0.5*( XFLD[:,1:ni-1] + XFLD[:,2:ni] ) - Umean[:,2:ni]
        ly[2:nj,:] = 0.5*( YFLD[1:nj-1,:] + YFLD[2:nj,:] ) - Vmean[2:nj,:]
        Xplot[:,:] = ( lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] ) * XMSK[:,:]

    del lx, ly
    print '... '+cv_out+' computed!\n'
        
    del XFLD,YFLD


    if l_save_nc:
        cf_out = 'nc/'+cv_out+'_NEMO_'+CNEMO+'-'+CRUN+'_'+CBOX+'_'+cday+'_'+chour+'_'+cpal_fld+'.nc'
        print ' Saving in '+cf_out
        bnc.dump_2d_field(cf_out, Xplot, xlon=Xlon, xlat=Xlat, name=cv_out)
        print ''


    

    print "Ploting"

    plt.axis([ 0, ni, 0, nj])

    Xplot[idx_land] = nmp.nan
    
    cf = plt.imshow(Xplot[:,:], cmap = pal_fld, norm = norm_fld, interpolation='none')


    #LOLO: rm ???
    if l_show_lsm or l_add_topo_land:
        if l_add_topo_land:
            clsm = plt.imshow(nmp.ma.masked_where(XLSM>0.0001, xtopo), cmap = pal_lsm, norm = norm_lsm, interpolation='none')
            plt.contour(XLSM, [0.9], colors='k', linewidths=0.5)
        else:
            clsm = plt.imshow(nmp.ma.masked_where(XLSM>0.0001, XLSM), cmap = pal_lsm, norm = norm_lsm, interpolation='none')

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

    if l_get_name_of_run and l_show_exp:
        xl = float(x_exp)/rfact_zoom
        yl = float(y_exp)/rfact_zoom
        ax.annotate('Experiment: '+CNEMO+'-'+CRUN, xy=(1, 4), xytext=(xl,yl), **cfont_exp)


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

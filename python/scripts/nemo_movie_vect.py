#!/usr/bin/env python3
#
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#     CLIMPORN
#
#  Prepare 2D maps (monthly) that will later become a movie!
#  NEMO output and observations needed
#
#    L. Brodeau, November 2019
##################################################################

import sys
from os import path, getcwd, makedirs
import argparse as ap
import numpy as np

from netCDF4 import Dataset, num2date

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.image as image
import matplotlib.cbook as cbook

from calendar import isleap
import datetime

from re import split

import climporn as cp
from climporn import fig_style_mov as fsm

import warnings
warnings.filterwarnings("ignore")


cwd = getcwd()

l_smooth = False

vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

fig_type='png'
rDPI = 240
color_top = 'white'
color_top_cb = 'white'
#color_top = 'k'

cv_out = 'unknown'

jt0 = 0

jk=0
l_show_lsm = True
l_log_field = False
l_pow_field = False


rof_log = 150.
rof_dpt = 0.

l_3d_field = False

l_save_nc = False ; # save the field we built in a netcdf file !!!

romega = 7.292115083046062E-005 # Coriolis [1/s] (same as in NEMO 3.6 / #romega = 2.*np.pi/86400.0)

cb_extend = 'both' ;#colorbar extrema

# Normally logos should be found there:
dir_logos = str.replace( path.dirname(path.realpath(__file__)) , 'climporn/python', 'climporn/misc/logos' )
print("\n --- logos found into : "+dir_logos+" !\n")


################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-u', '--fiu' , required=True, help='NEMO netCDF U file to read from...')
requiredNamed.add_argument('-v', '--fiv' , required=True, help='NEMO netCDF V file to read from...')
requiredNamed.add_argument('-x', '--fldx', required=True, help='name of the NEMO U field in U file')
requiredNamed.add_argument('-y', '--fldy', required=True, help='name of the NEMO V field in V file')
requiredNamed.add_argument('-w', '--what', required=True, help='field/diagnostic to plot (ex: CSPEED,CURLOF,ect.)')

parser.add_argument('-C', '--conf', default="none",           help='name of NEMO config (ex: eNATL60) (defined into `nemo_hboxes.py`)')
parser.add_argument('-E', '--nexp', default=None,             help='name of experiment (shows up in figure name)')
parser.add_argument('-b', '--box' , default="ALL",            help='extraction box name (ex: ALL) (defined into `nemo_hboxes.py`)')
parser.add_argument('-m', '--fmm' , default="mesh_mask.nc",   help='NEMO mesh_mask file (ex: mesh_mask.nc)')
#parser.add_argument('-s', '--sd0' , default=None,             help='initial date as <YYYY-MM-DD>')
parser.add_argument('-l', '--lev' , type=int, default=0,      help='level to use if 3D field (default: 0 => 2D)')
parser.add_argument('-z', '--zld' ,                           help='topography netCDF file to use (field="z")')
parser.add_argument('-t', '--tstep',  default="1h",           help='time step ("1h","2h",..,up to "1d") in input file')
parser.add_argument('-N', '--oname',  default="",             help='a name that overides `CONF` on the plot...')
parser.add_argument('-o', '--outdir', default="./figs",       help='path to directory where to save figures')
parser.add_argument('-f', '--fignm',  default="",             help='common string in name of figure to create')
parser.add_argument('-S', '--sign',   default="",             help='sign the image with some text')
parser.add_argument('-L', '--ltr' ,   default=None,           help='letter as in "a)" for figures in paper')
#
# Colorbar or not ?
parser.add_argument('--cb',    action='store_true')
parser.add_argument('--no-cb', dest='cb', action='store_false')
parser.set_defaults(cb=True)


args = parser.parse_args()

CNEMO = args.conf
CEXP  = args.nexp
CBOX  = args.box
CWHAT = args.what
cfx_in = args.fiu
cfy_in = args.fiv
cvx_in = args.fldx
cvy_in = args.fldy
cf_mm = args.fmm
#csd0  = args.sd0
jk    = args.lev
cf_topo_land = args.zld
cdt    = args.tstep  ; # time step, in the form "1h", "2h", ..., "12h", ..., "1m", ..., "6m", ... "1y", ..., etc
CONAME = args.oname
cd_out = args.outdir
cn_fig = args.fignm
CSIGN  = args.sign
lcb    = args.cb
cltr   = args.ltr

print('')
print(' *** CNEMO = ', CNEMO)
print(' *** CBOX  = ', CBOX)
print(' *** CWHAT = ', CWHAT)
print(' *** cfx_in = ', cfx_in)
print(' *** cvx_in = ', cvx_in)
print(' *** cfy_in = ', cfy_in)
print(' *** cvy_in = ', cvy_in)
print(' *** cf_mm = ', cf_mm)
print(' *** show colorbar:', lcb)


#lForceD0 = False
#if csd0:
#    lForceD0 = True
#    print(' *** csd0 = ', csd0)
#    yyyy0 = csd0[0:4]

print(' ***   jk  = ', jk)
if CONAME != "": print(' *** CONAME = ', CONAME)
l_add_topo_land = False
if args.zld != None:
    print(' *** cf_topo_land = ', cf_topo_land)
    l_add_topo_land = True
l3d = False
if jk > 0:
    l3d=True
else:
    jk=0

l_add_sign = ( CSIGN != '' )

##########################################################################################

cdir_figs = cd_out+'/'+CWHAT
makedirs( cdir_figs, exist_ok=True )

if l_save_nc and not path.exists('nc'): mkdir('nc')

CRUN = ''
if CNEMO == 'none':
    # Name of RUN:
    vv = split('-|_', path.basename(cfx_in))
    if vv[0] != CNEMO:
        print('ERROR: your file name is not consistent with '+CNEMO+' !!! ('+vv[0]+')'); sys.exit(0)
    CRUN = vv[1]
    print('\n Run is called: ''+CRUN+'' !\n')
    CRUN='-'+CRUN

#---------------------------------------------------------------

nemo_box = cp.nemo_hbox(CNEMO,CBOX)

(Ni0,Nj0) = nemo_box.size()
print(' '+CNEMO+': Ni0,Nj0 => ', Ni0,Nj0)

(i1,j1, i2,j2) = nemo_box.idx()
print(' i1,j1, i2,j2 => ', i1,j1, i2,j2,'\n')

if nemo_box.l_show_name:  (x_name,y_name)   = nemo_box.name
if nemo_box.l_show_clock: (x_clock,y_clock) = nemo_box.clock
if nemo_box.l_show_exp:   (x_exp,y_exp)     = nemo_box.exp
if nemo_box.l_add_logo:   (x_logo,y_logo)   = nemo_box.logo
if l_add_sign and nemo_box.l_show_sign:  (x_sign,y_sign)   = nemo_box.sign

#---------------------------------------------------------------




l_do_crl=False
l_do_cof=False
l_do_cspd=False
l_do_tke=False
l_do_eke=False

if   CWHAT == 'CURL':
    cpal_fld = 'on2' ; tmin=-1.2 ;  tmax=-tmin ;  df = 0.1 ; cb_jump = 2
    l_do_crl = True  ; # do curl (relative-vorticity) !!!

elif CWHAT == 'CURLOF':
    l_do_cof = True  ; # do curl/f
    cpal_fld = 'on2' ; tmin=-1.2 ;  tmax=-tmin ;  df = 0.1 ; cb_jump = 2
    cunit = r'$\zeta/f$'
    if CBOX ==     'ALLC': tmin=-1. ;  tmax=-tmin ; df = 0.1 ; cb_jump = 2 ; cpal_fld='RdBu_r'; #; cpal_fld = 'bone'
    #if CBOX ==     'ALLFR': tmin=-1. ;  tmax=-tmin ; df = 0.1 ; cb_jump = 2 ; cpal_fld='RdBu_r'; #; cpal_fld = 'bone'
    #if CBOX ==      'Med': tmin=-1. ;  tmax=-tmin ; df = 0.1 ; cb_jump = 2 ; cpal_fld='RdBu_r'
    if CBOX ==      'Med': tmin=-1. ;  tmax=-tmin ; df = 0.1 ; cb_jump = 2
    if CBOX in ['AzoresP','MeddiesW','ALLFR']: tmin=-0.8 ;  tmax=-tmin ;  df = 0.1 ; cb_jump = 1
    if CBOX == 'BlackSea': tmin=-0.6 ;  tmax=-tmin  ; df = 0.1 ; cb_jump = 1
    if CBOX ==  'EATLcom': tmin=-1. ;  tmax=-tmin ; df = 0.1 ; cb_jump = 2
    if CNEMO in ['TROPICO05_NST','CALEDO10']: tmin=-0.5 ;  tmax=-tmin ; df = 0.05 ; cb_jump = 2    
    if CNEMO in ['CALEDO60']: tmin=-1. ;  tmax=-tmin ; df = 0.1 ; cb_jump = 2    

elif CWHAT == 'CURLOF_1000':
    l_do_cof = True  ; # do curl/f
    cpal_fld = 'on2' ; tmin=-0.6 ;  tmax=-tmin ;  df = 0.1 ; cb_jump = 1
    cunit = r'$\zeta/f$ at 1000 m'
    
elif CWHAT == 'CSPEED':
    l_do_cspd = True  ; # do current speed
    tmin=0. ; tmax=1.8 ; df = 0.2 ; cpal_fld = 'on3' ; # Poster full-res!!!
    cunit = 'Surface current speed [m/s]' ; cb_jump = 1
    cb_extend = 'max' ;#colorbar extrema
    if CBOX ==    'Med'  : tmax=1.3 ; df = 0.1
    if CBOX ==  'AzoresS': tmax=1.2 ; df = 0.2 ; cb_jump = 1
    if CBOX == 'BlackSea': tmax=0.8 ; df = 0.1 ; cb_jump = 1
    if CBOX in ['Manche','Bretagne']:  tmax=2.4 ; df = 0.1 ; cb_jump = 2
    if CBOX in ['Brest']:  tmax=2.4 ; df = 0.1 ; cb_jump = 2
    if CNEMO in ['CALEDO60']: tmax=0.8 ; df = 0.1 ; cb_jump = 1

elif CWHAT == 'sivelo-t':
    l_do_cspd = True  ; # do current speed
    tmin=0. ; tmax=0.6 ; df=0.1 ; cb_extend = 'max' ; cb_jump = 1
    cpal_fld='cmocean_dense' ; color_top_cb='w'
    imask_no_ice_pc = 5 ; # we hide the field where A<10%                                                                                                                                   
    cunit = 'Sea-ice velocity [m/s]'

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
    print(' ERROR: unknow diagnostic ''+CWHAT+'' !!!')
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
#    print('ERROR: we do not know cvx_in and cvy_in! (''+cvx_in+'', ''+cvy_in+'')'
#    sys.exit(0)


if l3d and not l_3d_field:
    print('ERROR: you cannot extract a level is the field is not 3D!!!')
    sys.exit(0)



    
rfz   = nemo_box.rfact_zoom
fontr = nemo_box.font_rat

print('\n================================================================')
print('\n rfact_zoom = ', rfz)
print(' font_rat = ', fontr, '\n')


print(' i1,i2,j1,j2 =>', i1,i2,j1,j2)

nx_res = i2-i1
ny_res = j2-j1

print(' *** nx_res, ny_res =', nx_res, ny_res)

yx_ratio = float(nx_res+1)/float(ny_res+1)

rnxr = rfz*nx_res ; # widt image (in pixels)
rnyr = rfz*ny_res ; # height image (in pixels)

# Target resolution for figure:
rh_fig = round(rnyr/float(rDPI),3) ; # width of figure
rw_fig = round(rh_fig*yx_ratio      ,3) ; # height of figure
rh_img = rh_fig*float(rDPI)
rw_img = rw_fig*float(rDPI)
while rw_img < round(rnxr,0):
    rw_fig = rw_fig + 0.01/float(rDPI)
    rw_img = rw_fig*float(rDPI)
while rh_img < round(rnyr,0):
    rh_fig = rh_fig + 0.01/float(rDPI)
    rh_img = rh_fig*float(rDPI)
    print(' *** size figure =>', rw_fig, rh_fig, '\n')
    print(' *** Forecasted dimension of image =>', rw_img, rh_img)

print('\n================================================================\n\n\n')

# Field-specific information (colormaps, bounds, etc): lilo
#fa = cp.field_aspect( CWHAT, cbox=CBOX )

#col_cb = fa.color_top_cb

l_3d_field = False







cp.chck4f(cf_mm)
cp.chck4f(cfx_in)
cp.chck4f(cfy_in)

l_notime=False

with Dataset(cfx_in) as id_f:
    list_var = id_f.variables.keys()
    #
    if not cvx_in  in list_var:
        print(' PROBLEM [nemo_movie_vect.py]: variable `'+cvx_in+'` not present in X input file! => exiting!'); exit(0)
        #
    if 'time_instant' in list_var:
        cv_time = 'time_instant'
    elif 'time_counter' in list_var:
        cv_time = 'time_counter'
    elif 'time' in list_var:
        cv_time = 'time'
    else:
        l_notime=True
        print('Did not find a time variable! Assuming no time and Nt=1')
        #if not lForceD0:
        #    print('    ==> then use the `-s` switch to specify an initial date!!!')
        exit(0)
    vtime = id_f.variables[cv_time][:]
    time_units = id_f.variables[cv_time].units
    time_caldr = id_f.variables[cv_time].calendar

with Dataset(cfy_in) as id_f:
    list_var = id_f.variables.keys()
    #
    if not cvy_in  in list_var:
        print(' PROBLEM [nemo_movie_vect.py]: variable `'+cvy_in+'` not present in Y input file! => exiting!'); exit(0)



Nt = 1
if not l_notime: Nt = len(vtime)

if l_show_lsm or l_do_crl or l_do_cof or l_do_cspd or l_do_tke or l_do_eke:
    cv_msk = 'tmask'
    print('\nReading record metrics in '+cf_mm)
    id_lsm = Dataset(cf_mm)
    nb_dim = len(id_lsm.variables[cv_msk].dimensions)
    print(' The mesh_mask has '+str(nb_dim)+' dimmensions!')
    if l_show_lsm:
        if l_do_crl or l_do_cof:  cv_msk = 'fmask'
        print('\n Reading mask as ''+cv_msk+'' in:'); print(cv_msk,'\n')
        if nb_dim==4: XMSK = id_lsm.variables[cv_msk][0,jk,j1:j2,i1:i2]
        if nb_dim==3: XMSK = id_lsm.variables[cv_msk][0,j1:j2,i1:i2]
        if nb_dim==2: XMSK = id_lsm.variables[cv_msk][j1:j2,i1:i2]
        (nj,ni) = np.shape(XMSK)
    if l_do_crl or l_do_cof:
        # e2v, e1u, e1f, e2f
        e2v = id_lsm.variables['e2v'][0,j1:j2,i1:i2]
        e1u = id_lsm.variables['e1u'][0,j1:j2,i1:i2]
        e1f = id_lsm.variables['e1f'][0,j1:j2,i1:i2]
        e2f = id_lsm.variables['e2f'][0,j1:j2,i1:i2] # 
    if l_do_cof:
        ## Coriolis Parameter:
        ff  = id_lsm.variables['gphif'][0,j1:j2,i1:i2]
        ff[:,:] = 2.*romega*np.sin(ff[:,:]*np.pi/180.0)
        (nj,ni) = np.shape(XMSK)

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

    print('Shape Arrays => ni,nj =', ni,nj)

    print('Done!\n')



idx_land = np.where( XMSK <  0.5 )
#idx_ocea = np.where( XMSK >= 0.5 )

#(vjj_land, vji_land) = idx_land
#print(idx_land
#print(vjj_land
#print(vji_land
#print(len(vjj_land), len(vji_land)
#sys.exit(0)


XLSM = np.zeros((nj,ni)) ; # will define real continents not NEMO mask...

if l_add_topo_land:
    cp.chck4f(cf_topo_land)
    id_top = Dataset(cf_topo_land)
    print(' *** Reading "z" into:\n'+cf_topo_land)
    xtopo = id_top.variables['z'][0,j1:j2,i1:i2]
    id_top.close()
    if np.shape(xtopo) != (nj,ni):
        print('ERROR: topo and mask do not agree in shape!'); sys.exit(0)
    xtopo = xtopo*(1. - XMSK)
    #cp.dump_2d_field('topo_'+CBOX+'.nc', xtopo, name='z')    
    if l3d: xtopo = xtopo + rof_dpt
    xtopo[np.where( XMSK > 0.01)] = np.nan
    if nemo_box.l_fill_holes_k and not l3d:
        XLSM[np.where( xtopo < 0.0)] = 1
        xtopo[np.where( xtopo < 0.0)] = np.nan

XLSM[np.where( XMSK > 0.5)] = 1


# Font style:
kk = fsm( fontr, clr_top=color_top, clr_top_cb=color_top_cb )


# Colormaps for fields:
pal_fld = cp.chose_colmap(cpal_fld)
if   l_log_field:
    norm_fld = colors.LogNorm(                   vmin=tmin, vmax=tmax, clip=False)
elif l_pow_field:
    norm_fld = colors.PowerNorm(gamma=pow_field, vmin=tmin, vmax=tmax, clip=False)
else:
    norm_fld = colors.Normalize(                 vmin=tmin, vmax=tmax, clip=False)


if l_show_lsm or l_add_topo_land:
    if l_add_topo_land:
        xtopo = np.log10(xtopo+rof_log)
        pal_lsm = cp.chose_colmap('gray_r')
        #norm_lsm = colors.Normalize(vmin = np.log10(min(-100.+rof_dpt/3.,0.) + rof_log), vmax = np.log10(4000.+rof_dpt + rof_log), clip = False)
        norm_lsm = colors.Normalize(vmin = np.log10(-100. + rof_log), vmax = np.log10(4000.+rof_dpt + rof_log), clip = False)
    else:
        #pal_lsm = cp.chose_colmap('land_dark')
        pal_lsm = cp.chose_colmap('landm')
        norm_lsm = colors.Normalize(vmin=0., vmax=1., clip=False)


#if lForceD0:
#    cyr0=csd0[0:4]
#    cmn0=csd0[5:7]
#    cdd0=csd0[8:10]
#else:
csd0 = num2date(vtime[0], time_units, time_caldr) ; #LOLOsolution
cyr0=str(csd0.year)
cmn0 = '%2.2i'%(int(csd0.month))
cdd0 = '%2.2i'%(int(csd0.day))
print(' ==> csd0, cyr0, cmn0, cdd0 = ',csd0, cyr0, cmn0, cdd0,'\n')



# Time step as a string
if not len(cdt)==2:
    print('ERROR: something is wrong with the format of your time step => '+cdt+' !'); sys.exit(0)
if cdt=='1d':
    dt = 24 ; ntpd = 1
elif cdt[1]=='h':
    dt = int(cdt[0]) ; ntpd = 24 / dt
else:
    print('ERROR: something is wrong with the format of your time step => '+cdt+' !'); sys.exit(0)

vm = vmn
if isleap(int(cyr0)): vm = vml
#print(' year is ', vm, np.sum(vm)

jd = int(cdd0) - 1
jm = int(cmn0)








if l_do_eke:
    Umean = np.zeros((nj,ni))
    Vmean = np.zeros((nj,ni))
    # Need to go for the mean first!!!
    print('\n Goint for '+str(Nt)+' snaphots to compute the mean first!!!')
    for jt in range(jt0,Nt):
        print('Reading record #'+str(jt)+'...')

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
        
        print('Done!')
        # On T-points:
        Umean[:,2:ni] = Umean[:,2:ni] + 0.5*( XFLD[:,1:ni-1] + XFLD[:,2:ni] )/float(Nt)
        Vmean[2:nj,:] = Vmean[2:nj,:] + 0.5*( YFLD[1:nj-1,:] + YFLD[2:nj,:] )/float(Nt)
    print(' time averaging done...\n\n')



    

vm = vmn
if isleap(int(cyr0)): vm = vml
#print(' year is ', vm, np.sum(vm)

jd = int(cdd0) - 1
jm = int(cmn0)








Xplot = np.zeros((nj,ni))
if nemo_box.l_add_quiver:
    #VX = np.arange(0,ni,10)
    #VY = np.arange(0,nj,10)
    VX = np.arange(ni)
    VY = np.arange(nj)
    XU = np.zeros((nj,ni))
    XV = np.zeros((nj,ni))



# Opening netCDF files:
id_fx = Dataset(cfx_in)
id_fy = Dataset(cfy_in)


for jt in range(jt0,Nt):

    # Calendar stuff:
    itime = int( id_f.variables[cv_time][jt] )
    cdt = num2date(itime, time_units, time_caldr)
    cyr, cmo, cdd = str(cdt.year), '%2.2i'%(int(cdt.month)), '%2.2i'%(int(cdt.day))
    chh, cmn, csc = '%2.2i'%(int(cdt.hour)), '%2.2i'%(int(cdt.minute)), '%2.2i'%(int(cdt.second))
    cdate = cyr+cmo+cdd+'_'+chh+cmn            ; # as in output figure name
    cdats = cyr+'-'+cmo+'-'+cdd+' '+chh+':'+cmn #+':'+csc  ; # to display on figure

    # Name of figure to generate:
    cstr = CBOX
    if CEXP: cstr += '_'+CEXP
    if cn_fig != "":
        cfig = cdir_figs+'/'+cv_out+'_'+cstr+'_'+cn_fig+'_'+cdate+'.'+fig_type
    else:
        if l3d:
            cfig = cdir_figs+'/'+cv_out+'_'+cstr+'_'+CNEMO+CRUN+'_lev'+str(jk)+'_'+CBOX+'_'+cdate+'.'+fig_type
        else:
            cfig = cdir_figs+'/'+cv_out+'_'+cstr+'_'+CNEMO+CRUN+'_'+CBOX+'_'+cdate+'.'+fig_type


    if not path.exists(cfig):
    ###### FIGURE ##############

        fig = plt.figure(num = 1, figsize=(rw_fig, rh_fig), dpi=None, facecolor='w', edgecolor='0.5')
    
        ax  = plt.axes([0., 0., 1., 1.], facecolor = '0.85') # missing seas will be in 'facecolor' !
    
        vc_fld = np.arange(tmin, tmax + df, df)
    
    
        print('Reading record #'+str(jt)+' of '+cvx_in+' in '+cfx_in)
        if l3d: print('            => at level #'+str(jk)+' ('+cdepth+')!')

        if not l_3d_field:
            XFLD  = id_fx.variables[cvx_in][jt,j1:j2,i1:i2] ; # t, y, x
        else:
            print('j1:j2 =', j1,j2)
            print('i1:i2 =', i1,i2)
            XFLD  = id_fx.variables[cvx_in][jt,jk,j1:j2,i1:i2] ; # t, y, x
        print('Done!')
    
        print('Reading record #'+str(jt)+' of '+cvy_in+' in '+cfy_in)
        if l3d: print('            => at level #'+str(jk)+' ('+cdepth+')!')
        if not l_3d_field:
            YFLD  = id_fy.variables[cvy_in][jt,j1:j2,i1:i2] ; # t, y, x
        else:
            YFLD  = id_fy.variables[cvy_in][jt,jk,j1:j2,i1:i2] ; # t, y, x
        print('Done!')
    
    
        lx = np.zeros((nj,ni))
        ly = np.zeros((nj,ni))
    
        if l_do_crl or l_do_cof:
            print('\nComputing curl...')
            lx[:,1:ni-1] =   e2v[:,2:ni]*YFLD[:,2:ni] - e2v[:,1:ni-1]*YFLD[:,1:ni-1]
            ly[1:nj-1,:] = - e1u[2:nj,:]*XFLD[2:nj,:] + e1u[1:nj-1,:]*XFLD[1:nj-1,:]
            if l_do_cof: Xplot[:,:] = ( lx[:,:] + ly[:,:] )*XMSK[:,:] / ( e1f[:,:]*e2f[:,:]*ff[:,:] ) # Relative Vorticity...
            if l_do_crl: Xplot[:,:] = ( lx[:,:] + ly[:,:] )*XMSK[:,:] / ( e1f[:,:]*e2f[:,:] ) * 1000. # Curl...
    
        if l_do_cspd:
            print('\nComputing current speed at T-points ...')
            lx[:,2:ni] = 0.5*( XFLD[:,1:ni-1] + XFLD[:,2:ni] )
            ly[2:nj,:] = 0.5*( YFLD[1:nj-1,:] + YFLD[2:nj,:] )
            Xplot[:,:] = np.sqrt( lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] ) * XMSK[:,:]
            if nemo_box.l_add_quiver:
                XU[:,:] = XFLD[:,:]
                XV[:,:] = YFLD[:,:]
    
        if l_do_tke:
            print('\nComputing TKE at T-points ...')
            lx[:,2:ni] = 0.5*( XFLD[:,1:ni-1] + XFLD[:,2:ni] )
            ly[2:nj,:] = 0.5*( YFLD[1:nj-1,:] + YFLD[2:nj,:] )
            Xplot[:,:] = ( lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] ) * XMSK[:,:]
    
    
        if l_do_eke:
            print('\nComputing TKE at T-points ...')
            lx[:,2:ni] = 0.5*( XFLD[:,1:ni-1] + XFLD[:,2:ni] ) - Umean[:,2:ni]
            ly[2:nj,:] = 0.5*( YFLD[1:nj-1,:] + YFLD[2:nj,:] ) - Vmean[2:nj,:]
            Xplot[:,:] = ( lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] ) * XMSK[:,:]
    
        del lx, ly
        print('... '+cv_out+' computed!\n')
            
        del XFLD,YFLD
    
        print('')
        if not l_show_lsm and jt == jt0: ( nj , ni ) = np.shape(Xplot)
        print('  *** dimension of array => ', ni, nj, np.shape(Xplot))
        

        print('Ploting')
    
        plt.axis([ 0, ni, 0, nj])
        
        if nemo_box.c_imshow_interp == 'none':
            Xplot[idx_land] = np.nan
        else:
            print(" *** drowning array `Xplot`....\n")
            cp.drown(Xplot, XMSK, k_ew=-1, nb_max_inc=10, nb_smooth=10)
    
        if l_save_nc:
            if l3d:
                cf_out = 'nc/'+CWHAT+'_NEMO_'+CNEMO+CRUN+'_lev'+str(jk)+'_'+CBOX+'_'+cdate+'.nc'
            else:
                cf_out = 'nc/'+CWHAT+'_NEMO_'+CNEMO+CRUN+'_'+CBOX+'_'+cdate+'.nc'
            print(' Saving in '+cf_out)
            cp.dump_2d_field(cf_out, Xplot, xlon=Xlon, xlat=Xlat, name=CWHAT)
            print('')
    

        cf = plt.imshow( Xplot[:,:], cmap=pal_fld, norm=norm_fld, interpolation=nemo_box.c_imshow_interp )
    
        if nemo_box.l_add_quiver:
            idx_high = np.where(Xplot[:,:]>tmax)
            XU[idx_high] = np.nan ; XV[idx_high] = np.nan
            idx_miss = np.where(XLSM[:,:]<0.5)
            XU[idx_miss] = np.nan ; XV[idx_miss] = np.nan
            #XU = XU*XLSM ; XV = XV*XLSM 
            nss=nemo_box.n_subsamp_qvr
            cq = plt.quiver( VX[:ni:nss], VY[:nj:nss], XU[:nj:nss,:ni:nss], XV[:nj:nss,:ni:nss], scale=50, color='w', width=0.001, linewidth=0.1 )
    
        #LOLO: rm ???
        if l_show_lsm or l_add_topo_land:
            if l_add_topo_land:
                clsm = plt.imshow(np.ma.masked_where(XLSM>0.0001, xtopo), cmap = pal_lsm, norm = norm_lsm, interpolation='none')
                if nemo_box.c_imshow_interp == 'none':
                    plt.contour(XLSM, [0.9], colors='k', linewidths=0.5)
            else:
                clsm = plt.imshow(np.ma.masked_where(XLSM>0.0001, XLSM), cmap = pal_lsm, norm = norm_lsm, interpolation='none')
    
        ##### COLORBAR ######
        if nemo_box.l_show_cb:
            ax2 = plt.axes(nemo_box.vcb)
            clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend=cb_extend)
            cb_labs = []
            #if cb_jump > 1:
            cpt = 0
            for rr in vc_fld:
                if cpt % cb_jump == 0:
                    if df >= 1.: cb_labs.append(str(int(rr)))
                    if df <  1.: cb_labs.append(str(round(rr,int(np.ceil(np.log10(1./df)))+1) ))
                else:
                    cb_labs.append(' ')
                cpt = cpt + 1
            #else:
            #    for rr in vc_fld: cb_labs.append(str(round(rr,int(np.ceil(np.log10(1./df)))+1) ))
    
            clb.ax.set_xticklabels(cb_labs, **fsm.cfont_clb_tcks)
            clb.set_label(cunit, **fsm.cfont_clb)
            clb.ax.yaxis.set_tick_params(color=color_top_cb) ; # set colorbar tick color
            clb.outline.set_edgecolor(color_top_cb) ; # set colorbar edgecolor
            clb.ax.tick_params(which = 'minor', length = 2, color = color_top_cb )
            clb.ax.tick_params(which = 'major', length = 4, color = color_top_cb )
    
    
        if nemo_box.l_show_clock:
            xl = float(x_clock)/rfz
            yl = float(y_clock)/rfz
            ax.annotate(cdats, xy=(1, 4), xytext=(xl,yl), **fsm.cfont_clock)
    
        if nemo_box.l_show_exp:
            xl = float(x_exp)/rfz
            yl = float(y_exp)/rfz
            ax.annotate('Experiment: '+CNEMO+CRUN, xy=(1, 4), xytext=(xl,yl), **fsm.cfont_exp)

        if l_add_sign and nemo_box.l_show_sign:
            xl = float(x_sign)/rfz
            yl = float(y_sign)/rfz
            ax.annotate(CSIGN, xy=(1, 4), xytext=(xl,yl), zorder=150, **fsm.cfont_sign)

        if nemo_box.l_show_name:
            if nemo_box.loc_name=='land':
                fsm.cfont_titl['color'] = 'w' ; # update font color in font dictionary            
            cbla = CNEMO
            if CONAME != "": cbla = CONAME
            xl = float(x_name)/rfz
            yl = float(y_name)/rfz
            ax.annotate(cbla, xy=(1, 4), xytext=(xl, yl), zorder=52,  **fsm.cfont_titl)

        if nemo_box.l_add_logo:
            datafile = cbook.get_sample_data(dir_logos+'/'+nemo_box.cf_logo_on, asfileobj=False)
            im = image.imread(datafile)
            #im[:, :, -1] = 0.5  # set the alpha channel
            fig.figimage(im, x_logo, y_logo, zorder=9)
            del datafile, im
            #
            if nemo_box.l_add_logo_ige:
                datafile = cbook.get_sample_data(dir_logos+'/'+nemo_box.cf_logo_ige, asfileobj=False)
                im = image.imread(datafile)
                fig.figimage(im, x_logo+144, y_logo-150., zorder=9)
                del datafile, im
                #
            if nemo_box.l_add_logo_prc:
                datafile = cbook.get_sample_data(dir_logos+'/'+nemo_box.cf_logo_prc, asfileobj=False)
                im = image.imread(datafile)
                fig.figimage(im, x_logo-77, y_logo-140., zorder=9)
                del datafile, im
    
        plt.savefig(cfig, dpi=rDPI, orientation='portrait', facecolor='k')
        print(cfig+' created!\n')
        plt.close(1)

        if l_show_lsm: del clsm
        del cf, fig, ax
        if nemo_box.l_show_cb: del clb

    else:
        print('\n Figure '+cfig+' already there!\n')


id_fx.close()
id_fy.close()


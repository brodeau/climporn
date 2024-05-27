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

from netCDF4 import Dataset

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


vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

fig_type='png'
rDPI = 240

jt0 = 0

jk=0
l_mask_no_ice = False


rof_log = 150.
rof_dpt = 0.

grav = 9.80665    # acceleration of gravity, [m/s^2] (same as in NEMO 3.6)

l_save_nc = False ; # save the field we built in a netcdf file !!!

romega = 7.292115083046062E-005 # Coriolis [1/s] (same as in NEMO 3.6 / #romega = 2.*np.pi/86400.0)


# Normally logos should be found there:
dir_logos = str.replace( path.dirname(path.realpath(__file__)) , 'climporn/python', 'climporn/misc/logos' )
print("\n --- logos found into : "+dir_logos+" !\n")



def moduleOfVelocity( cvar_u, cvar_v, pSSU, pSSV ):
    #
    print(' *** Computing module of surface velocity "'+cvar_u+','+cvar_v+'" !')
    #
    (nj,ni) = np.shape( pSSU )
    #
    zu, zv     = np.zeros((nj,ni)), np.zeros((nj,ni))
    #
    zu[:,1:] = 0.5 * ( pSSU[:,1:] + pSSU[:,:ni-1] ) ; # U @T
    zv[1:,:] = 0.5 * ( pSSV[1:,:] + pSSV[:nj-1,:] ) ; # U @T
    #
    return np.sqrt( zu*zu + zv*zv )



def comp_LapOfSSH( cvar, pe1t, pe2t, pe1u, pe2u, pe1v, pe2v, pSSH ):
    #
    print(' *** Computing curvilinear Laplacian of SSH "'+cvar+'" !')
    (nj,ni) = np.shape( pSSH )
    #
    zL         = np.zeros((nj,ni))
    zx, zy     = np.zeros((nj,ni)), np.zeros((nj,ni))
    dFdx, dFdy = np.zeros((nj,ni)), np.zeros((nj,ni))
    #
    dFdx[:,:ni-1] = 1.E9*( pSSH[:,1:ni] - pSSH[:,:ni-1] ) / pe1u[:,:ni-1]  ; # => at U-point `i`
    dFdy[:nj-1,:] = 1.E9*( pSSH[1:nj,:] - pSSH[:nj-1,:] ) / pe2v[:nj-1,:]  ; # => at V-point `j`
    #
    zx[:,:] = 1.E6*pe2u[:,:]/pe1u[:,:]*dFdx[:,:]
    zy[:,:] = 1.E6*pe1v[:,:]/pe2v[:,:]*dFdy[:,:]
    #
    zL[1:nj,1:ni] = (  (zx[1:nj,1:ni]-zx[1:nj,:ni-1])/pe1t[1:nj,1:ni] \
                     + (zy[1:nj,1:ni]-zy[:nj-1,1:ni])/pe2t[1:nj,1:ni] ) / ( pe1t[1:nj,1:ni]*pe2t[1:nj,1:ni] )
    #
    return zL
















################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--fin' , required=True,                help='NEMO netCDF file to read from...')
requiredNamed.add_argument('-w', '--what', required=True, help='field/diagnostic to plot (ex: CSPEED,CURLOF,ect.)')

parser.add_argument('-C', '--conf', default="none",           help='name of NEMO config (ex: eNATL60) (defined into `nemo_hboxes.py`)')
parser.add_argument('-E', '--nexp', default=None,             help='name of experiment (shows up in figure name)')
parser.add_argument('-b', '--box' , default="ALL",            help='extraction box name (ex: ALL) (defined into `nemo_hboxes.py`)')
parser.add_argument('-m', '--fmm' , default="mesh_mask.nc",   help='NEMO mesh_mask file (ex: mesh_mask.nc)')
parser.add_argument('-s', '--sd0' , default=None,             help='initial date as <YYYYMMDD>')
parser.add_argument('-l', '--lev' , type=int, default=0,      help='level to use if 3D field (default: 0 => 2D)')
parser.add_argument('-z', '--zld' ,                           help='topography netCDF file to use (field="z")')
parser.add_argument('-t', '--tstep',  default="1h",           help='time step ("1h","2h",..,up to "1d") in input file')
parser.add_argument('-N', '--oname',  default="",             help='a name that overides `CONF` on the plot...')
parser.add_argument('-o', '--outdir', default="./figs",       help='path to directory where to save figures')
parser.add_argument('-f', '--fignm',  default="",             help='common string in name of figure to create')
parser.add_argument('-T', '--addSST', default="",             help='show the SST over open ocean!')
parser.add_argument('-H', '--addSSH', default="",             help='show the SSH over open ocean!')
parser.add_argument('-U', '--addSSU', default="",             help='show the module of ocean surface velocity over open ocean!')
parser.add_argument('-V', '--addVOR', default="",             help='show the vorticity (laplacian of SSH) over open ocean!')
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
cf_in = args.fin
cf_mm = args.fmm
csd0  = args.sd0
jk    = args.lev
cf_topo_land = args.zld
cdt    = args.tstep  ; # time step, in the form "1h", "2h", ..., "12h", ..., "1m", ..., "6m", ... "1y", ..., etc
CONAME = args.oname
cd_out = args.outdir
cn_fig = args.fignm
caSST  = args.addSST
caSSH  = args.addSSH
caSSU  = args.addSSU
caVOR  = args.addVOR
CSIGN  = args.sign
lcb    = args.cb
cltr   = args.ltr

print('')
print(' *** CNEMO = ', CNEMO)
print(' *** CBOX  = ', CBOX)
print(' *** CWHAT = ', CWHAT)
print(' *** cf_in = ', cf_in)
print(' *** cf_mm = ', cf_mm)
print(' *** show colorbar:', lcb)


lForceD0 = False
if csd0:
    lForceD0 = True
    print(' *** csd0 = ', csd0)

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

l_add_SST_to_ice_field = False
if caSST != "":
    l_add_SST_to_ice_field = True
    print(' *** We shall add following SST field below ice: ', caSST)

l_add_SSH_to_ice_field = False
if caSSH != "":
    l_add_SSH_to_ice_field = True
    print(' *** We shall add following SSH field below ice: ', caSSH)

l_add_SSU_to_ice_field = False
if caSSU != "":
    l_add_SSU_to_ice_field = True
    caSSV = str.replace(caSSU, 'u', 'v')
    print(' *** We shall add following SSU field below ice: ', caSSU,caSSV)

l_add_VOR_to_ice_field = False
if caVOR != "":
    l_add_VOR_to_ice_field = True
    print(' *** We shall show the vorticity over open ocean based on Laplacian of : ', caVOR)



i_add_something = int(l_add_SST_to_ice_field) + int(l_add_SSH_to_ice_field) + int(l_add_SSU_to_ice_field) + int(l_add_VOR_to_ice_field)
    
if i_add_something > 1:
    #if l_add_SST_to_ice_field and l_add_VOR_to_ice_field:
    print(' ERROR: chose between the `-T`, `-H`, & `-V` options, only 1 !!!'); exit(0)


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

#print(CNEMO,CBOX)
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
#if cltr:                  (x_ltr,y_ltr)     = nemo_box.ltr


#---------------------------------------------------------------




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
fa = cp.field_aspect( CWHAT, cbox=CBOX )



# SST over open ocean:
if l_add_SST_to_ice_field:
    cpal_sst = 'YlGnBu_r'
    rmin_sst = -2. ; rmax_sst = 14. ; dsst = 2.
    pal_sst = cp.chose_colmap(cpal_sst)
    norm_sst = colors.Normalize(vmin=rmin_sst, vmax=rmax_sst , clip = False)

# SSH over open ocean:
if l_add_SSH_to_ice_field:
    #cpal_ssh = 'cmocean_deep'
    cpal_ssh = 'cmocean_matter'
    rmin_ssh = -1.4 ; rmax_ssh = 0.9 ; dssh = 0.1
    pal_ssh = cp.chose_colmap(cpal_ssh)
    norm_ssh = colors.Normalize(vmin=rmin_ssh, vmax=rmax_ssh , clip = False)

# SSU over open ocean:
if l_add_SSU_to_ice_field:
    #cpal_ssu = 'cmocean_amp'
    #cpal_ssu = 'cmocean_deep_r'
    cpal_ssu = 'cmocean_tempo'
    #cpal_ssu = 'cmocean_ice_r'
    if CWHAT in ['sivolu','siconc']: cpal_ssu = 'cmocean_thermal'
    rmin_ssu = 0. ; rmax_ssu = 1. ; dssu = 0.1
    pal_ssu = cp.chose_colmap(cpal_ssu)
    norm_ssu = colors.Normalize(vmin=rmin_ssu, vmax=rmax_ssu , clip = False)

# VOR over open ocean:
if l_add_VOR_to_ice_field:
    #cpal_vor = 'cmocean_curl'
    cpal_vor = 'PiYG_r'
    #cpal_vor = 'BrBG_r'
    rmin_vor = -2500. ; rmax_vor = -rmin_vor ; dvor = 250. ; # NANUK36
    pal_vor = cp.chose_colmap(cpal_vor)
    norm_vor = colors.Normalize(vmin=rmin_vor, vmax=rmax_vor , clip = False)

    

# Ice over ocean field:
if fa.l_show_ice:
    cf_ice = str.replace(cf_in, 'gridT-2D', 'icemod')
    rmin_ice = 0.25
    #cpal_ice = 'ncview_bw'
    #cpal_ice = 'Blues_r'
    cpal_ice = 'bone'
    vcont_ice = np.arange(rmin_ice, 1.05, 0.05)
    #
    pal_ice = cp.chose_colmap(cpal_ice)
    norm_ice = colors.Normalize(vmin = rmin_ice, vmax = 0.95, clip = False)


# Ice
l_i_need_A = ( fa.l_show_ice or fa.imask_no_ice_pc>0 or i_add_something>0 )


#print('LOLO: cv_in =',fa.cv_in,' fa.imask_no_ice_pc =',fa.imask_no_ice_pc)
#exit(0)

rt_io = 0.0
#if fa.imask_no_ice_pc > 0:
if l_i_need_A:
    # => going to read variable `fa.nm_ice_conc` into the same file
    cf_ice = cf_in
    rt_io = float(fa.imask_no_ice_pc)/100.





if fa.l_show_ice: cp.chck4f(cf_ice)

cp.chck4f(cf_mm)
cp.chck4f(cf_in)

l_notime=False

id_f = Dataset(cf_in)
list_var = id_f.variables.keys()
if 'time_instant' in list_var:
    cv_time = 'time_instant'
elif 'time_counter' in list_var:
    cv_time = 'time_counter'
elif 'time' in list_var:
    cv_time = 'time'
else:
    l_notime=True
    print('Did not find a time variable! Assuming no time and Nt=1')
    if not lForceD0:
        print('    ==> then use the `-s` switch to specify an initial date!!!')
        exit(0)
vtime = id_f.variables[cv_time][:]
id_f.close()

Nt = 1
if not l_notime: Nt = len(vtime)

if fa.l_show_lsm:
    print('\nReading record metrics in '+cf_mm)
    id_lsm = Dataset(cf_mm)
    nb_dim = len(id_lsm.variables[fa.cv_msk].dimensions)
    print(' The mesh_mask has '+str(nb_dim)+' dimmensions!')
    print(' *** Reading '+fa.cv_msk+' !')
    if nb_dim==4: XMSK = id_lsm.variables[fa.cv_msk][0,jk,j1:j2,i1:i2]
    if nb_dim==3: XMSK = id_lsm.variables[fa.cv_msk][jk,j1:j2,i1:i2]
    if nb_dim==2: XMSK = id_lsm.variables[fa.cv_msk][j1:j2,i1:i2]
    (nj,ni) = np.shape(XMSK)

    if fa.l_apply_lap or l_add_VOR_to_ice_field:
        print(' *** Reading e1t and e2t !')
        XE1T  = id_lsm.variables['e1t'][0,j1:j2,i1:i2]
        XE2T  = id_lsm.variables['e2t'][0,j1:j2,i1:i2]
        XE1U  = id_lsm.variables['e1u'][0,j1:j2,i1:i2]
        XE2U  = id_lsm.variables['e2u'][0,j1:j2,i1:i2]
        XE1V  = id_lsm.variables['e1v'][0,j1:j2,i1:i2]
        XE2V  = id_lsm.variables['e2v'][0,j1:j2,i1:i2]
        #XE1T2 = XE1T*XE1T
        #XE2T2 = XE2T*XE2T
    if fa.l_apply_hgrad or fa.l_apply_geov:
        print(' *** Reading e1u and e2v !')
        e1u = id_lsm.variables['e1u'][0,j1:j2,i1:i2]
        e2v = id_lsm.variables['e2v'][0,j1:j2,i1:i2]
        print(' *** Reading umask and vmask !')
        if nb_dim==4:
            UMSK = id_lsm.variables['umask'][0,jk,j1:j2,i1:i2]
            VMSK = id_lsm.variables['vmask'][0,jk,j1:j2,i1:i2]
        if nb_dim==3:
            UMSK = id_lsm.variables['umask'][jk,j1:j2,i1:i2]
            VMSK = id_lsm.variables['vmask'][jk,j1:j2,i1:i2]
        if nb_dim==2:
            UMSK = id_lsm.variables['umask'][j1:j2,i1:i2]
            VMSK = id_lsm.variables['vmask'][j1:j2,i1:i2]

    if fa.l_apply_geov or fa.l_apply_lap or l_add_VOR_to_ice_field:
        ## Coriolis Parameter:
        ff  = id_lsm.variables['gphif'][0,j1:j2,i1:i2]
        ff[:,:] = 2.*romega*np.sin(ff[:,:]*np.pi/180.0)

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
    if l3d: xtopo = xtopo + rof_dpt
    xtopo[np.where( xtopo <= 0 )] = 0
    xtopo = np.log10(xtopo+rof_log)
    xtopo[np.where(XMSK > 0.01)] = np.nan
    #
    #if nemo_box.l_fill_holes_k and not l3d:
    #    XLSM[np.where( xtopo < 0.0)] = 1
    #    #xtopo[np.where( xtopo < 0.0)] = np.nan
    #cp.dump_2d_field('topo_'+CBOX+'.nc', xtopo, name='z')

XLSM[np.where( XMSK > 0.5)] = 1


# Font style:
kk = fsm( fontr, clr_top=fa.color_top, clr_top_cb=fa.color_top_cb )


# Colormaps for fields:
pal_fld = cp.chose_colmap(fa.cpal_fld)
if   fa.l_log_field:
    norm_fld = colors.LogNorm(                   vmin=fa.tmin, vmax=fa.tmax, clip=False)
elif fa.l_pow_field:
    norm_fld = colors.PowerNorm(gamma=fa.pow_field, vmin=fa.tmin, vmax=fa.tmax, clip=False)
else:
    norm_fld = colors.Normalize(                 vmin=fa.tmin, vmax=fa.tmax, clip=False)


if fa.l_show_lsm or l_add_topo_land:
    if l_add_topo_land:
        pal_lsm = cp.chose_colmap('gray_r')
        #norm_lsm = colors.Normalize(vmin = np.log10(min(-100.+rof_dpt/3.,0.) + rof_log), vmax = np.log10(4000.+rof_dpt + rof_log), clip = False)
        norm_lsm = colors.Normalize(vmin = np.log10(-100. + rof_log), vmax = np.log10(4000.+rof_dpt + rof_log), clip = False)
    else:
        #pal_lsm = cp.chose_colmap('land_dark')
        pal_lsm = cp.chose_colmap('landm')
        norm_lsm = colors.Normalize(vmin=0., vmax=1., clip=False)

if not lForceD0:
    csd0 = cp.epoch2clock( int(vtime[jt0]), frmt='nodash' )
#
cyr0=csd0[0:4]
cmn0=csd0[4:6]
cdd0=csd0[6:8]

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


# Opening netCDF files:
id_f = Dataset(cf_in)
id_f.set_auto_mask(False) ;# prevent the application of `valid_min` and `valid_max`....

for jt in range(jt0,Nt):

    itime = int( id_f.variables[cv_time][jt] )

    cdate = cp.epoch2clock(itime, precision='h')
    cdats = cp.epoch2clock(itime)



    # Name of figure to generate:
    cstr = CBOX
    if CEXP: cstr += '_'+CEXP
    if cn_fig != "":
        cfig = cdir_figs+'/'+fa.cv_out+'_'+cstr+'_'+cn_fig+'_'+cdate+'.'+fig_type
    else:
        if l3d:
            cfig = cdir_figs+'/'+fa.cv_out+'_'+cstr+'_'+CNEMO+CRUN+'_lev'+str(jk)+'_'+CBOX+'_'+cdate+'.'+fig_type
        else:
            cfig = cdir_figs+'/'+fa.cv_out+'_'+cstr+'_'+CNEMO+CRUN+'_'+CBOX+'_'+cdate+'.'+fig_type


    if not path.exists(cfig):
    ###### FIGURE ##############

        fig = plt.figure(num = 1, figsize=(rw_fig, rh_fig), dpi=rDPI, facecolor='w', edgecolor='0.5')

        ax  = plt.axes([0., 0., 1., 1.], facecolor=fa.color_missing) # missing seas will be in 'facecolor' !

        if len(fa.vc_fld_force)>0:
            vc_fld = fa.vc_fld_force
        else:
            vc_fld = np.arange(fa.tmin, fa.tmax + fa.df, fa.df)



        print('Reading record #'+str(jt)+' of '+fa.cv_in+' in '+cf_in)
        if l3d: print('            => at level #'+str(jk)+' ('+cdepth+')!')

        if l_notime:
            if l3d:
                Xplot  = id_f.variables[fa.cv_in][jk,j1:j2,i1:i2]
            else:
                Xplot  = id_f.variables[fa.cv_in][j1:j2,i1:i2]
        else:
            if l3d:
                Xplot  = id_f.variables[fa.cv_in][jt,jk,j1:j2,i1:i2] ; # t, y, x
            else:
                Xplot  = id_f.variables[fa.cv_in][jt,j1:j2,i1:i2] ; # t, y, x


        if l_i_need_A:
            print('Reading record #'+str(jt)+' of '+fa.nm_ice_conc+' in '+cf_ice)
            id_ice = Dataset(cf_ice)
            XICE  = id_ice.variables[fa.nm_ice_conc][jt,j1:j2,i1:i2] ; # t, y, x
            id_ice.close()
            print('Done!\n')

        if l_add_SST_to_ice_field:
            Xpsst = id_f.variables[caSST][jt,j1:j2,i1:i2] ; # t, y, x

        if l_add_SSH_to_ice_field:
            Xpssh = id_f.variables[caSSH][jt,j1:j2,i1:i2] ; # t, y, x

        if l_add_SSU_to_ice_field:
            Xpssu = id_f.variables[caSSU][jt,j1:j2,i1:i2] ; # t, y, x
            Xpssv = id_f.variables[caSSV][jt,j1:j2,i1:i2] ; # t, y, x
            Xpssu = moduleOfVelocity( caSSU, caSSV, Xpssu, Xpssv )
            del Xpssv

        if l_add_VOR_to_ice_field:
            Xpssh = id_f.variables[caVOR][jt,j1:j2,i1:i2] ; # t, y, x
            Xpvor = comp_LapOfSSH( caVOR, XE1T, XE2T, XE1U, XE2U, XE1V, XE2V, Xpssh )
            Xpvor[:,:] = Xpvor[:,:] / ff[:,:]
            del Xpssh

        print('Done!\n')

        if fa.rmult != 1.: Xplot[:,:] = fa.rmult * Xplot[:,:]


        if fa.l_apply_lap:
            Xplot = comp_LapOfSSH( fa.cv_in, XE1T, XE2T, XE1U, XE2U, XE1V, XE2V, Xplot )

        if fa.l_apply_hgrad:
            print(' *** Computing gradient of "'+fa.cv_in+'"!')
            lx = np.zeros((nj,ni))
            ly = np.zeros((nj,ni))

            if fa.l_smooth: cp.smoother(Xplot, XMSK, nb_smooth=fa.nb_smooth)

            # Zonal gradient on T-points:
            lx[:,1:ni-1] = (Xplot[:,2:ni] - Xplot[:,0:ni-2]) / (e1u[:,1:ni-1] + e1u[:,0:ni-2]) * UMSK[:,1:ni-1] * UMSK[:,0:ni-2]
            lx[:,:] = XMSK[:,:]*lx[:,:]

            # Meridional gradient on T-points:
            ly[1:nj-1,:] = (Xplot[2:nj,:] - Xplot[0:nj-2,:]) / (e2v[1:nj-1,:] + e2v[0:nj-2,:]) * VMSK[1:nj-1,:] * VMSK[0:nj-2,:]
            ly[:,:] = XMSK[:,:]*ly[:,:]

            Xplot[:,:] = 0.0
            # Modulus of vector gradient:
            Xplot[:,:] = np.sqrt(  lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] )

            del lx, ly


        if fa.l_apply_geov:
            print(' *** Computing gradient of "'+fa.cv_in+'"!')
            lx = np.zeros((nj,ni))
            ly = np.zeros((nj,ni))

            # Zonal gradient on T-points:
            lx[:,1:ni-1] = (Xplot[:,2:ni] - Xplot[:,0:ni-2]) / (e1u[:,1:ni-1] + e1u[:,0:ni-2]) * UMSK[:,1:ni-1] * UMSK[:,0:ni-2]
            lx[:,:] = XMSK[:,:]*lx[:,:]

            # Meridional gradient on T-points:
            ly[1:nj-1,:] = (Xplot[2:nj,:] - Xplot[0:nj-2,:]) / (e2v[1:nj-1,:] + e2v[0:nj-2,:]) * VMSK[1:nj-1,:] * VMSK[0:nj-2,:]
            ly[:,:] = XMSK[:,:]*ly[:,:]

            Xplot[:,:] = 0.0
            # Modulus of vector gradient:
            Xplot[:,:] = grav/ff * np.sqrt( lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] )

            del lx, ly

        if fa.lSubstractMean:
            zmean = np.sum( Xplot[:,:]*XMSK[:,:])/np.sum(XMSK[:,:])
            Xplot[:,:] = Xplot[:,:] - zmean



        print('')
        if not fa.l_show_lsm and jt == jt0: ( nj , ni ) = np.shape(Xplot)
        print('  *** dimension of array => ', ni, nj, np.shape(Xplot))


        print('Ploting')

        plt.axis([ 0, ni, 0, nj])


        cinterp = nemo_box.c_imshow_interp ; # should be 'none' normally...
        if rfz<1.:
            cinterp = 'antialiased'
        #cintdat = 'rgba' ; # data/rgba
        cintdat = 'data' ; # data/rgba

        if not l_add_topo_land:            
            if cinterp == 'none':
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


        # if not i_add_something !!
        #lili!
        #if not l_add_topo_land and l_i_need_A and rfz>=1.:
        #    # masking field where too litle sea-ice:
        #    Xplot = np.ma.masked_where( XICE <= rt_io, Xplot )
                    
        #cf = plt.figimage( np.flipud(Xplot[:,:]), cmap=pal_fld, norm=norm_fld,zorder=10 )
        #cf = plt.imshow( Xplot[:,:], cmap=pal_fld, norm=norm_fld, interpolation=cinterp, zorder=10 )
        if rfz<1.:
            print(' *** Using `imshow()` with `interpolation='+cinterp+'` & `interpolation_stage='+cintdat+'`')
            # Field is subsampled, better to use a dealiasing => ... interpolation=interp, interpolation_stage=data/rgba,
            cf = plt.imshow(     Xplot[:,:], cmap=pal_fld, norm=norm_fld, interpolation=cinterp, interpolation_stage=cintdat, zorder=10 )
        else:
            print(' *** Using `pcolormesh()` with `interpolation='+cinterp+'` & `interpolation_stage='+cintdat+'`')
            cf = plt.pcolormesh( Xplot[:,:], cmap=pal_fld, norm=norm_fld, zorder=10 )

            
        # Add SST over open ocean:
        if l_add_SST_to_ice_field:
            psst = np.ma.masked_where(XICE > rt_io, Xpsst)
            ct   = plt.imshow(psst, cmap=pal_sst, norm=norm_sst, interpolation=cinterp, interpolation_stage=cintdat, zorder=12)
            del psst, ct

        # Add SSH over open ocean:
        if l_add_SSH_to_ice_field:
            pssh = np.ma.masked_where(XICE > rt_io, Xpssh)
            ct   = plt.imshow(pssh, cmap=pal_ssh, norm=norm_ssh, interpolation=cinterp, interpolation_stage=cintdat, zorder=12)
            del pssh, ct

        # Add SSU over open ocean:
        if l_add_SSU_to_ice_field:
            pssu = np.ma.masked_where(XICE > rt_io, Xpssu)
            ct   = plt.imshow(pssu, cmap=pal_ssu, norm=norm_ssu, interpolation=cinterp, interpolation_stage=cintdat, zorder=12)
            del pssu, ct

        # Add vorticity aka Laplacian of SSH over open ocean:
        if l_add_VOR_to_ice_field:
            pvor = np.ma.masked_where(XICE > rt_io, Xpvor)
            ct   = plt.imshow(pvor, cmap=pal_vor, norm=norm_vor, interpolation=cinterp, interpolation_stage=cintdat, zorder=12)
            del pvor, ct

        # Add Sea-Ice onto a open ocean field:
        if fa.l_show_ice:
            pice = np.ma.masked_where(XICE < rmin_ice, XICE)
            ci = plt.imshow(pice, cmap=pal_ice, norm=norm_ice, interpolation=cinterp, interpolation_stage=cintdat, zorder=12) ; del pice, ci

        if fa.l_show_lsm or l_add_topo_land:
            if l_add_topo_land:
                if rfz<1.:
                    clsm = plt.imshow(     np.ma.masked_where(XLSM>0.0001, xtopo), cmap=pal_lsm, norm=norm_lsm, interpolation=cinterp, interpolation_stage=cintdat,zorder=50 )
                else:
                    clsm = plt.pcolormesh( np.ma.masked_where(XLSM>0.0001, xtopo), cmap=pal_lsm, norm=norm_lsm, zorder=50 )
                #if cinterp == 'none':
                plt.contour(XLSM, [0.5], colors='k', linewidths=rfz*0.6, zorder=100)
            else:
                pmsk = np.ma.masked_where(XLSM[:,:] > 0.2, XLSM[:,:])
                clsm = plt.imshow( pmsk, cmap=pal_lsm, norm=norm_lsm, interpolation=cinterp, interpolation_stage=cintdat )
                del pmsk

        ##### COLORBAR ######
        if nemo_box.l_show_cb and lcb:
            ax2 = plt.axes(nemo_box.vcb)
            if fa.l_pow_field or fa.l_log_field:
                clb = mpl.colorbar.ColorbarBase(ax=ax2, ticks=fa.vc_fld_powlog, cmap=pal_fld, norm=norm_fld,
                                                orientation='horizontal', extend=fa.cb_extend)
                cb_labs = np.array( [ str(i) for i in fa.vc_fld_powlog ], dtype='U16' )
                cb_labs[(cb_labs=='0.0')]  = '0'
                cb_labs[(cb_labs=='1.0')]  = '1'
            else:
                clb = mpl.colorbar.ColorbarBase(ax=ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld,
                                                orientation='horizontal', extend=fa.cb_extend)
                if len(fa.vc_fld_force)>0:
                    cb_labs = np.array( [ str(i) for i in fa.vc_fld_force ], dtype='U16' )
                    cb_labs[(cb_labs=='0.0')]  = '0'
                    cb_labs[(cb_labs=='1.0')]  = '1'
                else:
                    cpt = 0 ; cb_labs = []
                    for rr in vc_fld:
                        if cpt % fa.cb_jump == 0:
                            if fa.df >= 1.: cb_labs.append(str(int(rr)))
                            if fa.df <  1.: cb_labs.append(str(round(rr,int(np.ceil(np.log10(1./fa.df)))+1) ))
                        else:
                            cb_labs.append(' ')
                        cpt = cpt + 1
            #
            clb.set_label(fa.cunit, **fsm.cfont_clb)
            clb.ax.set_xticklabels(cb_labs, **fsm.cfont_clb_tcks)
            clb.ax.yaxis.set_tick_params(color=fa.color_top_cb) ; # set colorbar tick color
            clb.outline.set_edgecolor(fa.color_top_cb) ; # set colorbar edgecolor
            clb.outline.set_linewidth(0.3*rfz)
            clb.ax.tick_params(which = 'minor', length=1*rfz, width=0.3*rfz, color=fa.color_top_cb )
            clb.ax.tick_params(which = 'major', length=2*rfz, width=0.3*rfz, color=fa.color_top_cb )

        if nemo_box.l_show_clock:
            xl = float(x_clock)/rfz
            yl = float(y_clock)/rfz
            #ax.annotate('Date: '+cdats, xy=(1, 4), xytext=(xl,yl), **fsm.cfont_clock)
            ax.annotate(cdats, xy=(1, 4), xytext=(xl,yl), zorder=150, **fsm.cfont_clock)

        if nemo_box.l_show_exp:
            xl = float(x_exp)/rfz
            yl = float(y_exp)/rfz
            ax.annotate('Experiment: '+CNEMO+CRUN, xy=(1, 4), xytext=(xl,yl), **fsm.cfont_exp)

        if l_add_sign and nemo_box.l_show_sign:
            xl = float(x_sign)/rfz
            yl = float(y_sign)/rfz
            ax.annotate(CSIGN, xy=(1, 4), xytext=(xl,yl), zorder=150, **fsm.cfont_sign)

        if nemo_box.l_show_name:
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

        if cltr:
            ax.annotate(cltr+')', xy=(0.02, 0.93), xycoords='figure fraction', **fsm.cfont_letter )


        plt.savefig(cfig, dpi=rDPI, orientation='portrait', facecolor='k')
        print(cfig+' created!\n')
        plt.close(1)

        if fa.l_show_lsm: del clsm
        del cf, fig, ax
        if nemo_box.l_show_cb and lcb: del clb

    else:
        print('\n Figure '+cfig+' already there!\n')

id_f.close()


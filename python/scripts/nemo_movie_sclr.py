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


################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--fin' , required=True,                help='NEMO netCDF file to read from...')
requiredNamed.add_argument('-w', '--what', required=True, help='field/diagnostic to plot (ex: CSPEED,CURLOF,ect.)')

parser.add_argument('-C', '--conf', default="none",           help='name of NEMO config (ex: eNATL60) (defined into `nemo_hboxes.py`)')
parser.add_argument('-E', '--exp',  default="none",           help='name of experiment')
parser.add_argument('-b', '--box' , default="ALL",            help='extraction box name (ex: ALL) (defined into `nemo_hboxes.py`)')
parser.add_argument('-m', '--fmm' , default="mesh_mask.nc",   help='NEMO mesh_mask file (ex: mesh_mask.nc)')
parser.add_argument('-s', '--sd0' , default=None,       help='initial date as <YYYYMMDD>')
parser.add_argument('-l', '--lev' , type=int, default=0,      help='level to use if 3D field (default: 0 => 2D)')
parser.add_argument('-z', '--zld' ,                           help='topography netCDF file to use (field="z")')
parser.add_argument('-t', '--tstep',  default="1h",           help='time step ("1h","2h",..,up to "1d") in input file')
parser.add_argument('-N', '--oname',  default="",             help='a name that overides `CONF` on the plot...')
parser.add_argument('-o', '--outdir', default="./figs",       help='path to directory where to save figures')
parser.add_argument('-f', '--fignm',  default="",             help='common string in name of figure to create')
parser.add_argument('-T', '--addSST', default="",             help='add this SST field if showing a sea-ice field')
parser.add_argument('-S', '--sign',   default="",             help='sign the image with some text')

args = parser.parse_args()

CNEMO = args.conf
CEXP  = args.exp
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
CSIGN  = args.sign


print('')
print(' *** CNEMO = ', CNEMO)
print(' *** CBOX  = ', CBOX)
print(' *** CWHAT = ', CWHAT)
print(' *** cf_in = ', cf_in)
print(' *** cf_mm = ', cf_mm)

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



# SST under sea-ice field:
if l_add_SST_to_ice_field:
    cpal_sst = 'YlGnBu_r'
    rmin_sst = -2. ; rmax_sst = 22. ; dsst = 2.
    pal_sst = cp.chose_colmap(cpal_sst)
    norm_sst = colors.Normalize(vmin=rmin_sst, vmax=rmax_sst , clip = False)

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

if fa.imask_no_ice_pc > 0:
    # => going to read variable `fa.nm_ice_conc` into the same file 
    cf_ice = cf_in


    



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
    
    if fa.l_apply_lap:
        print(' *** Reading e1t and e2t !')
        XE1T2 = id_lsm.variables['e1t'][0,j1:j2,i1:i2]
        XE2T2 = id_lsm.variables['e2t'][0,j1:j2,i1:i2]
        XE1T2 = XE1T2*XE1T2
        XE2T2 = XE2T2*XE2T2
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
    
    if fa.l_apply_geov:        
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
    
    if cn_fig != "":
        cfig = cdir_figs+'/'+fa.cv_out+'_'+cn_fig+'_'+cdate+'.'+fig_type
    else:
        if l3d:
            cfig = cdir_figs+'/'+fa.cv_out+'_'+CNEMO+CRUN+'_lev'+str(jk)+'_'+CBOX+'_'+cdate+'.'+fig_type
        else:
            cfig = cdir_figs+'/'+fa.cv_out+'_'+CNEMO+CRUN+'_'+CBOX+'_'+cdate+'.'+fig_type


    if not path.exists(cfig):
    ###### FIGURE ##############

        fig = plt.figure(num = 1, figsize=(rw_fig, rh_fig), dpi=rDPI, facecolor='w', edgecolor='0.5')
    
        ax  = plt.axes([0., 0., 1., 1.], facecolor=fa.color_missing) # missing seas will be in 'facecolor' !
    
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

        if l_add_SST_to_ice_field:
            Xpsic = id_f.variables['siconc'][jt,j1:j2,i1:i2] ; # t, y, x  => we need it to separate open ocean from sea-ice !
            Xpsst = id_f.variables[caSST][jt,j1:j2,i1:i2] ; # t, y, x        
                
        print('Done!\n')

        if fa.rmult != 1.: Xplot[:,:] = fa.rmult * Xplot[:,:]
    
    
        # Ice
        if fa.l_show_ice or fa.imask_no_ice_pc>0:
            print('Reading record #'+str(jt)+' of '+fa.nm_ice_conc+' in '+cf_ice)
            id_ice = Dataset(cf_ice)
            XICE  = id_ice.variables[fa.nm_ice_conc][jt,j1:j2,i1:i2] ; # t, y, x
            id_ice.close()
            print('Done!\n')
    
    
        if fa.l_apply_lap:
            print(' *** Computing Laplacian of "'+fa.cv_in+'"!')
            lx = np.zeros((nj,ni))
            ly = np.zeros((nj,ni))
            lx[:,1:ni-1] = 1.E9*(Xplot[:,2:ni] -2.*Xplot[:,1:ni-1] + Xplot[:,0:ni-2])/XE1T2[:,1:ni-1]
            ly[1:nj-1,:] = 1.E9*(Xplot[2:nj,:] -2.*Xplot[1:nj-1,:] + Xplot[0:nj-2,:])/XE2T2[1:nj-1,:]
            Xplot[:,:] = lx[:,:] + ly[:,:]
            del lx, ly
    
        if fa.l_apply_hgrad:
            print(' *** Computing gradient of "'+fa.cv_in+'"!')
            lx = np.zeros((nj,ni))
            ly = np.zeros((nj,ni))
    
            if fa.l_smooth: cp.smoother(Xplot, XMSK, nb_smooth=fa.nb_smooth)
            
            # Zonal gradient on T-points:
            lx[:,1:ni-1] = (Xplot[:,2:ni] - Xplot[:,0:ni-2]) / (e1u[:,1:ni-1] + e1u[:,0:ni-2]) * UMSK[:,1:ni-1] * UMSK[:,0:ni-2]
            lx[:,:] = XMSK[:,:]*lx[:,:]
            #cp.dump_2d_field('dsst_dx_gridT.nc', lx, xlon=Xlon, xlat=Xlat, name='dsst_dx')
            # Meridional gradient on T-points:
            ly[1:nj-1,:] = (Xplot[2:nj,:] - Xplot[0:nj-2,:]) / (e2v[1:nj-1,:] + e2v[0:nj-2,:]) * VMSK[1:nj-1,:] * VMSK[0:nj-2,:]
            ly[:,:] = XMSK[:,:]*ly[:,:]
            #cp.dump_2d_field('dsst_dy_gridT.nc', ly, xlon=Xlon, xlat=Xlat, name='dsst_dy')
            Xplot[:,:] = 0.0
            # Modulus of vector gradient:        
            Xplot[:,:] = np.sqrt(  lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] )
            #cp.dump_2d_field('mod_grad_sst.nc', Xplot, xlon=Xlon, xlat=Xlat, name='dsst')
            del lx, ly
    
    
        if fa.l_apply_geov:
            print(' *** Computing gradient of "'+fa.cv_in+'"!')
            lx = np.zeros((nj,ni))
            ly = np.zeros((nj,ni))
    
            # Zonal gradient on T-points:
            lx[:,1:ni-1] = (Xplot[:,2:ni] - Xplot[:,0:ni-2]) / (e1u[:,1:ni-1] + e1u[:,0:ni-2]) * UMSK[:,1:ni-1] * UMSK[:,0:ni-2]
            lx[:,:] = XMSK[:,:]*lx[:,:]
            #cp.dump_2d_field('dsst_dx_gridT.nc', lx, xlon=Xlon, xlat=Xlat, name='dsst_dx')
            # Meridional gradient on T-points:
            ly[1:nj-1,:] = (Xplot[2:nj,:] - Xplot[0:nj-2,:]) / (e2v[1:nj-1,:] + e2v[0:nj-2,:]) * VMSK[1:nj-1,:] * VMSK[0:nj-2,:]
            ly[:,:] = XMSK[:,:]*ly[:,:]
            #cp.dump_2d_field('dsst_dy_gridT.nc', ly, xlon=Xlon, xlat=Xlat, name='dsst_dy')
            Xplot[:,:] = 0.0
            # Modulus of vector gradient:        
            Xplot[:,:] = grav/ff * np.sqrt( lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] )
            #cp.dump_2d_field('mod_grad_sst.nc', Xplot, xlon=Xlon, xlat=Xlat, name='dsst')
            del lx, ly

        if fa.lSubstractMean:
            zmean = np.sum( Xplot[:,:]*XMSK[:,:])/np.sum(XMSK[:,:])
            Xplot[:,:] = Xplot[:,:] - zmean
            
    
            
        print('')
        if not fa.l_show_lsm and jt == jt0: ( nj , ni ) = np.shape(Xplot)
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


        if fa.imask_no_ice_pc > 0:
            # masking field where too litle sea-ice:
            Xplot = np.ma.masked_where( XICE < float(fa.imask_no_ice_pc)/100., Xplot )
            

        cf = plt.imshow( Xplot[:,:], cmap=pal_fld, norm=norm_fld, interpolation=nemo_box.c_imshow_interp )
        #cf = plt.pcolormesh( Xplot[:,:], cmap=pal_fld, norm=norm_fld )

        # Add SST onto a sea-ice field:
        if l_add_SST_to_ice_field:
            psst = np.ma.masked_where(Xpsic > 0.05, Xpsst)
            ct   = plt.imshow(psst, cmap=pal_sst, norm=norm_sst, interpolation='none')
            del psst, ct
            
        # Add Sea-Ice onto a open ocean field:
        if fa.l_show_ice:
            pice = np.ma.masked_where(XICE < rmin_ice, XICE)
            ci = plt.imshow(pice, cmap=pal_ice, norm=norm_ice, interpolation='none') ; del pice, ci
    
        if fa.l_show_lsm or l_add_topo_land:
            if l_add_topo_land:
                clsm = plt.imshow( np.ma.masked_where(XLSM>0.0001, xtopo), cmap=pal_lsm, norm=norm_lsm, interpolation='none' )
                if nemo_box.c_imshow_interp == 'none':
                    plt.contour(XLSM, [0.9], colors='k', linewidths=0.5)
            else:
                pmsk = np.ma.masked_where(XLSM[:,:] > 0.2, XLSM[:,:])
                clsm = plt.imshow( pmsk, cmap=pal_lsm, norm=norm_lsm, interpolation='none' )                
                del pmsk
    
        ##### COLORBAR ######
        if nemo_box.l_show_cb:
            ax2 = plt.axes(nemo_box.vcb)
            if fa.l_pow_field or fa.l_log_field:
                clb = mpl.colorbar.ColorbarBase(ax=ax2, ticks=fa.vc_fld_powlog, cmap=pal_fld, norm=norm_fld,
                                                orientation='horizontal', extend=fa.cb_extend)
                #cb_labs = clb.ax.get_xticklabels()
                cb_labs = np.array( [ str(i) for i in fa.vc_fld_powlog ], dtype='U16' )
                #print(" fa.vc_fld_powlog = ", fa.vc_fld_powlog )
                #print(" cb_labs =", cb_labs) ; exit(0)
                cb_labs[(cb_labs=='0.0')]  = '0'
                cb_labs[(cb_labs=='1.0')]  = '1'                
            else:
                clb = mpl.colorbar.ColorbarBase(ax=ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld,
                                                orientation='horizontal', extend=fa.cb_extend)
                cb_labs = []
                cpt = 0
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
            ax.annotate('Date: '+cdats, xy=(1, 4), xytext=(xl,yl), **fsm.cfont_clock)
    
        if nemo_box.l_show_exp:
            xl = float(x_exp)/rfz
            yl = float(y_exp)/rfz
            ax.annotate('Experiment: '+CNEMO+CRUN, xy=(1, 4), xytext=(xl,yl), **fsm.cfont_exp)

        if l_add_sign and nemo_box.l_show_sign:
            xl = float(x_sign)/rfz
            yl = float(y_sign)/rfz
            ax.annotate(CSIGN, xy=(1, 4), xytext=(xl,yl), **fsm.cfont_sign)
    
        if nemo_box.l_show_name:
            cbla = CNEMO
            if CONAME != "": cbla = CONAME
            xl = float(x_name)/rfz
            yl = float(y_name)/rfz
            ax.annotate(cbla, xy=(1, 4), xytext=(xl, yl), **fsm.cfont_titl)
    
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

        if fa.l_show_lsm: del clsm
        del cf, fig, ax
        if nemo_box.l_show_cb: del clb

    else:
        print('\n Figure '+cfig+' already there!\n')

id_f.close()


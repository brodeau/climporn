#!/usr/bin/env python
#
#     CLIMPORN
#
#  Show SST + sea-ice concentration in the Arcti on a polar stereographic projection!
#  NEMO output + mesh_mask needed.
#
#    L. Brodeau, November 2019
#
import sys
from os import path, getcwd, mkdir
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

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import shiftgrid

from calendar import isleap
import datetime

from re import split
import warnings
warnings.filterwarnings("ignore")

import barakuda_colmap as bcm
import barakuda_tool as bt
import barakuda_ncio as bnc

ldrown = True

vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

fig_type='png'
rDPI = 100
color_top = 'white'
color_top_cb = 'white'
#color_top = 'k'

cv_out = 'sst_+_ice'

jt0 = 0


# SST
rmin_sst=-2.
rmax_sst=14.
dsst = 1. ; cb_jump = 1
#cpal_sst = 'ncview_nrl'
cpal_sst = 'on3'
cunit = r'$^{\circ}$C'

cdt = '6h'
l_get_name_of_run = True

# Ice:
rmin_ice=0.2
rmax_ice=0.98
#cpal_ice = 'ncview_bw'
#cpal_ice = 'Blues_r'
cpal_ice = 'bone'

vp =  ['nanuk', 'stere', -60., 40., 122., 57.,     75.,  -12., 10., 'h' ]  # North Pole


################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--fin' , required=True,                help='specify the NEMO netCDF file to read from...')
#requiredNamed.add_argument('-w', '--what', required=True, default="sst", help='specify the field/diagnostic to plot (ex: sst)')

parser.add_argument('-C', '--conf', default="NANUK025",     help='specify NEMO config (ex: eNATL60)')
#parser.add_argument('-b', '--box' , default="nanuk",        help='specify extraction box name (ex: ALL)')
parser.add_argument('-m', '--fmm' , default="mesh_mask.nc", help='specify the NEMO mesh_mask file (ex: mesh_mask.nc)')
parser.add_argument('-s', '--sd0' , default="19950101",     help='specify initial date as <YYYYMMDD>')
parser.add_argument('-l', '--lev' , type=int, default=0,    help='specify the level to use if 3D field (default: 0 => 2D)')
#parser.add_argument('-z', '--zld' ,                         help='specify the topography netCDF file to use (field="z")')

args = parser.parse_args()

CNEMO = args.conf
#CBOX  = args.box
#CWHAT = args.what
cf_in = args.fin
cf_mm = args.fmm
csd0  = args.sd0
jk    = args.lev

print ''
print ' *** CNEMO = ', CNEMO
#print ' *** CBOX  = ', CBOX
#print ' *** CWHAT = ', CWHAT
print ' *** cf_in = ', cf_in
print ' *** cf_mm = ', cf_mm
print ' *** csd0 = ', csd0
print ' ***   jk  = ', jk
l3d = False
if jk > 0:
    l3d=True
else:
    jk=0
    ###############################################################################################################################################

if not path.exists("figs"): mkdir("figs")
cdir_figs = './figs'
if not path.exists(cdir_figs): mkdir(cdir_figs)


# Name of RUN:
CRUN = ''
if l_get_name_of_run:
    vv = split('-|_', path.basename(cf_in))
    if vv[0] != CNEMO:
        print 'ERROR: your file name is not consistent with "'+CNEMO+'" !!! ('+vv[0]+')' ; sys.exit(0)
    CRUN = vv[1]
    print '\n Run is called: "'+CRUN+'" !\n'
#---------------------------------------------------------------
# Test - beta development:
#cf_in = '/data/gcm_output/CREG025/NANUK025-ILBOXE50_6h_gridT-2D_199506-199506.nc'
#cf_mm = '/data/gcm_setup/CREG025/CREG025-I/mesh_mask_CREG025_3.6_NoMed.nc'
if CNEMO == 'NANUK025':
    jk=0
    j1=0 ; j2=603
    i1=0 ; i2=528
else:
    print('ERRO: unknow conf '+CNEMO)
    ###############################


bt.chck4f(cf_mm)

l_notime=False
bt.chck4f(cf_in)
id_in = Dataset(cf_in)
list_var = id_in.variables.keys()
if 'time_counter' in list_var:
    vtime = id_in.variables['time_counter'][:]
elif 'time' in list_var:
    vtime = id_in.variables['time'][:]
else:
    l_notime=True
    print 'Did not find a time variable! Assuming no time and Nt=1'
    id_in.close()
    Nt = 1
if not l_notime: Nt = len(vtime)
print(' Nt ='), Nt



id_lsm = Dataset(cf_mm)
nb_dim = len(id_lsm.variables['tmask'].dimensions)
print(' The mesh_mask has '+str(nb_dim)+' dimmensions!')
print(' *** Reading '+'tmask'+' !')
if nb_dim==4: XMSK = id_lsm.variables['tmask'][0,jk,j1:j2,i1:i2]
if nb_dim==3: XMSK = id_lsm.variables['tmask'][jk,j1:j2,i1:i2]
if nb_dim==2: XMSK = id_lsm.variables['tmask'][j1:j2,i1:i2]
Xlon = id_lsm.variables['glamt'][0,j1:j2,i1:i2]
Xlat = id_lsm.variables['gphit'][0,j1:j2,i1:i2]
id_lsm.close()

(nj,ni) = nmp.shape(XMSK)  ; print('Shape Arrays => ni,nj ='), ni,nj

print('Done!\n')



fontr=1.2
params = { 'font.family':'Helvetica Neue',
           'font.weight':    'normal',
           'font.size':       int(9.*fontr),
           'legend.fontsize': int(9.*fontr),
           'xtick.labelsize': int(9.*fontr),
           'ytick.labelsize': int(9.*fontr),
           'axes.labelsize':  int(9.*fontr) }
mpl.rcParams.update(params)
cfont_clb  =  { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(8.*fontr), 'color':color_top_cb}
cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(9.*fontr), 'color':color_top }
#cfont_exp= { 'fontname':'Open Sans'  , 'fontweight':'light', 'fontsize':int(9.*fontr), 'color':color_top }
#cfont_mail  =  { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*fontr), 'color':'0.8'}
cfont_titl1 = { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(18.*fontr), 'color':color_top }
cfont_titl2 = { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(14.*fontr), 'color':color_top }


# for colorbar:
vc_sst = nmp.arange(rmin_sst, rmax_sst + dsst, dsst)



# Avec rebords:
#vfig_size = [ 10., 10.8 ]
#vsporg = [0.03, 0.1, 1., 0.8]

# For movie
vfig_size = [ 7.54, 7.2 ]
#vsporg = [0.001, 0.0011, 0.997, 0.999]
vsporg = [0., 0., 1., 1.]


vcbar = [0.05, 0.065, 0.92, 0.03]

cyr0=csd0[0:4]
cmn0=csd0[4:6]
cdd0=csd0[6:8]

if cdt == '6h':
    dt = 6
elif cdt == '3h':
    dt = 3
elif cdt == '1h':
    dt = 1
else:
    print 'ERROR: unknown dt!'

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
    cd = '%3.3i'%(jd)
    cm = '%2.2i'%(jm)
    ct = str(datetime.datetime.strptime(cyr0+'-'+cm+'-'+cd+' '+ch, '%Y-%m-%j %H'))
    ct=ct[:5]+cm+ct[7:] #lolo bug !!! need to do that to get the month and not "01"
    print ' ct = ', ct
    cday  = ct[:10]   ; print ' *** cday  :', cday
    chour = ct[11:13] ; print ' *** chour :', chour

    if l3d:
        cfig = cdir_figs+'/'+cv_out+'_'+CNEMO+'-'+CRUN+'_lev'+str(jk)+'_'+cday+'_'+chour+'.'+fig_type
    else:
        cfig = cdir_figs+'/'+cv_out+'_'+CNEMO+'-'+CRUN+'_'+cday+'_'+chour+'.'+fig_type

    # Getting SST and sea-ice concentration at time record "jt":
    id_in = Dataset(cf_in)
    XSST = id_in.variables['sst']   [jt,j1:j2,i1:i2]
    XIFR = id_in.variables['siconc'][jt,j1:j2,i1:i2]
    id_in.close()

    cjt = '%4.4i'%(jt)

    #######################################################################################################
    fig = plt.figure(num = 1, figsize=(vfig_size), dpi=None, facecolor='k', edgecolor='k')
    ax  = plt.axes(vsporg, facecolor = 'k')

    pal_sst = bcm.chose_colmap(cpal_sst)
    nrm_sst = colors.Normalize(vmin=rmin_sst, vmax=rmax_sst, clip=False)

    pal_ice = bcm.chose_colmap(cpal_ice)
    nrm_ice = colors.Normalize(vmin=rmin_ice, vmax=rmax_ice, clip = False)

    carte = Basemap(llcrnrlon=vp[2], llcrnrlat=vp[3], urcrnrlon=vp[4], urcrnrlat=vp[5], \
                    resolution=vp[9], area_thresh=1000., projection='stere', \
                    lat_0=vp[6], lon_0=vp[7])

    x0,y0 = carte(Xlon,Xlat)

    #f = carte.pcolor(x0, y0, XIFR, cmap = pal_ice, norm = nrm_ice) #, interpolation='none')

    if ldrown:
        print(' Drowning...')
        bt.drown(XIFR, XMSK, k_ew=-1, nb_max_inc=10, nb_smooth=10)
        bt.drown(XSST, XMSK, k_ew=-1, nb_max_inc=10, nb_smooth=10)

    idx_oce = nmp.where(XIFR  < rmin_ice)
    idx_ice = nmp.where(XIFR >= rmin_ice)
    XSST[idx_ice] = nmp.nan
    XIFR[idx_oce] = nmp.nan

    ft = carte.pcolormesh(x0, y0, XSST, cmap = pal_sst, norm = nrm_sst )
    fi = carte.pcolormesh(x0, y0, XIFR, cmap = pal_ice, norm = nrm_ice )
    #fi = carte.contour(x0, y0, XIFR, [ rmax_ice ], colors='k', linewidths=0.5 )

    carte.drawcoastlines(linewidth=0.5)
    carte.fillcontinents(color='grey')
    #carte.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    #carte.drawmapboundary()

    carte.drawmeridians(nmp.arange(-180,180,20), labels=[0,0,0,1])
    carte.drawparallels(nmp.arange( -90, 90,10), labels=[1,0,0,0])

    ax2 = plt.axes([0.64, 0.965, 0.35, 0.018])
    clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_sst, cmap=pal_sst, norm=nrm_sst, orientation='horizontal', extend='both')
    cb_labs = []
    #if cb_jump > 1:                                                                                                                                  
    cpt = 0
    for rr in vc_sst:
        if cpt % cb_jump == 0:
            if dsst >= 1.: cb_labs.append(str(int(rr)))
            if dsst <  1.: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./dsst)))+1) ))
        else:
            cb_labs.append(' ')
        cpt = cpt + 1
    clb.ax.set_xticklabels(cb_labs, **cfont_clb)
    clb.set_label(cunit, **cfont_clb)
    clb.ax.yaxis.set_tick_params(color=color_top_cb) ; # set colorbar tick color                                                                      
    clb.outline.set_edgecolor(color_top_cb) ; # set colorbar edgecolor                                                                                
    clb.ax.tick_params(which = 'minor', length = 2, color = color_top_cb )
    clb.ax.tick_params(which = 'major', length = 4, color = color_top_cb )


    
    ax.annotate('Date: '+cday+' '+chour+':00', xy=(0.785, 0.015), xycoords='figure fraction', **cfont_clock)
    ax.annotate('OPA - neXtSIM', xy=(0.02, 0.81 ), xycoords='figure fraction', **cfont_titl1)
    ax.annotate(' (CREG025) '  , xy=(0.07, 0.76 ), xycoords='figure fraction', **cfont_titl2)

    print(' Saving figure: '+cfig)
    plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)
    plt.close(1)

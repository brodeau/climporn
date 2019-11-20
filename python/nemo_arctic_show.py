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

import warnings
warnings.filterwarnings("ignore")



import barakuda_colmap as bcm

import barakuda_tool as bt
import barakuda_ncio as bnc

ldrown = True

fig_type='png'
rDPI = 100
color_top = 'white'
color_top_cb = 'white'
#color_top = 'k'


# SST
rmin_sst=-1.
rmax_sst=15.
#cpal_sst = 'ncview_nrl'
cpal_sst = 'on3'


# Ice:
rmin_ice=0.2
rmax_ice=0.98
#cpal_ice = 'ncview_bw'
#cpal_ice = 'Blues_r'
cpal_ice = 'bone'

vp =  ['nanuk', 'stere', -60., 40., 122., 57.,     75.,  -12., 10., 'h' ]  # North Pole



# Test - beta development:
cf_in = '/data/gcm_output/CREG025/NANUK025-ILBOXE50_6h_gridT-2D_199506-199506.nc'
cf_mm = '/data/gcm_setup/CREG025/CREG025-I/mesh_mask_CREG025_3.6_NoMed.nc'
jk=0
j1=0 ; j2=603
i1=0 ; i2=528
###############################

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



# Avec rebords:
#vfig_size = [ 10., 10.8 ]
#vsporg = [0.03, 0.1, 1., 0.8]

# For movie
vfig_size = [ 7.54, 7.2 ]
#vsporg = [0.001, 0.0011, 0.997, 0.999]
vsporg = [0., 0., 1., 1.]


vcbar = [0.05, 0.065, 0.92, 0.03]




for jt in range(Nt):

    # Getting SST and sea-ice concentration at time record "jt":
    id_in = Dataset(cf_in)
    XSST = id_in.variables['sst']   [jt,j1:j2,i1:i2]
    XIFR = id_in.variables['siconc'][jt,j1:j2,i1:i2]
    id_in.close()

    cjt = '%4.4i'%(jt)
    
    #######################################################################################################
    fig = plt.figure(num = 1, figsize=(vfig_size), dpi=None, facecolor='w', edgecolor='k')
    ax  = plt.axes(vsporg, facecolor = 'w')

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
    
    plt.savefig('test_'+cjt+'.png', dpi=rDPI, orientation='portrait', transparent=False)

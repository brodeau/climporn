#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#     CLIMPORN
#
#  Plot only the trajectories of all the buoys, even those who disapear
#
#  INPUT DATA: a `npz` file created with `trackice/scripts/traj2npz.py` (conversion from CSV to NPZ)
#
#    L. Brodeau, August 2022
#
# TO DO: use `nemo_box = cp.nemo_hbox(CNEMO,CBOX)` !!!
#
#
#  ABOUT input `npz` file:
#   * Name: should be of the form `NANUK4_ICE-BBM00_6h_19960101_19961031(_xxx).npz`
#
##################################################################

from sys import argv, exit
from os import path, mkdir
import numpy as np

from re import split

from netCDF4 import Dataset

import csv

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from calendar import isleap
import datetime

import climporn as cp

CBOX = 'ALL' ; # what box of `CCONF` ???
#CBOX = 'EastArctic' ; # what box of `CCONF` ???


idebug = 1

color_top = 'w'
color_top_cb = 'k'
clr_yellow = '#ffed00'

rDPI = 200

# Defaults:
l_scientific_mode = False

l_show_msh = False

fig_type='png'


narg = len(argv)
if not narg in [3,4]:
    print('Usage: '+argv[0]+' <file_trj.npz> <LSM_file> (iTsubsampl)')
    exit(0)

cf_npz = argv[1]
cf_lsm = argv[2]
# Subsampling in time...
itsubs = 1
if narg == 4 :
    itsubs = int(argv[3])

cp.chck4f(cf_npz)
cp.chck4f(cf_lsm)



cnfig  = str.replace( path.basename(cf_npz), '.npz', '.'+fig_type )


# Getting time info and time step from input npz file which is should look like NEMO output file:
vv = split('-|_', path.basename(cf_npz))

print(vv)


CCONF = vv[0]
print('\n *** CONF = '+CCONF)

dir_conf = path.dirname(cf_npz)
if dir_conf == '':  dir_conf = '.'
print('\n *** dir_conf =',dir_conf,'\n')


# Getting setup for NEMO CONFIG/BOX:
HBX = cp.nemo_hbox(CCONF,CBOX)
(i1,j1,i2,j2) = HBX.idx()

# About fields:
l_log_field = False
l_pow_field = False
cextend='both'
l_hide_cb_ticks = False
tmin=0. ; tmax=1. ; df=0.01
cb_jump = 1


bgclr = 'w'   ; # background color for ocean in figure

cfield = 'duration'
#fixme: end of october!
#tmin=0. ;  tmax=NbDays-1  ;  df = 30 ; # Arctic!
#cpal_fld = 'ncview_ice'
cpal_fld = 'viridis_r'
cunit = 'Time from seed (days)'
bgclr = 'w'   ; # background color for ocean in figure





if not path.exists(cf_npz):
    print('\n *** '+cf_npz+' NOT found !!!   :(  ')
    exit(0)

else:
    print('\n *** '+cf_npz+' was found !!!   :D  ')


if not path.exists("figs"): mkdir("figs")


#############################3
print('\n *** Reading into '+cf_npz+' !!!')
with np.load(cf_npz) as data:
    NrTraj = data['NrTraj']
    vtime   = data['time']
    xmask   = data['mask']
    xIDs    = data['IDs']
    xJIs    = data['JIs']
    xJJs    = data['JJs']
    xFFs    = data['FFs']



NbDays = int( (vtime[-1] - vtime[0]) / (3600.*24.) )
cdt1 = cp.epoch2clock(vtime[0] )
cdt2 = cp.epoch2clock(vtime[-1])

print('\n *** Start and End dates => '+cdt1+' -- '+cdt2)
print('      ==> nb of days =', NbDays)
    

tmin=0. ;  tmax=NbDays-1  ;  df = int(NbDays/10) ; # colorbar min max...

print('\n\n *** Trajectories contain '+str(NrTraj)+' records...')

cnmsk = 'tmask'
print('\n *** Reading "'+cnmsk+'" in meshmask file...')
id_lsm = Dataset(cf_lsm)
nb_dim = len(id_lsm.variables[cnmsk].dimensions)
Ni = id_lsm.dimensions['x'].size
Nj = id_lsm.dimensions['y'].size
#if i2 == 0: i2 = Ni
#if j2 == 0: j2 = Nj
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
#
nxr = int(HBX.rfact_zoom*nx_res) ; # widt image (in pixels)
nyr = int(HBX.rfact_zoom*ny_res) ; # height image (in pixels)
rh  = float(nxr)/float(rDPI) ; # width of figure as for figure...

print('\n *** width and height of image to create:', nxr, nyr, '\n')

pmsk    = np.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)
idx_oce = np.where(XMSK[:,:] > 0.5)


# Configuring fonts:
kk = cp.fig_style( HBX.font_rat, clr_top=color_top, clr_top_cb=color_top_cb )

# Colormaps for fields:
pal_fld = cp.chose_colmap(cpal_fld)
if l_log_field:
    norm_fld = colors.LogNorm(  vmin = tmin, vmax = tmax, clip = False)
if l_pow_field:
    norm_fld = colors.PowerNorm(gamma=pow_field, vmin = tmin, vmax = tmax, clip = False)
else:
    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

pal_lsm = cp.chose_colmap('land')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

vc_fld = np.arange(tmin, tmax + df, df)





# Only 1 figure
cfig = './figs/'+cnfig


fig = plt.figure(num=1, figsize=(rh,rh*yx_ratio), dpi=rDPI, facecolor='k', edgecolor='k')

if l_scientific_mode:
    ax1 = plt.axes([0.09, 0.09, 0.9, 0.9], facecolor = 'r')
else:
    ax1 = plt.axes([0., 0., 1., 1.],     facecolor = bgclr)




#ccm = plt.imshow(pmsk, cmap=pal_lsm, norm=norm_lsm, interpolation='none')
ccm = plt.pcolormesh( pmsk, cmap=pal_lsm, norm=norm_lsm )

plt.axis([ 0,i2-i1,0,j2-j1])




# If we show backward:
#print(xJIs[:,NrTraj-1])
#exit(0)

# Loop over time records:
icpt = -1
for jtt in range(NrTraj):
    #for jtt in np.arange(NrTraj-1,0,-1):
    icpt = icpt+1
    #    if jtt%itsubs == 0:
    print('jtt =',jtt)

    rfade = float(icpt)/float(NrTraj-1)*NbDays

    # Showing trajectories:
    csct = plt.scatter(xJIs[:,jtt]-i1, xJJs[:,jtt]-j1, c=xJIs[:,jtt]*0.+rfade, cmap=pal_fld , norm=norm_fld, marker='o', s=0.01) ; #, alpha=rfade ) ;#s=HBX.pt_sz_track ) ; # c=xFFs[:,jtt],






if l_scientific_mode:
    plt.xlabel('i-points', **cfont_axis)
    plt.ylabel('j-points', **cfont_axis)

if HBX.l_show_cb:
    ax2 = plt.axes(HBX.vcb)
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
    clb.set_label(cunit, **cp.fig_style.cfont_clb)
    clb.ax.yaxis.set_tick_params(color=color_top_cb) ; # set colorbar tick color
    #clb.outline.set_edgecolor(color_top_cb) ; # set colorbar edgecolor
    if l_hide_cb_ticks: clb.ax.tick_params(axis=u'both', which=u'both',length=0) ; # remove ticks!
    plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_top_cb) ; # set colorbar ticklabels


if HBX.l_show_name:  ax1.annotate(CCONF,          xy=(1, 4), xytext=HBX.name,  **cp.fig_style.cfont_ttl)

if HBX.l_show_exp:   ax1.annotate(CCONF,          xy=(1, 4), xytext=HBX.exp,   **cp.fig_style.cfont_ttl)

#if HBX.l_show_clock: ax1.annotate('Date: '+cdats, xy=(1, 4), xytext=HBX.clock, **cfont_clock)




#plt.savefig(cfig, dpi=rDPI, orientation='portrait', facecolor='b', transparent=True)
plt.savefig(cfig, dpi=rDPI, orientation='portrait') #, transparent=True)
print(cfig+' created!\n')
plt.close(1)

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
#  ABOUT input `npz` file:
#   * Name: should be of the form `NANUK4_ICE-BBM00_6h_19960101_19961031(_xxx).npz`
#
##################################################################

from sys import exit
from os import path, mkdir
import argparse as ap
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

#CBOX = 'ALL' ; # what box of `CCONF` ???
CBOX = 'ZoomArctic1' ; # what box of `CCONF` ???


idebug = 1

l_show_mod_field = False

color_top = 'w'
color_top_cb = 'k'
clr_yellow = '#ffed00'

rDPI = 200

# Defaults:
l_scientific_mode = False

l_show_msh = False

fig_type='png'

###################################################################################

################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--fin' , required=True,       help='`npz` file containing the csv comversion of TrackIce output')
requiredNamed.add_argument('-m', '--fmm' , required=True,       help='`mesh_mask` file of the NEMO config used...')
#
parser.add_argument(       '-j', '--fmd' , default="",          help='NEMO/SI3 output file that TrackIce used to build the npz file...')
parser.add_argument(       '-v', '--var' , default="siconc",    help='NEMO/SI3 name of the variable of field saved into npz file, to use for coloring buoys')
parser.add_argument(       '-s', '--tss' , type=int, default=1, help='temporal subsampling for generating figures')
#
args = parser.parse_args()
#
cf_npz = args.fin
cf_lsm = args.fmm
cf_mod = args.fmd
l_show_mod_field = ( cf_mod != "" )
cv_mod = args.var
itsubs = args.tss
################################################################################################

cp.chck4f(cf_npz)
if l_show_mod_field: cp.chck4f(cf_mod)
cp.chck4f(cf_lsm)

cnfig  = str.replace( path.basename(cf_npz), '.npz', '' )


# Getting time info and time step from input npz file which is should look like NEMO output file:
vv = split('-|_', path.basename(cf_npz))

print(vv)


CCONF = vv[0]
print('\n *** CONF = '+CCONF)

cdt   = vv[3]
print('\n *** Time frequency = '+cdt)

cdt1, cdt2 = vv[4], vv[5]
print('\n *** Start and End dates => '+cdt1+' -- '+cdt2)

#cyear = cdt1[0:4]
#print('\n *** Year = '+cyear+'\n')

NbDays = cp.Dates2NbDays(cdt1,cdt2)
print('     ==> number of days =', NbDays)

if cdt[-1]=='h':
    dt = int(cdt[0:-1]) ; ntpd = 24 / dt
    NbRecs = int(24/int(cdt[0:-1])*NbDays)
    if cdt=='1h': NbRecs = NbRecs+24 ; # fixme: why???? not fotr '6h' ???
else:
    print('ERROR: please adapt unit frequency: "'+cdt[-1]+'" !!!'); exit(0)
print('     ==> expected number of records from file name =', NbRecs)

dir_conf = path.dirname(cf_npz)
if dir_conf == '':  dir_conf = '.'
print('\n *** dir_conf =',dir_conf,'\n')

if l_show_mod_field: print('\n Field to show in background: "'+cv_mod+'" of file "'+cf_mod+'" !\n')

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

if l_show_mod_field: print('\n *** Model field to show in bacground: = '+cv_mod)

bgclr = 'w'   ; # background color for ocean in figure

if   cv_mod in ['sosstsst','tos']:
    cfield = 'SST'
    tmin=-3. ;  tmax=2.   ;  df = 0.1 ; # Arctic!
    #tmin=14. ;  tmax=26.   ;  df = 1.
    cpal_fld = 'inferno'
    cunit = r'SST ($^{\circ}$C)'

elif cv_mod in ['sosaline','sos']:
    cfield = 'SSS'
    tmin=32. ;  tmax=36.   ;  df = 1.
    cpal_fld = 'viridis'
    cunit = r'SSS (PSU)'

elif cv_mod in ['siconc']:
    cfield = cv_mod
    tmin=0. ;  tmax=1.   ;  df = 0.1 ; # Arctic!
    cpal_fld = 'ncview_ice'
    cunit = 'Ice concentration'
    bgclr = 'k'   ; # background color for ocean in figure

else:
    print('ERROR: variable '+cv_mod+' is not known yet...'); exit(0)


if not path.exists("figs"): mkdir("figs")


#############################3
print('\n *** Reading into '+cf_npz+' !!!')
with np.load(cf_npz) as data:
    vtime  = data['time'] ; # calendar in epoch time...
    Nrec = data['NbRec']
    xmask  = data['mask']
    xIDs   = data['IDs']
    xJIs   = data['JIs']
    xJJs   = data['JJs']
    xFFs   = data['FFs']

#for rt in vtime:    print(cp.epoch2clock(rt))


if Nrec != NbRecs-1:
    print('Warning: Nrec != NbRecs-1 !!!',Nrec,NbRecs-1)

print('\n\n *** Trajectories contain '+str(Nrec)+' records...')

if l_show_mod_field:
    with Dataset(cf_mod) as id_f_mod:
        Nt_mod = id_f_mod.dimensions['time_counter'].size - 1
    print('   => and '+str(Nt_mod)+' records of field '+cv_mod+' in NEMO file '+cf_mod+' !')
    if not Nrec%Nt_mod == 0:
        print('==> ERROR: they are not a multiple of each other!'); exit(0)
    nsubC = Nrec//Nt_mod
    print('    ==> number of subcycles for trajectories w.r.t model data =>', nsubC)

else:
    Nt_mod = Nrec
    nsubC  = 1


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
elif l_pow_field:
    norm_fld = colors.PowerNorm(gamma=pow_field, vmin = tmin, vmax = tmax, clip = False)
else:
    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

pal_lsm = cp.chose_colmap('land_dark')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

pal_filled = cp.chose_colmap('gray_r')
norm_filled = colors.Normalize(vmin = 0., vmax = 0.1, clip = False)

vc_fld = np.arange(tmin, tmax + df, df)


if l_show_mod_field:
    print('\n *** Opening file '+cf_mod)
    id_f_mod = Dataset(cf_mod)


    
# Loop over time records:

jtm = -1 ; # time record to use for model
l_read_mod = True

for jtt in range(Nrec):

    if jtt%nsubC == 0:
        jtm = jtm+1
        l_read_mod = True
    else:
        l_read_mod = False

    print( ' ### jtt, jtm = ',jtt, jtm )


    if jtt%itsubs == 0:

        cdate = cp.epoch2clock(vtime[jtt])

        # Only going if image not already present:
        cfig = './figs/'+cnfig+'_'+cdate+'.'+fig_type

        if not path.exists(cfig):

            fig = plt.figure(num=1, figsize=(rh,rh*yx_ratio), dpi=rDPI, facecolor='k', edgecolor='k')

            if l_scientific_mode:
                ax  = plt.axes([0.09, 0.09, 0.9, 0.9], facecolor = 'r')
            else:
                ax  = plt.axes([0., 0., 1., 1.],     facecolor = bgclr)


            if l_show_mod_field and l_read_mod:
                print('    => Reading record #'+str(jtm)+' of '+cv_mod+' in '+cf_mod)
                XFLD  = id_f_mod.variables[cv_mod][jtm-1,j1:j2,i1:i2] ; # t, y, x
                print('          Done!\n')

                if jtm == 0:
                    if XMSK[:,:].shape != XFLD.shape:
                        print('\n PROBLEM: field and mask do not agree in shape!')
                        print(XMSK.shape , XFLD.shape)
                        exit(0)
                    print('  *** Shape of field and mask => ', np.shape(XFLD))

            l_add_true_filled = False

            if l_show_mod_field:
                cf = plt.pcolormesh(XFLD[:,:], cmap=pal_fld, norm=norm_fld )

            if l_show_msh:
                ccx = plt.contour(Xlon[:,:], 60, colors='k', linewidths=0.5)
                ccy = plt.contour(Xlat[:,:], 30, colors='k', linewidths=0.5)


            #ccm = plt.imshow(pmsk, cmap=pal_lsm, norm=norm_lsm, interpolation='none')
            ccm = plt.pcolormesh( pmsk, cmap=pal_lsm, norm=norm_lsm )

            if l_add_true_filled: del pfilled


            plt.axis([ 0,i2-i1,0,j2-j1])

            # Showing trajectories:
            # Points have the color of the field!
            csct = plt.scatter(xJIs[:,jtt]-i1, xJJs[:,jtt]-j1, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=HBX.pt_sz_track )


            if l_scientific_mode:
                plt.xlabel('i-points', **cp.fig_style.cfont_axis)
                plt.ylabel('j-points', **cp.fig_style.cfont_axis)

            if l_show_mod_field and HBX.l_show_cb:
                ax2 = plt.axes(HBX.vcb)
                if l_pow_field or l_log_field:
                    clb = mpl.colorbar.ColorbarBase(ax=ax2,               cmap=pal_fld, norm=norm_fld,
                                                    orientation='horizontal', extend='neither')
                else:
                    clb = mpl.colorbar.ColorbarBase(ax=ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld,
                                                    orientation='horizontal', extend=cextend)
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
                clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color
                #clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor
                if l_hide_cb_ticks: clb.ax.tick_params(axis=u'both', which=u'both',length=0) ; # remove ticks!
                plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_top) ; # set colorbar ticklabels

                del csct

            if HBX.l_show_name:  ax.annotate(CCONF,          xy=(1, 4), xytext=HBX.name,  **cp.fig_style.cfont_ttl)

            if HBX.l_show_exp:   ax.annotate(CCONF,          xy=(1, 4), xytext=HBX.exp,   **cp.fig_style.cfont_ttl)

            if HBX.l_show_clock: ax.annotate('Date: '+cdate, xy=(1, 4), xytext=HBX.clock, **cp.fig_style.cfont_clck)





            #plt.savefig(cfig, dpi=rDPI, orientation='portrait', facecolor='b', transparent=True)
            plt.savefig(cfig, dpi=rDPI, orientation='portrait') #, transparent=True)
            print(cfig+' created!\n')
            plt.close(1)

            del ccm, fig, ax

        else:
            print('   ----- Figure '+cfig+' already there! -----\n')

# END OF LOOP !!!

if l_show_mod_field: id_f_mod.close()

#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#     CLIMPORN
#
#  Plot only the trajectories of all the buoys, even those who disapear
#
#    L. Brodeau, August 2022
#
# TO DO: use `nemo_box = cp.nemo_hbox(CNEMO,CBOX)` !!!
##################################################################

import sys
from os import path, mkdir
import numpy as nmp

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


i_field_plot = 7 ; # (C) column index of field to plot

idebug = 1

l_show_mod_field = False

color_top = 'w'
clr_yellow = '#ffed00'

rDPI = 200

# Defaults:
l_scientific_mode = False

l_show_msh = False

fig_type='png'

vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

#CBOX = 'ALL' ; # what box of `CCONF` ???
CBOX = 'EastArctic' ; # what box of `CCONF` ???

narg = len(sys.argv)
if not narg in [6,7]:
    print('Usage: '+sys.argv[0]+' <CONF> <file_trj.csv> <file_mod,var> <name_fig> <LSM_file> (iTsubsampl)')
    sys.exit(0)

CCONF  = sys.argv[1]
cf_trj = sys.argv[2]
vv = split(',',sys.argv[3])
cf_mod = vv[0]
cv_mod = vv[1]
cnfig  = sys.argv[4]
#
cf_lsm = sys.argv[5]
#
# Subsampling in time...
itsubs = 1
if narg == 7 :
    itsubs = int(sys.argv[6])



# Getting time info and time step from input model file:


vv = split('-|_', path.basename(cf_mod))

cyear = vv[-2]
cdt   = vv[-3]

print('\n *** Year = '+cyear)
print('\n *** time increment = '+cdt)

cyr0=cyear
yr0=int(cyear)
cmn0='01'
cdd0='01'

# Time step (h) as a string
if not len(cdt)==2:
    print('ERROR: something is wrong with the format of your time step => '+cdt+' !'); sys.exit(0)
if cdt=='1d':
    dt = 24 ; ntpd = 1
elif cdt[1]=='h':
    dt = int(cdt[0]) ; ntpd = 24 / dt
else:
    print('ERROR: something is wrong with the format of your time step => '+cdt+' !'); sys.exit(0)

vm = vmn
if isleap(yr0): vm = vml
#print(' year is ', vm, nmp.sum(vm)

jd = int(cdd0) - 1
jm = int(cmn0)

print('     ==> dt = ', dt,'h')


dir_conf = path.dirname(cf_trj)
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
    cfield = 'siconc'
    tmin=0. ;  tmax=1.   ;  df = 0.1 ; # Arctic!
    cpal_fld = 'ncview_ice'
    cunit = 'Ice concentration'
    bgclr = 'k'   ; # background color for ocean in figure

else:
    print('ERROR: variable '+cv_mod+' is not known yet...'); sys.exit(0)


if not path.exists("figs"): mkdir("figs")

#######################################################################################
# Testing, then reading CSV file
#######################################################################################

# A/ Scan the begining of the file to see how many trajectories have been initiated (before possible "killing" occurs...)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('\n *** Scanning #1 !')
with open(cf_trj, 'r') as ftxt:
    NbTrajInit = -1
    for line in csv.reader(ftxt, delimiter=','):
        iID = int(line[0])    # ID of current trajectory as an integer
        if iID < NbTrajInit: break # Then we are starting a new time record
        NbTrajInit = iID
print('      ===> number of initiated trajectories: = ',NbTrajInit)



# B/ Scan the end of the file to identify the `NbTrajEnd` trajectories (by their IDs) that made it to the end
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('\n *** Scanning #2 !')
LastStandIDs = []
with open(cf_trj, 'r') as ftxt:
    iID_o = 999999999
    for line in reversed(list(csv.reader(ftxt, delimiter=','))):
        iID = int(line[0])    # ID of current trajectory as an integer
        if iID > iID_o: break # Then we are starting a new time record
        LastStandIDs.append(iID)
        iID_o = iID
NbTrajEnd = len(LastStandIDs)
print('      ===> number of remaining trajectories at the end: = ',NbTrajEnd)
LastStandIDs = nmp.flipud(LastStandIDs) ; # back to increasing order + numpy array
#print('        ==> LastStandIDs = ', LastStandIDs[:], len(LastStandIDs) )

# C/ Scan the entire file to see how many time records are present
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('\n *** Scanning #3 !')
with open(cf_trj, 'r') as ftxt:
    iID_o = 999999999999
    NrTraj=0 ; # time record for the trajectories
    for line in csv.reader(ftxt, delimiter=','):
        iID = int(line[0])     ; # ID of current trajectory as an integer
        if iID < iID_o: NrTraj = NrTraj + 1 ; # This is the start of a new time record
        iID_o = iID
print('      ===> number of time records for the trajectories: = ', NrTraj)


# D/ Scan the entire file to store how many buoys are still alive at each of the NrTraj records
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('\n *** Scanning #4 !')
NbAlive = nmp.zeros(           NrTraj , dtype=int ) ; # number of buoys alive at each record
xIDs = nmp.zeros( (NbTrajInit, NrTraj), dtype=int ) ; # will mask dead IDs...
xJIs = nmp.zeros( (NbTrajInit, NrTraj), dtype=nmp.float32 ) ; # coordinates in terms of `ji` as float
xJJs = nmp.zeros( (NbTrajInit, NrTraj), dtype=nmp.float32 ) ; # coordinates in terms of `jj` as float
xFFs = nmp.zeros( (NbTrajInit, NrTraj), dtype=nmp.float32 ) ; # field #1 at position column #8 (7 in C)

with open(cf_trj, 'r') as ftxt:
    ID_o = -1
    jrec = 0
    Nliv = 0
    IDsR = []
    JIsR = []
    JJsR = []
    FFsR = []
    #
    for line in csv.reader(ftxt, delimiter=','):
        iID =   int(line[0])    # ID of current trajectory as an integer
        rJI = float(line[1])
        rJJ = float(line[2])
        rFF = float(line[i_field_plot])
        IDsR.append(iID)
        JIsR.append(rJI)
        JJsR.append(rJJ)
        FFsR.append(rFF)
        Nliv = Nliv + 1
        #
        if iID < ID_o: # Then we are starting a new time record
            Nliv = Nliv-1      ; # do not count this one as this is the first of the new record!
            ip = len(IDsR)-1   ; # index of last element
            iID_bkp = IDsR[ip] ; # we have appended one extra but we must remember for next record!!!
            rJI_bkp = JIsR[ip] ; # we have appended one extra but we must remember for next record!!!
            rJJ_bkp = JJsR[ip] ; # we have appended one extra but we must remember for next record!!!
            rFF_bkp = FFsR[ip] ; # we have appended one extra but we must remember for next record!!!
            #
            NbAlive[jrec] = Nliv
            #
            xIDs[0:Nliv,jrec] = nmp.array(IDsR[:-1])
            xJIs[0:Nliv,jrec] = nmp.array(JIsR[:-1])
            xJJs[0:Nliv,jrec] = nmp.array(JJsR[:-1])
            xFFs[0:Nliv,jrec] = nmp.array(FFsR[:-1])
            #
            # Preparin for next record, we have already everything for first value:
            jrec = jrec + 1    ; # we are already next record
            Nliv = 1           ; # That's number #1 of this "next" record
            IDsR = [ iID_bkp ]
            JIsR = [ rJI_bkp ]
            JJsR = [ rJJ_bkp ]
            FFsR = [ rFF_bkp ]
            #
        ID_o = iID
    #
#
# Values for last record:
NbAlive[   NrTraj-1] = NbTrajEnd     ;
xIDs[:Nliv,NrTraj-1] = LastStandIDs[:]
xJIs[:Nliv,NrTraj-1] = JIsR[:]
xJJs[:Nliv,NrTraj-1] = JJsR[:]
xFFs[:Nliv,NrTraj-1] = FFsR[:]


xJIs[:Nliv,NrTraj-1] = JIsR[:]
xJJs[:Nliv,NrTraj-1] = JJsR[:]

# Masking dead buoys:
xIDs = nmp.ma.masked_where(xIDs[:,:]==0, xIDs[:,:])
xJIs = nmp.ma.masked_where(xIDs[:,:]==0, xJIs[:,:])
xJJs = nmp.ma.masked_where(xIDs[:,:]==0, xJJs[:,:])
xFFs = nmp.ma.masked_where(xIDs[:,:]==0, xFFs[:,:])

if idebug > 1:
    print('\n Content of NbAlive:')
    for jrec in range(NrTraj):
        nba = NbAlive[jrec]
        print('###  Rec. # '+str(jrec)+' ==> '+str(nba)+' buoys alive!')
        print(' * IDs:\n',                       xIDs[:nba,jrec],'\n')
        print(' * Longitudes `ji` positions:\n', xJIs[:nba,jrec],'\n')
        print(' * Latitudes  `jj` positions:\n', xJJs[:nba,jrec],'\n')
        print(' * Values of field #1:\n',        xFFs[:nba,jrec],'\n')

print('   => done!\n')

# Time record stuff...
if l_show_mod_field:
    cp.chck4f(cf_mod)
    id_f_mod = Dataset(cf_mod)
    list_var = id_f_mod.variables.keys()
    if 'time_counter' in list_var:
        vtime = id_f_mod.variables['time_counter'][:]
        Nt_mod = len(vtime)
        print('\n There is a "time_counter" in file '+cf_mod+' !')
        print('   => '+str(Nt_mod)+' snapshots!')
    else:
        print('\nWARNING: there is NO "time_counter" in file '+cf_mod+' !')
        print('   ==> setting Nt_mod = 0 !\n')
        Nt_mod = 0
    id_f_mod.close()


print('\n\n *** Trajectories contain '+str(NrTraj)+' records each in CSV file')

if l_show_mod_field:
    print('   => and '+str(Nt_mod)+' records of field '+cv_mod+' in NEMO file '+cf_mod+' !')
    if not NrTraj%Nt_mod == 0:
        print('==> ERROR: they are not a multiple of each other!'); sys.exit(0)
    nsubC = NrTraj//Nt_mod
    print('    ==> number of subcycles for trajectories w.r.t model data =>', nsubC)

else:
    Nt_mod = NrTraj
    nsubC  = 1



cp.chck4f(cf_lsm)
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

pmsk    = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)
idx_oce = nmp.where(XMSK[:,:] > 0.5)

#font_rat
#params = { 'font.family':'Ubuntu',
params = { 'font.family':'Open Sans',
           'font.weight':    'normal',
           'font.size':       int(12.*HBX.font_rat),
           'legend.fontsize': int(22.*HBX.font_rat),
           'xtick.labelsize': int(18.*HBX.font_rat),
           'ytick.labelsize': int(18.*HBX.font_rat),
           'axes.labelsize':  int(15.*HBX.font_rat) }
mpl.rcParams.update(params)
cfont_clb   = { 'fontname':'Open Sans',   'fontweight':'medium', 'fontsize':int(18.*HBX.font_rat), 'color':color_top }
cfont_cnfn  = { 'fontname':'Open Sans',   'fontweight':'light' , 'fontsize':int(35.*HBX.font_rat), 'color':color_top }
cfont_axis  = { 'fontname':'Open Sans',   'fontweight':'medium', 'fontsize':int(18.*HBX.font_rat), 'color':color_top }
cfont_ttl   = { 'fontname':'Open Sans',   'fontweight':'medium', 'fontsize':int(25.*HBX.font_rat), 'color':color_top }
cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(19.*HBX.font_rat), 'color':color_top }

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

vc_fld = nmp.arange(tmin, tmax + df, df)




if l_show_mod_field:
    print('\n *** Opening file '+cf_mod)
    id_f_mod = Dataset(cf_mod)


# Loop over time records:

jtm = -1 ; # time record to use for model
l_read_mod = True
for jtt in range(NrTraj):

    if jtt%nsubC == 0:
        jtm = jtm+1
        l_read_mod = True
    else:
        l_read_mod = False

    print( ' ### jtt, jtm = ',jtt, jtm )

    ct   = '%4.4i'%(jtt+1)



    #---------------------- Calendar stuff --------------------------------------------
    jh   = (jtt*dt)%24
    #rjh  = ((float(jtt)+0.5)*dt)%24
    rjh  = ( float(jtt)*dt )%24
    if jtt%ntpd == 0: jd = jd + 1
    if jd == vm[jm-1]+1 and (jtt)%ntpd == 0 :
        jd = 1
        jm = jm + 1
        if jm==13:
            yr0 = yr0+1
            cyr0 = str(yr0)
            jm = 1
            vm = vmn
            if isleap(yr0): vm = vml

    ch  = '%2.2i'%(jh)
    crh = '%2.2i'%(rjh)
    cd  = '%3.3i'%(jd)
    cm  = '%2.2i'%(jm)
    #
    jhou = int(rjh)
    jmin = int((rjh-jhou)*60)
    chou = '%2.2i'%(jhou)
    cmin = '%2.2i'%(jmin)
    #
    ct  = str(datetime.datetime.strptime(cyr0+'-'+cm+'-'+cd+' '+ch, '%Y-%m-%j %H'))
    ct  = ct[:5]+cm+ct[7:] #lolo bug !!! need to do that to get the month and not '01
    cday  = ct[:10]   ; #print(' *** cday  :', cday
    if dt >= 24:
        cdate = cday
        cdats = cday
    else:
        chour = ct[11:13] ; #print(' *** chour :', chour
        cdate = cday+'_'+chour
        if jmin==0:
            cdats = cday+' '+chour+':00'
        else:
            cdats = cday+' '+chou+':'+cmin
    print('   ==> current date = ', cdats+' !')
    #-----------------------------------------------------------------------------------


    if jtt%itsubs == 0:


        # Only going if image not already present: lilo

        
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
                        sys.exit(0)
                    print('  *** Shape of field and mask => ', nmp.shape(XFLD))
    
            l_add_true_filled = False
    
            if l_show_mod_field:
                cf = plt.pcolormesh(XFLD[:,:], cmap=pal_fld, norm=norm_fld )
    
            if l_show_msh:
                ccx = plt.contour(Xlon[:,:], 60, colors='k', linewidths=0.5)
                ccy = plt.contour(Xlat[:,:], 30, colors='k', linewidths=0.5)
    
    
            #cm = plt.imshow(pmsk, cmap=pal_lsm, norm=norm_lsm, interpolation='none')
            cm = plt.pcolormesh( pmsk, cmap=pal_lsm, norm=norm_lsm )
    
            if l_add_true_filled: del pfilled
    
    
            plt.axis([ 0,i2-i1,0,j2-j1])
    
            # Showing trajectories:
            ct = plt.scatter(xJIs[:,jtt]-i1, xJJs[:,jtt]-j1, c=xFFs[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=HBX.pt_sz_track )
    
    
    
    
    
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
                clb.set_label(cunit, **cfont_clb)
                clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color
                #clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor
                if l_hide_cb_ticks: clb.ax.tick_params(axis=u'both', which=u'both',length=0) ; # remove ticks!
                plt.setp(plt.getp(clb.ax.axes, 'xticklabels'), color=color_top) ; # set colorbar ticklabels
    
                del ct
    
            if HBX.l_show_name:  ax.annotate(CCONF,          xy=(1, 4), xytext=HBX.name,  **cfont_cnfn)
    
            if HBX.l_show_exp:   ax.annotate(CCONF,          xy=(1, 4), xytext=HBX.exp,   **cfont_ttl)
    
            if HBX.l_show_clock: ax.annotate('Date: '+cdats, xy=(1, 4), xytext=HBX.clock, **cfont_clock)
    
    
            
    
    
            #plt.savefig(cfig, dpi=rDPI, orientation='portrait', facecolor='b', transparent=True)
            plt.savefig(cfig, dpi=rDPI, orientation='portrait') #, transparent=True)
            print(cfig+' created!\n')
            plt.close(1)
        
            del cm, fig, ax

        else:
            print('   ----- Figure '+cfig+' already there! -----\n')
            
# END OF LOOP !!!

if l_show_mod_field: id_f_mod.close()

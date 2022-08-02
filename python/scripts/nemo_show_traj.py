#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
##################################################################
#     CLIMPORN
#
#  Prepare 2D maps (monthly) that will later become a movie!
#  NEMO output and observations needed
#
#    L. Brodeau, August 2022
##################################################################

import sys
from os import path
import numpy as nmp

from re import split

from netCDF4 import Dataset

import csv

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import climporn as cp

idebug = 1
#followIDs = [ 1, 2, 3 ] ;# fixme # Trajectories to follow:
followIDs = [ 1, 2, 3, 4, 5, 6, 7 ] ;# fixme # Trajectories to follow:

l_show_mod_field = False




bathy_max = 5000. # m

color_top = 'w'
clr_yellow = '#ffed00'

# Defaults:
lknown = True
rfact_zoom = 1.
l_show_cb = False
l_show_nm = False
l_scientific_mode = False
l_show_ttl = False
vcb = [0.15, 0.96, 0.8, 0.02]
font_rat = 1.



l_show_msh = False
    
pt_sz_track = 30

fig_type='png'

narg = len(sys.argv)
if not narg in [6]:
    print('Usage: '+sys.argv[0]+' <CONF> <file_run_csv> <file_run,var> <name_fig> <LSM_file>')
    sys.exit(0)
#
# `file_run` => just to plot a field
    
CCONF  = sys.argv[1]
cf_csv = sys.argv[2]
vv = split(',',sys.argv[3])
cf_run = vv[0]
cv_run = vv[1]
cnfig  = sys.argv[4]
#
cf_lsm = sys.argv[5]

#cpnt = 't'
#if narg == 8 :
#    cpnt = sys.argv[7]

#if not cpnt in ['t','f','u','v']:
#    print('ERROR: what to do with C-grid "'+cpnt+'" point!?')
#    sys.exit(0)


    
print('\n Field to show in background: "'+cv_run+'" of file "'+cf_run+'" !\n')

dir_conf = path.dirname(cf_csv)
if dir_conf == '':  dir_conf = '.'
print('\n *** dir_conf =',dir_conf,'\n')

#######################################################################################
# Testing, then reading CSV file
#######################################################################################
ft = open( cf_csv, newline='' )
truc = csv.reader(ft, delimiter=',')

jl=0
jt=0
for line in truc:
    jl = jl + 1
    iID = int(line[0])      ; # ID of trajectory as an integer
    ctraj = '%4.4i'%(iID)  ; #       "      as a 4-char long string
    if iID == 1:
        jt = jt + 1
        if idebug > 1: print('\n record #',jt,':')
        # This is a new time record as we are dealing with first trajectory (again)
    
    if idebug > 1: print( 'Line:', jl, ' => traj #',iID, ' LINE => ', line )
    #print( line[0] , '\n')
Nrec_traj = jt
ft.close()


NbTraj = len(followIDs)

print('\n Number of trajectories to follow:', NbTraj)
print('   => with IDs:', followIDs )
print('   => along ',Nrec_traj,' integration time steps!')

ITRID = nmp.zeros( (NbTraj, Nrec_traj), dtype=int ) ; #trajectory ID (integer)
COORX = nmp.zeros( (NbTraj, Nrec_traj) )
COORY = nmp.zeros( (NbTraj, Nrec_traj) )
FLDO1 = nmp.zeros( (NbTraj, Nrec_traj) ) ; # first field at position column #8 (7 in C)

ft = open( cf_csv, newline='' )
truc = csv.reader(ft, delimiter=',')
print('\n')
jt=0
for line in truc:    
    iID = int(line[0])      ; # ID of trajectory as an integer
    if iID == 1:
        # This is a new time record as we are dealing with first trajectory (again)
        jt = jt + 1
        jtraj=0        
    if iID in followIDs:
        ITRID[jtraj,jt-1] = iID
        COORX[jtraj,jt-1] = float(line[1])
        COORY[jtraj,jt-1] = float(line[2])
        FLDO1[jtraj,jt-1] = float(line[7])
        jtraj = jtraj+1   # itteration of 1 trajectory for this particular record
ft.close()
        
# Debug, checking trajectories:
if idebug > 0:
    for jt in range(Nrec_traj):
        print('\n *** Record #',jt,':')
        for jtraj in range(NbTraj):
            print(' Traj. #',jtraj+1,' ==> ID =', ITRID[jtraj,jt], ' | x =', COORX[jtraj,jt], ' | y =', COORY[jtraj,jt] )
        
#######################################################################################
#######################################################################################

Nt_traj = Nrec_traj + 1 ; #lilo            


i2=0
j2=0

if CCONF == 'ORCA1':
    i1 = 0 ; j1 = 0 ; i2 = 362 ; j2 = 292 ; rfact_zoom = 3. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 1
    
else:
    print('\n WARNING [nemo_imshow_2d_field.py]: "'+CCONF+'" is an unknown config!\n     ==> falling back on default setup')
    lknown = False
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0



laplacian = False
l_log_field = False
l_pow_field = False
cextend='both'
l_hide_cb_ticks = False
tmin=0. ; tmax=1. ; df=0.01
cb_jump = 1


print(' cv_run = '+cv_run)

if   cv_run in ['sosstsst','tos']:
    cfield = 'SST'
    #tmin=-2. ;  tmax=32.   ;  df = 1.
    tmin=14. ;  tmax=26.   ;  df = 1.
    cpal_fld = 'ncview_nrl'    
    cunit = r'SST ($^{\circ}$C)'

elif cv_run in ['sosaline','sos']:
    cfield = 'SSS'
    tmin=32. ;  tmax=36.   ;  df = 1.
    cpal_fld = 'viridis'    
    cunit = r'SSS (PSU)'
    
else:
    print('ERROR: variable '+cv_run+' is not known yet...'); sys.exit(0)




# Time record stuff...
cp.chck4f(cf_run)
id_fld = Dataset(cf_run)
list_var = id_fld.variables.keys()
if 'time_counter' in list_var:    
    vtime = id_fld.variables['time_counter'][:]
    Nt_mod = len(vtime)
    print('\n There is a "time_counter" in file '+cf_run+' !')
    print('   => '+str(Nt_mod)+' snapshots!')
else:
    print('\nWARNING: there is NO "time_counter" in file '+cf_run+' !')
    print('   ==> setting Nt_mod = 0 !\n')
    Nt_mod = 0
id_fld.close()


print('\n\n *** So, we have '+str(Nt_traj)+' records for trajectories in CSV file')
print('   => and '+str(Nt_mod)+' records of field '+cv_run+' in NEMO file '+cf_run+' !')
#lilo

if not Nt_traj%Nt_mod == 0:
    print('==> ERROR: they are not a multiple of each other!'); sys.exit(0)
nsubC = Nt_traj//Nt_mod

print('    ==> number of subcycles for trajectories =', nsubC)

ibath=1

cp.chck4f(cf_lsm)
cnmsk = 'tmask'
print('\n *** Reading "'+cnmsk+'" in meshmask file...')
id_lsm = Dataset(cf_lsm)
nb_dim = len(id_lsm.variables[cnmsk].dimensions)
Ni = id_lsm.dimensions['x'].size
Nj = id_lsm.dimensions['y'].size
if i2 == 0: i2 = Ni
if j2 == 0: j2 = Nj
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
if not lknown:
    rfact_zoom = round(1000./float(ny_res),1)
nxr = int(rfact_zoom*nx_res) ; # widt image (in pixels)
nyr = int(rfact_zoom*ny_res) ; # height image (in pixels)
rDPI = 400
rh  = float(nxr)/float(rDPI) ; # width of figure as for figure...

print('\n *** width and height of image to create:', nxr, nyr, '\n')

pmsk    = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)
idx_oce = nmp.where(XMSK[:,:] > 0.5)

#font_rat
#params = { 'font.family':'Ubuntu',
params = { 'font.family':'Open Sans',
           'font.weight':    'normal',
           'font.size':       int(12.*font_rat),
           'legend.fontsize': int(22.*font_rat),
           'xtick.labelsize': int(18.*font_rat),
           'ytick.labelsize': int(18.*font_rat),
           'axes.labelsize':  int(15.*font_rat) }
mpl.rcParams.update(params)
cfont_clb  = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*font_rat), 'color':color_top }
cfont_date = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(12.*font_rat), 'color':'w' }
cfont_mail = { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*font_rat), 'color':'0.8'}
cfont_cnfn = { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(35.*font_rat), 'color':'w' }
cfont_axis  = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*font_rat), 'color':color_top }
cfont_ttl = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(25.*font_rat), 'color':color_top }


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

fsize  = ( rh, rh*yx_ratio )
vc_fld = nmp.arange(tmin, tmax + df, df)





print('\n *** Opening file '+cf_run)
id_fld = Dataset(cf_run)


# Loop over time records:

jtm = 0 ; # time record to use for model
l_read_mod = True
for jtt in range(Nt_traj):

    
    if jtt%nsubC == 0:
        jtm = jtm+1
        l_read_mod = True
    else:
        l_read_mod = False
        

    print( ' ### jtt, jtm = ',jtt, jtm )
        
    ct   = '%4.4i'%(jtt+1)
    cfig = cnfig+'_'+ct+'.'+fig_type

    fig = plt.figure(num = 1, figsize=fsize, dpi=rDPI, facecolor='k', edgecolor='k')

    if l_scientific_mode:
        ax  = plt.axes([0.09, 0.09, 0.9, 0.9], facecolor = 'r')
    else:
        ax  = plt.axes([0., 0., 1., 1.],     facecolor = '0.9')


    if l_show_mod_field and l_read_mod:
        print('    => Reading record #'+str(jtm)+' of '+cv_run+' in '+cf_run)
        XFLD  = id_fld.variables[cv_run][jtm-1,j1:j2,i1:i2] ; # t, y, x
        print('          Done!\n')

        if jtm == 0:
            if XMSK[:,:].shape != XFLD.shape:
                print('\n PROBLEM: field and mask do not agree in shape!')
                print(XMSK.shape , XFLD.shape)
                sys.exit(0)
            print('  *** Shape of field and mask => ', nmp.shape(XFLD))
    
    
    
    l_add_true_filled = False
    
    #if cfield == 'Bathymetry':
    #    (idy_nan,idx_nan) = nmp.where( nmp.isnan(XFLD) )
    #    #
    #    # LSM with different masking for true lsm and filled lsm...
    #    cf_mask_lbc = dir_conf+'/lsm_LBC_'+CCONF+'.nc'
    #    if path.exists(cf_mask_lbc):
    #        print('\n *** '+cf_mask_lbc+' found !!!')
    #        l_add_true_filled = True 
    #        id_filled = Dataset(cf_mask_lbc)
    #        xtmp = id_filled.variables['lsm'][j1:j2,i1:i2]
    #        id_filled.close()
    #        pfilled = nmp.ma.masked_where(xtmp[:,:] != -1., xtmp[:,:]*0.+40.)
    #        del xtmp
    #
    #        print('  => done filling "pfilled" !\n')
    
        
    #if cv_run == 'track':
    #    
    #    XFLD[nmp.where(nmp.isnan(XFLD))] = -1000
    #    indx = nmp.where( XFLD > 0 )
    #    (idy,idx) = indx
    #    
    #    tmin=nmp.amin(XFLD[indx]) ;  tmax=nmp.amax(XFLD[indx])
    #    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)
    # 
    #    cf = plt.scatter(idx, idy, c=XFLD[indx], cmap=pal_fld, norm=norm_fld, alpha=0.5, marker='.', s=pt_sz_track )
    #
    #else:
    
    #cf = plt.imshow(XFLD[:,:], cmap=pal_fld, norm=norm_fld, interpolation='none')
    if l_show_mod_field:
        cf = plt.pcolormesh(XFLD[:,:], cmap=pal_fld, norm=norm_fld )
            
    if l_show_msh:
        ccx = plt.contour(Xlon[:,:], 60, colors='k', linewidths=0.5)
        ccy = plt.contour(Xlat[:,:], 30, colors='k', linewidths=0.5)
                
    
    #cm = plt.imshow(pmsk, cmap=pal_lsm, norm=norm_lsm, interpolation='none')
    cm = plt.pcolormesh( pmsk, cmap=pal_lsm, norm=norm_lsm )
    
    if l_add_true_filled: del pfilled
    
    
    plt.axis([ 0, Ni, 0, Nj])

    # Showing trajectories:
    ##ct = plt.scatter([Ni/2,Ni/3], [Nj/2,Nj/3], c=XFLD[indx], cmap=pal_fld, norm=norm_fld, alpha=0.5, marker='.', s=pt_sz_track )
    #ct = plt.scatter([Ni/2,Ni/3], [Nj/2,Nj/3], cmap=pal_fld, norm=norm_fld, alpha=0.5, marker='.', s=pt_sz_track )

    #ct = plt.scatter(COORX[:,jtt], COORY[:,jtt], color='r', marker='.', s=pt_sz_track )
    ct = plt.scatter(COORX[:,jtt], COORY[:,jtt], c=FLDO1[:,jtt], cmap=pal_fld, norm=norm_fld, marker='.', s=pt_sz_track )



    
    if l_scientific_mode:
        plt.xlabel('i-points', **cfont_axis)
        plt.ylabel('j-points', **cfont_axis)
    
    if l_show_cb:
    
        ax2 = plt.axes(vcb)
    
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
            
        del cf
        
    if l_show_nm:  ax.annotate(CCONF, xy=(1, 4), xytext=(x_cnf, y_cnf), **cfont_cnfn)
    
    if l_show_ttl: ax.annotate(CCONF, xy=(1, 4), xytext=(x_ttl, y_ttl), **cfont_ttl)
    
    
    
    #plt.savefig(cfig, dpi=rDPI, orientation='portrait', facecolor='b', transparent=True)
    plt.savefig(cfig, dpi=rDPI, orientation='portrait') #, transparent=True)
    print(cfig+' created!\n')
    plt.close(1)



    del cm, fig, ax
# END OF LOOP !!!

id_fld.close()

#!/usr/bin/env python

#       B a r a K u d a
#
#  Prepare 2D maps (monthly) that will later become a GIF animation!
#  NEMO output and observations needed
#
#    L. Brodeau, May 2018

import sys
import os
from string import replace
import numpy as nmp

from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import warnings
warnings.filterwarnings("ignore")

import clprn_colmap as bcm

import clprn_tool as bt

bathy_max = 5000. # m

color_top = 'white'
clr_yellow = '#ffed00'


l_show_cb = True
l_show_nm = True
l_scientific_mode = False
l_show_ttl = False

l_show_msh = False
    
pt_sz_track = 20

l_read_lsm=False


fig_type='png'

narg = len(sys.argv)
if not narg in [5,6]:
    print 'Usage: '+sys.argv[0]+' <CONF> <file> <variable> <snapshot> (<LSM_file>)'
    sys.exit(0)
CNEMO  = sys.argv[1]
cf_fld = sys.argv[2]
cv_in  = sys.argv[3]
jt=int(sys.argv[4])

if narg ==6 :
    l_read_lsm=True
    cf_lsm = sys.argv[5]


l_bathy_var = [ 'Bathymetry', 'elevation' ]
    

if not l_read_lsm and ( not cv_in in l_bathy_var):
    print "It's only for bathymetric fields that you can skip providing the mesh_mask file!"
    sys.exit(0)




    
i2=0
j2=0

if CNEMO == 'NATL60':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.25 ; vcb = [0.6, 0.1, 0.39, 0.025] ; font_rat = 5.*rfact_zoom
    x_cnf = 160. ; y_cnf = 2300. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False
    
elif CNEMO == 'NANUK1':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 4. ; vcb = [0.5, 0.875, 0.49, 0.02] ; font_rat = 0.16*rfact_zoom

elif CNEMO == 'NANUK025':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.5, 0.875, 0.49, 0.02] ; font_rat = 0.5*rfact_zoom
    x_cnf = 30. ; y_cnf = 540. ; # where to put label of conf on Figure...

elif CNEMO == 'ROALD12':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.63, 0.95, 0.36, 0.02] ; font_rat = 1.*rfact_zoom
    x_cnf = 50. ; y_cnf = 1250. ; # where to put label of conf on Figure...
    
elif CNEMO in [ 'CREG025' ] :
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2.
    vcb = [0.6, 0.975, 0.38, 0.02] ; font_rat = 0.5*rfact_zoom
    x_cnf = 20. ; y_cnf = 560. ; # where to put label of conf on Figure...

elif CNEMO in [ 'CREG4' ] :
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2.
    vcb = [0.14, 0.05, 0.8, 0.02] ; font_rat = 0.5*rfact_zoom
    x_ttl = 210. ; y_ttl = 620. ; # where to put label of conf on Figure...
    l_show_nm = False ; l_show_msh = True
    l_scientific_mode = True ; l_show_ttl = True
    color_top = 'k'

elif CNEMO == 'eNATL4':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.25 ; vcb = [0.6, 0.11, 0.39, 0.025] ; font_rat = 0.5*rfact_zoom
    x_cnf = 20. ; y_cnf = 560. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False

elif CNEMO == 'eNATL36':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.3 ; vcb = [0.6, 0.1, 0.39, 0.025] ; font_rat = 2.5*rfact_zoom
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False

elif CNEMO == 'eNATL60':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 0.25 ; vcb = [0.6, 0.1, 0.39, 0.025] ; font_rat = 5.*rfact_zoom
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False

elif CNEMO == 'KANAK60':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.5, 0.1, 0.45, 0.025] ; font_rat = 1.*rfact_zoom
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False

elif CNEMO == 'eNATL1':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 6.
    vcb = [0.62, 0.11, 0.35, 0.025] ; font_rat = 0.12*rfact_zoom
    x_cnf = 4. ; y_cnf = 120. ; # where to put label of conf on Figure...

elif CNEMO == 'SouthPac':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.5, 0.875, 0.49, 0.02] ; font_rat = 0.16*rfact_zoom
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False
    bathy_max = 8000. # m
    
elif CNEMO == 'SWEPAC2':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 12. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 1./rfact_zoom*10.
    x_cnf = 160. ; y_cnf = 4000. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = False
    bathy_max = 6000. # m
    
elif CNEMO == 'TROPICO2':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 10. ; vcb = [0.35, 0.94, 0.6, 0.04] ; font_rat = 8./rfact_zoom
    x_cnf = 20. ; y_cnf = 3. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = True ; l_scientific_mode=False
    bathy_max = 6000. # m

elif CNEMO == 'TROPICO05':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 3. ; vcb = [0.02, 0.15, 0.4, 0.04] ; font_rat = 2./rfact_zoom
    x_cnf = 280. ; y_cnf = 135. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = True ; l_scientific_mode=False
    bathy_max = 6000. # m
    
elif CNEMO == 'TROPICO12':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.02, 0.15, 0.4, 0.04] ; font_rat = 2./rfact_zoom
    x_cnf = 1200. ; y_cnf = 800. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = True ; l_scientific_mode=False
    bathy_max = 6000. # m

elif CNEMO == 'GEBCO':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 1. ; vcb = [0.2, 0.06, 0.6, 0.03] ; font_rat = 1.5/rfact_zoom
    x_cnf = 900. ; y_cnf = 1350. ; # where to put label of conf on Figure...
    l_show_cb = True ; l_show_nm = True ; l_scientific_mode=False
    bathy_max = 5000. # m
    
elif CNEMO == 'Azores':
    i1 = 0 ; j1 = 0 ; i2 = 0 ; j2 = 0 ; rfact_zoom = 2. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 1./rfact_zoom*10.
    x_cnf = 0. ; y_cnf = 0. ; # where to put label of conf on Figure...
    l_show_cb = False ; l_show_nm = False
    bathy_max = 6000. # m

elif CNEMO == 'GulfS':
    i1 = 0 ; j1 = 420 ; i2 = 900 ; j2 = 980 ; rfact_zoom = 1. ; vcb = [0.15, 0.96, 0.8, 0.02] ; font_rat = 0.1
    l_show_cb = False ; l_show_nm = False
    pt_sz_track = 3
    
else:
    print '\n PROBLEM: "'+CNEMO+'" is an unknown config!!!'
    sys.exit(0)








laplacian = False
l_log_field = False
l_pow_field = False
cextend='both'
l_hide_cb_ticks = False

print ' cv_in = '+cv_in

if cv_in in ['sosstsst','tos']:
    cfield = 'SST'
    tmin=0. ;  tmax=28.   ;  df = 1.
    cpal_fld = 'ncview_nrl'    
    cunit = r'SST ($^{\circ}$C)'
    cb_jump = 1

elif cv_in in l_bathy_var:   # 
    cfield = 'Bathymetry'       # 
    tmin=0. ;  tmax=bathy_max   ;  df = 100.
    #cpal_fld = 'ocean'
    #cpal_fld = 'Blues'
    #cpal_fld = 'PuBu'
    cpal_fld = 'on2_r'
    #cpal_fld = 'ncview_ssec'
    #cpal_fld = 'ncview_hotres'
    #cpal_fld = 'ncview_helix'
    cunit = r'Bathymetry (m)'
    cb_jump = 10
    l_pow_field = True
    pow_field = 1.5
    #l_log_field = False
    cextend='max'
    l_hide_cb_ticks=True
    if cv_in == 'elevation':
        #cpal_fld = 'ncview_hotres' ; l_pow_field = False
        cpal_fld = 'ncview_nrl' ; l_pow_field = True
        #cpal_fld = 'ocean' ; l_pow_field = False
        

    
elif cv_in == 'sossheig':
    cfield = 'SSH'
    tmin=-0.3 ;  tmax=0.3   ;  df = 0.1
    cpal_fld = 'ncview_jaisnc'    
    cunit = r'SSH (m)'
    cb_jump = 1

elif cv_in == 'somxl010':
    cfield == 'MLD'
    tmin=50. ;  tmax=1500. ;  df = 50.
    cpal_fld = 'viridis_r'

elif cv_in == 'track':
    cfield = 'TRACK'
    tmin=140. ;  tmax=24800.   ;  df = 1.
    cpal_fld = 'viridis'    
    cunit = r'SST ($^{\circ}$C)'
    cb_jump = 1
    fig_type='svg'

else:
    print 'ERROR: variable '+cv_in+' is not known yet...'; sys.exit(0)




# Time record stuff...
bt.chck4f(cf_fld)
id_fld = Dataset(cf_fld)
list_var = id_fld.variables.keys()
if 'time_counter' in list_var:    
    vtime = id_fld.variables['time_counter'][:]
    Nt = len(vtime)
    print '\n There is a "time_counter" in file '+cf_fld+' !'
    print '   => '+str(Nt)+' snapshots!'
    if jt <= 0 or jt > Nt: print ' PROBLEM: your time record does not exist!', jt ; sys.exit(0)
else:
    print '\nWARNING: there is NO "time_counter" in file '+cf_fld+' !'
    print '   ==> setting Nt = 0 !\n'
    Nt = 0
id_fld.close()

ibath=1

if l_read_lsm:
    bt.chck4f(cf_lsm)
    print '\n *** Reading "tmask" in meshmask file...'
    id_lsm = Dataset(cf_lsm)
    nb_dim = len(id_lsm.variables['tmask'].dimensions)
    Ni = id_lsm.dimensions['x'].size
    Nj = id_lsm.dimensions['y'].size
    if i2 == 0: i2 = Ni
    if j2 == 0: j2 = Nj
    if nb_dim == 4: XMSK  = id_lsm.variables['tmask'][0,0,j1:j2,i1:i2] ; # t, y, x
    if nb_dim == 3: XMSK  = id_lsm.variables['tmask'][0,  j1:j2,i1:i2] ; # t, y, x
    if nb_dim == 2: XMSK  = id_lsm.variables['tmask'][    j1:j2,i1:i2] ; # t, y, x
    if l_show_msh:
        Xlon = id_lsm.variables['glamu'][0,j1:j2,i1:i2]
        Xlat = id_lsm.variables['gphiv'][0,j1:j2,i1:i2]
    id_lsm.close()
    print '      done.'

elif cv_in in l_bathy_var:
    bt.chck4f(cf_fld)
    print '\n *** Will build mask from "Bathymetry"...'
    id_fld = Dataset(cf_fld)
    list_dim = id_fld.dimensions.keys()
    nb_dim = len(id_fld.variables[cv_in].dimensions)
    for jd in range(nb_dim):
        if list_dim[jd] in ['x','lon','longitude','glamt']: cdim_x = list_dim[jd]
        if list_dim[jd] in ['y','lat','latitude' ,'gphit']: cdim_y = list_dim[jd]
    print ' *** x, y dims =>', cdim_x, cdim_y
    #
    Ni = id_fld.dimensions[cdim_x].size
    Nj = id_fld.dimensions[cdim_y].size
    if i2 == 0: i2 = Ni
    if j2 == 0: j2 = Nj
    if nb_dim == 3: XBATH = id_fld.variables[cv_in][0,  j1:j2,i1:i2] ; # t, y, x
    if nb_dim == 2: XBATH = id_fld.variables[cv_in][    j1:j2,i1:i2] ; # t, y, x
    if l_show_msh:
        Xlon = id_fld.variables['nav_lon'][j1:j2,i1:i2]
        Xlat = id_fld.variables['nav_lat'][j1:j2,i1:i2]
    id_fld.close()

    # Does bathymetry come as positive or negative in file???
    bt_max = nmp.max(XBATH[:,:])
    bt_min = nmp.min(XBATH[:,:])
    print 'min, max bathy =>', bt_max, bt_min

    if abs(bt_max) > 12000. or abs(bt_min) > 12000.:
        print 'PROBLEM #1: we dont know what do do with min and max of bathymetry!'; sys.exit(0)
    if abs(bt_min) > abs(bt_max):
        print ' *** Bathymetry seems negative!'
        ibath = -1
        XBATH[:,:] = XBATH[:,:]*ibath
    
    XMSK = nmp.zeros((Nj,Ni))
    idx_oce = nmp.where(XBATH > 0.2)
    XMSK[idx_oce] = 1.
    del XBATH
    print '      done.'

    
else:
    print 'PROBLEM #2'; sys.exit(0)


print '\n According to "tmask" the shape of the domain is Ni, Nj =', Ni, Nj


# Stuff for size of figure respecting pixels...
print '  *** we are going to show: i1,i2,j1,j2 =>', i1,i2,j1,j2, '\n'
nx_res = i2-i1
ny_res = j2-j1
yx_ratio = float(ny_res)/float(nx_res)
nxr = int(rfact_zoom*nx_res) ; # widt image (in pixels)
nyr = int(rfact_zoom*ny_res) ; # height image (in pixels)
dpi = 100
rh  = float(nxr)/float(dpi) ; # width of figure as for figure...


print '\n *** width and height of image to create:', nxr, nyr, '\n'




pmsk = nmp.ma.masked_where(XMSK[:,:] > 0.2, XMSK[:,:]*0.+40.)





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
cfont_cnfn = { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(50.*font_rat), 'color':'w' }
cfont_axis  = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(18.*font_rat), 'color':color_top }
cfont_ttl = { 'fontname':'Open Sans', 'fontweight':'medium', 'fontsize':int(25.*font_rat), 'color':color_top }


# Colormaps for fields:
pal_fld = bcm.chose_colmap(cpal_fld)
if l_log_field:
    norm_fld = colors.LogNorm(  vmin = tmin, vmax = tmax, clip = False)
if l_pow_field:
    norm_fld = colors.PowerNorm(gamma=pow_field, vmin = tmin, vmax = tmax, clip = False)
else:
    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)

pal_lsm = bcm.chose_colmap('land_dark')
norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)


if Nt == 0:
    cfig = cv_in+'_'+CNEMO+'_'+cpal_fld+'.'+fig_type    
else:
    cfig = 'snapshot_'+str(jt)+'_'+cv_in+'_'+CNEMO+'_'+cpal_fld+'.'+fig_type    


#rextra_height = 1.
#if l_scientific_mode: rextra_height = 1.12
#fig = plt.figure(num = 1, figsize=(rh,rh*yx_ratio*rextra_height), dpi=None, facecolor='w', edgecolor='0.5')

fig = plt.figure(num = 1, figsize=(rh,rh*yx_ratio), dpi=None, facecolor='w', edgecolor='k')

if l_scientific_mode:
    ax  = plt.axes([0.09, 0.09, 0.9, 0.9], axisbg = 'r')
else:
    ax  = plt.axes([0., 0., 1., 1.],     axisbg = 'k')

vc_fld = nmp.arange(tmin, tmax + df, df)


print '\n *** Opening file '+cf_fld
id_fld = Dataset(cf_fld)
if Nt > 0:
    print '    => Reading record #'+str(jt)+' of '+cv_in+' in '+cf_fld
    XFLD  = id_fld.variables[cv_in][jt-1,j1:j2,i1:i2] ; # t, y, x
else:
    print '    => Reading 2D field '+cv_in+' in '+cf_fld+' (no time records...)'
    XFLD  = id_fld.variables[cv_in][j1:j2,i1:i2] ; # t, y, x

id_fld.close()
print '          Done!\n'


if XMSK.shape != XFLD.shape:
    print '\n PROBLEM: field and mask do not agree in shape!'
    print XMSK.shape , XFLD.shape
    sys.exit(0)

print '  *** Shape of field and mask => ', nmp.shape(XFLD)

del XMSK


if cv_in in l_bathy_var and ibath==-1: XFLD = ibath*XFLD


(idy_nan,idx_nan) = nmp.where( nmp.isnan(XFLD) )


print 'Ploting'

if cv_in == 'track':
    indx = nmp.where( XFLD > 0 )
    (idy,idx) = indx
    cf = plt.scatter(idx, idy, c=XFLD[indx], cmap = pal_fld, norm = norm_fld, alpha=0.5, marker='.', s=pt_sz_track )
    #
else:
    cf = plt.imshow(XFLD[:,:], cmap = pal_fld, norm = norm_fld, interpolation='nearest' ) #, interpolation='none')
    if len(idy_nan) > 0:
        idd = nmp.where(idy_nan==1); idy_nan[idd] = int(10./rfact_zoom)/2  # just so the boundary line is not too thin on plot...
        plt.scatter(idx_nan, idy_nan, color=clr_yellow, marker='s', s=int(10./rfact_zoom))

if l_show_msh:
    ccx = plt.contour(Xlon[:,:], 60, colors='k', linewidths=0.5)
    ccy = plt.contour(Xlat[:,:], 30, colors='k', linewidths=0.5)

#vj = nmp.arange(0,Nj-20,20)
#for jj in vj:
#    ccx = plt.plot(nmp.arange(len(Xlat[jj,:])), 5.*Xlat[jj,:], 'k', linewidth=0.5)

#vi = nmp.arange(0,Ni-20,20)
#for ii in vi:
#    ccx = plt.plot(nmp.arange(len(Xlon[:,ii])), 5.*Xlon[:,ii], 'k', linewidth=0.5)


del XFLD
print 'Done!'




cm = plt.imshow(pmsk, cmap = pal_lsm, norm = norm_lsm, interpolation='none')

del pmsk


plt.axis([ 0, Ni, 0, Nj])

if l_scientific_mode:
    plt.xlabel('i-points', **cfont_axis)
    plt.ylabel('j-points', **cfont_axis)



#plt.title('NEMO: '+cfield+', coupled '+CNEMO+', '+cday+' '+chour+':00', **cfont_cnfne)


#ax2 = plt.axes([0.3, 0.08, 0.4, 0.025])

if l_show_cb:

    ax2 = plt.axes(vcb)
    
    clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend=cextend)
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
    



if l_show_nm:  ax.annotate(CNEMO, xy=(1, 4), xytext=(x_cnf, y_cnf), **cfont_cnfn)

if l_show_ttl: ax.annotate(CNEMO, xy=(1, 4), xytext=(x_ttl, y_ttl), **cfont_ttl)




    

#plt.savefig(cfig, dpi=dpi, orientation='portrait', facecolor='b', transparent=True)
plt.savefig(cfig, dpi=dpi, orientation='portrait', transparent=True)
print cfig+' created!\n'
plt.close(1)


del cm, fig, ax



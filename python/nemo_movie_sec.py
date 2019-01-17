#!/usr/bin/env python

#       B a r a K u d a
#
#  Prepare 2D maps of vertical section (monthly) that will later become a movie!
#  NEMO output and observations needed
#
#    L. Brodeau, January 2019

import sys
from os import path, getcwd
from string import replace
import numpy as nmp

from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.image as image
import matplotlib.cbook as cbook

import warnings
warnings.filterwarnings("ignore")

from calendar import isleap
import datetime

from re import split

import barakuda_colmap as bcm

import barakuda_plot as bp
import barakuda_tool as bt
import barakuda_ncio as bnc



cwd = getcwd()


l_smooth = False

vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

fig_type='png'
dpi = 110
color_top = 'white'
#color_top = 'k'

jt0 = 0


i2=0
j2=0
l_get_name_of_run = False
l_show_lsm = True
l_show_cb = False
l_annotate_name = False
l_show_clock = True


cdir_logos = cwd+'/logos'
l_add_logo_on = False
cf_logo_on  = cdir_logos+'/ocean-next_trans_white_281x191.png'
l_add_logo_ige = False
cf_logo_ige = cdir_logos+'/IGE_blanc_notext.png'
l_add_logo_prace = False
cf_logo_prace = cdir_logos+'/PRACE_blanc.png'

rof_log = 150.



l_save_nc = False ; # save the field we built in a netcdf file !!!

l_apply_lap   = False
l_apply_hgrad = False

cb_extend = 'both' ;#colorbar extrema

narg = len(sys.argv)
if narg != 7: print 'Usage: '+sys.argv[0]+' <NEMOCONF> <SEC> <WHAT (SST, SSH, MLD, SSU, SSV, LAP_SSH, GRAD_SST)> <file> <LSM_file> <YYYYMMDD (start)>'; sys.exit(0)
CNEMO  = sys.argv[1]
CSEC   = sys.argv[2]
CWHAT  = sys.argv[3]
cf_in = sys.argv[4]
cf_lsm=sys.argv[5]
cf_clock0=sys.argv[6]




x_logo  = 50 ; y_logo  = 50

cv_msk = 'tmask'

rfactor = 1.0

if   CWHAT == 'U':
    cv_in = 'vozocrtx'
    tmin=-0.4 ; tmax=-tmin ; df=0.1 ; cpal_fld='RdBu_r' ; cb_jump=1
    cunit=r'Zonal component of current speed [m/s]'
    cv_msk='umask' ; l_show_cb=True

elif CWHAT == 'V':
    cv_in = 'vomecrty'
    tmin=-0.4 ; tmax=-tmin ; df=0.1 ; cpal_fld='RdBu_r' ; cb_jump=1
    cunit=r'Meridional component of current speed [m/s]'
    cv_msk='vmask' ; l_show_cb=True

elif CWHAT == 'W':
    cv_in = 'vovecrtz'
    rfactor = 24.*3600. ; # => m/day
    tmin=-300. ; tmax=-tmin ; df=5. ; cpal_fld='RdBu_r' ; cb_jump=1
    cunit=r'Vertical current speed [m/day]'
    cv_msk='tmask' ; l_show_cb=True

elif CWHAT == 'T':
    cv_in = 'sosstsst' ; #in ['sosstsst','tos']:    
    tmin=-2 ;  tmax=30.   ;  df = 1. ; cpal_fld = 'ncview_nrl' ;     cb_jump = 2
    #tmin=0. ;  tmax=32.   ;  df = 2. ; cpal_fld = 'viridis'
    #tmin=4. ;  tmax=20.   ;  df = 1. ; cpal_fld = 'PuBu'
    cunit = r'T ($^{\circ}$C)'
    l_show_cb = True

elif CWHAT == 'GRAD_T':
    cv_in = 'sosstsst'
    l_apply_hgrad = True
    l_smooth = True ; nb_smooth  = 5
    tmin=0. ;  tmax=0.001 ;  df = 0.0001 ; cpal_fld = 'ncview_hotres' ; cb_jump = 1
    #tmin=0. ;  tmax=32.   ;  df = 2. ; cpal_fld = 'viridis'
    #tmin=4. ;  tmax=20.   ;  df = 1. ; cpal_fld = 'PuBu'
    cunit = r'$\left|\vec{\nabla}T\right|$ (K/m)'
    l_show_cb = True
    
else:
    print 'ERROR: we do not know variable "'+str(cv_in)+'" !'
    sys.exit(0)

    
if CNEMO == 'eNATL60':

    # Defaults:
    Ni0 = 8354
    Nj0 = 4729
    cdt = '1h'
    l_get_name_of_run = True

    # Boxes:
    if   CSEC == 'Azores':
        i1=4175 ; j1=1000 ; i2=i1 ; j2=3000 ; k_stop=294 ; x_min = 22.5 ; x_max = 49.0 ; dx=2.
        size_img_px=nmp.array([1920.,800.]) ; rfact_zoom=1. ; vcb=[0.4, 0.15, 0.5, 0.02]  ; font_rat=1.6
        l_show_clock=True ; x_clock=1600 ; y_clock=300
        l_save_nc=False

    else:
        print ' ERROR: unknow section "'+CSEC+'" for config "'+CNEMO+'" !!!'
        sys.exit(0)

elif CNEMO == 'eNATL4':
    # Defaults:
    Ni0 = 559
    Nj0 = 313
    l_add_logo_on = False
    x_clock = 1600 ; y_clock = 200 ; x_logo = 2200 ; y_logo  = 50
    cdt = '1h'

    print 'FIX ME!!!'; sys.exit(0)

    # Boxes:
    if  CSEC == 'ALL':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=3. ; vcb=[0.59, 0.1, 0.39, 0.018]  ; font_rat=0.7*rfact_zoom
        l_show_cb = True

    
else:
    print 'ERROR: we do not know NEMO config "'+str(CNEMO)+'" !'
    sys.exit(0)



CRUN = ''
if l_get_name_of_run:
    # Name of RUN:
    vv = split('-|_', path.basename(cf_in))
    if vv[0] != CNEMO:
        print 'ERROR: your file name is not consistent with "'+CNEMO+'" !!! ('+vv[0]+')' ; sys.exit(0)
    CRUN = vv[1]
    print '\n Run is called: "'+CRUN+'" !\n'

    


    

print '\n================================================================'
print '\n rfact_zoom = ', rfact_zoom
print ' font_rat = ', font_rat, '\n'


print ' i1,i2,j1,j2,k_stop =>', i1,i2,j1,j2, k_stop

l_zonal = False
l_merid = False

if i1==i2:
    l_merid = True
    nx_res = j2-j1
    i2=i2+1

if j1==j2:
    l_zonal = True
    nx_res = i2-i1
    j2=j2+1

if l_merid and l_zonal: print 'ERROR: cannot be zonal AND meridional!!!!' ; sys.exit(0)

ny_res = k_stop

print ' *** nx_res, ny_res =', nx_res, ny_res

yx_ratio = float(nx_res+1)/float(ny_res+1)

rnxr = rfact_zoom*nx_res ; # widt image (in pixels)
rnyr = rfact_zoom*ny_res ; # height image (in pixels)

# Target resolution for figure:
size_figure = size_img_px/float(dpi)

#rh_fig = round(rnyr/float(dpi),3) ; # width of figure
#rw_fig = round(rh_fig*yx_ratio      ,3) ; # height of figure
#rh_img = rh_fig*float(dpi)
#rw_img = rw_fig*float(dpi)
#while rw_img < round(rnxr,0):
#    rw_fig = rw_fig + 0.01/float(dpi)
#    rw_img = rw_fig*float(dpi)
#while rh_img < round(rnyr,0):
#    rh_fig = rh_fig + 0.01/float(dpi)
#    rh_img = rh_fig*float(dpi)
#    print ' *** size figure =>', rw_fig, rh_fig, '\n'
#    print ' *** Forecasted dimension of image =>', rw_img, rh_img
#
print '\n================================================================\n\n\n'

print ' size_figure =>', size_figure

cyr0=cf_clock0[0:4]
cmn0=cf_clock0[4:6]
cdd0=cf_clock0[6:8]


l_3d_field = False


bt.chck4f(cf_lsm)

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





bt.chck4f(cf_lsm)
id_lsm = Dataset(cf_lsm)
#
nb_dim = len(id_lsm.variables[cv_msk].dimensions)
print ' The mesh_mask has '+str(nb_dim)+' dimmensions!'
print ' *** Reading '+cv_msk+' !'
if nb_dim==4: XMSK = nmp.squeeze( id_lsm.variables[cv_msk][0,0:k_stop,j1:j2,i1:i2] )
if nb_dim==3: XMSK = nmp.squeeze( id_lsm.variables[cv_msk][  0:k_stop,j1:j2,i1:i2] )
if nb_dim==2: print 'ERROR: cannot be a 2D field!!!'; sys.exit(0)
print nmp.shape(XMSK)
(nj,ni) = nmp.shape(XMSK)

Vdepth = nmp.squeeze( id_lsm.variables['gdept_1d'][0,0:k_stop] )
Vlon   = nmp.squeeze( id_lsm.variables['glamt'][0,j1:j2,i1:i2] )
Vlat   = nmp.squeeze( id_lsm.variables['gphit'][0,j1:j2,i1:i2] )


print ' Shape Vlon:', nmp.shape(Vlon)
print ' Shape Vlat:', nmp.shape(Vlat)

Vx = nmp.zeros(ni)
if l_merid: Vx[:] = Vlat[:]
if l_zonal: Vx[:] = Vlon[:]

print ' Shape Vx:', nmp.shape(Vx)
print ' Shape Vdepth:', nmp.shape(Vdepth),'\n'

print " *** Vx = ", Vx[:],'\n'
print " *** Vdepth = ", Vdepth[:],'\n'
#sys.exit(0);#lolo


if l_apply_lap:
    print ' *** Reading e1t and e2t !'
    XE1T2 = id_lsm.variables['e1t'][0,j1:j2,i1:i2]
    XE2T2 = id_lsm.variables['e2t'][0,j1:j2,i1:i2]
    XE1T2 = XE1T2*XE1T2
    XE2T2 = XE2T2*XE2T2
if l_apply_hgrad:
    print ' *** Reading e1u and e2v !'
    XE1U = id_lsm.variables['e1u'][0,j1:j2,i1:i2]
    XE2V = id_lsm.variables['e2v'][0,j1:j2,i1:i2]
    print ' *** Reading umask and vmask !'
    if nb_dim==4:
        UMSK = id_lsm.variables['umask'][0,0:k_stop,j1:j2,i1:i2]
        VMSK = id_lsm.variables['vmask'][0,0:k_stop,j1:j2,i1:i2]
    if nb_dim==3:
        UMSK = id_lsm.variables['umask'][0:k_stop,j1:j2,i1:i2]
        VMSK = id_lsm.variables['vmask'][0:k_stop,j1:j2,i1:i2]

#
id_lsm.close()

print 'Shape Arrays => ni,nj =', ni,nj

print 'Done!\n'


idx_land = nmp.where( XMSK < 0.01)



params = { 'font.family':'Helvetica Neue',
           'font.weight':    'normal',
           'font.size':       int(9.*font_rat),
           'legend.fontsize': int(9.*font_rat),
           'xtick.labelsize': int(9.*font_rat),
           'ytick.labelsize': int(9.*font_rat),
           'axes.labelsize':  int(9.*font_rat) }
mpl.rcParams.update(params)
cfont_clb  =  { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(9.*font_rat), 'color':color_top}
cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(9.*font_rat), 'color':color_top }
cfont_mail =  { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*font_rat), 'color':'0.8'}
cfont_titl =  { 'fontname':'Helvetica Neue', 'fontweight':'light', 'fontsize':int(30.*font_rat), 'color':color_top }


# Colormaps for fields:
pal_fld = bcm.chose_colmap(cpal_fld)
norm_fld = colors.Normalize(vmin=tmin, vmax=tmax, clip=False)


if l_show_lsm:
    pal_lsm = bcm.chose_colmap('land_dark')
    norm_lsm = colors.Normalize(vmin=0., vmax=1., clip=False)


if cdt == '3h':
    dt = 3
elif cdt == '1h':
    dt = 1
else:
    print 'ERROR: unknown dt!'




ntpd = 24/dt

vm = vmn
if isleap(int(cyr0)): vm = vml
#print ' year is ', vm, nmp.sum(vm)


print '\n\n Opening file '+cf_in+' !'
id_fld = Dataset(cf_in)


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
    #cdc= '%3.3i'%(jdc)
    cd = '%3.3i'%(jd)
    cm = '%2.2i'%(jm)

    #print '\n\n *** jt, ch, cd, cm =>', jt, ch, cd, cm


    ct = str(datetime.datetime.strptime(cyr0+'-'+cm+'-'+cd+' '+ch, '%Y-%m-%j %H'))
    ct=ct[:5]+cm+ct[7:] #lolo bug !!! need to do that to get the month and not "01"
    print ' ct = ', ct
    cday  = ct[:10]   ; print ' *** cday  :', cday
    chour = ct[11:13] ; print ' *** chour :', chour



    cfig = 'figs/'+cv_in+'_'+CNEMO+'-'+CRUN+'_SECTION-'+CSEC+'_'+cday+'_'+chour+'_'+cpal_fld+'.'+fig_type

    ###### FIGURE ##############

    fig = plt.figure(num = 1, figsize=(size_figure[0], size_figure[1]), dpi=None) ###, facecolor='0.5', edgecolor='k')

    ax = plt.axes([0.05, 0.042, 0.94, 0.93], axisbg = '0.35') ; # '0.35' inside ploting region

    vc_fld = nmp.arange(tmin, tmax + df, df)


    print "Reading record #"+str(jt)+" of "+cv_in+" in "+cf_in
    if l_notime:
        Xplot  = nmp.squeeze(  rfactor*id_fld.variables[cv_in][0:k_stop,j1:j2,i1:i2] )
    else:
        Xplot  = nmp.squeeze(  rfactor*id_fld.variables[cv_in][jt,0:k_stop,j1:j2,i1:i2] ) ; # t, y, x        
    print "Done!\n"

    #print ' *** shape of Xplot => ', nmp.shape(Xplot)

    #if l_apply_lap:
    #    print ' *** Computing Laplacian of "'+cv_in+'"!'
    #    lx = nmp.zeros((nj,ni))
    #    ly = nmp.zeros((nj,ni))
    #    lx[:,1:ni-1] = 1.E9*(Xplot[:,2:ni] -2.*Xplot[:,1:ni-1] + Xplot[:,0:ni-2])/XE1T2[:,1:ni-1]
    #    ly[1:nj-1,:] = 1.E9*(Xplot[2:nj,:] -2.*Xplot[1:nj-1,:] + Xplot[0:nj-2,:])/XE2T2[1:nj-1,:]
    #    Xplot[:,:] = lx[:,:] + ly[:,:]
    #    del lx, ly

    #if l_apply_hgrad:
    #    print ' *** Computing gradient of "'+cv_in+'"!'
    #    lx = nmp.zeros((nj,ni))
    #    ly = nmp.zeros((nj,ni))
    #     if l_smooth: bt.smoother(Xplot, XMSK, nb_smooth=nb_smooth)
    #    # Zonal gradient on T-points:
    #    lx[:,1:ni-1] = (Xplot[:,2:ni] - Xplot[:,0:ni-2]) / (XE1U[:,1:ni-1] + XE1U[:,0:ni-2]) * UMSK[:,1:ni-1] * UMSK[:,0:ni-2]
    #    lx[:,:] = XMSK[:,:]*lx[:,:]
    #    #bnc.dump_2d_field('dsst_dx_gridT.nc', lx, xlon=Xlon, xlat=Xlat, name='dsst_dx')
    #    # Meridional gradient on T-points:
    #    ly[1:nj-1,:] = (Xplot[2:nj,:] - Xplot[0:nj-2,:]) / (XE2V[1:nj-1,:] + XE2V[0:nj-2,:]) * VMSK[1:nj-1,:] * VMSK[0:nj-2,:]
    #    ly[:,:] = XMSK[:,:]*ly[:,:]
    #    #bnc.dump_2d_field('dsst_dy_gridT.nc', ly, xlon=Xlon, xlat=Xlat, name='dsst_dy')
    #    Xplot[:,:] = 0.0
    #    # Modulus of vector gradient:        
    #    Xplot[:,:] = nmp.sqrt(  lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] )
    #    #bnc.dump_2d_field('mod_grad_sst.nc', Xplot, xlon=Xlon, xlat=Xlat, name='dsst')
    #    del lx, ly


    print ''
    print '  *** dimension of array => ', ni, nj, nmp.shape(Xplot)

    if l_save_nc:
        cf_out = 'nc/'+CWHAT+'_NEMO_'+CNEMO+'_'+CSEC+'_'+cday+'_'+chour+'_'+cpal_fld+'.nc'
        print ' Saving in '+cf_out
        bnc.dump_2d_field(cf_out, Xplot, xlon=Vx, xlat=Vdepth, name=CWHAT)
        bnc.dump_2d_field(cf_out, Xplot, name=CWHAT)
        print ''


    

    print "Ploting"

    plt.axis([ x_min, x_max, nmp.max(Vdepth),  nmp.min(Vdepth) ])

    #if l_show_lsm: Xplot[idx_land] = nmp.nan
    
    #cf = plt.imshow(Xplot[:,:], cmap = pal_fld, norm = norm_fld, interpolation='none')

    cf = plt.pcolormesh(Vx[:], Vdepth[:], Xplot[:,:], cmap=pal_fld, norm=norm_fld)

    del Xplot

    if l_merid: bp.__nice_latitude_axis__( ax, plt, x_min, x_max, dx, axt='x')
    if l_zonal: bp.__nice_longitude_axis__(ax, plt, x_min, x_max, dx, axt='x')
    bp.__nice_depth_axis__(ax, plt, nmp.min(Vdepth), nmp.max(Vdepth), l_log=False, l_z_inc=False, cunit='[m]') ; ###, cfont=cfont_clb)

    if l_show_lsm:
        clsm = plt.pcolormesh(Vx[:], Vdepth[:], nmp.ma.masked_where(XMSK>0.0001, XMSK), cmap=pal_lsm, norm=norm_lsm) ###, interpolation='none')

    if l_show_cb:
        ax2 = plt.axes(vcb)
        clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend=cb_extend)
        cb_labs = []
        if cb_jump > 1:
            cpt = 0
            for rr in vc_fld:
                if cpt % cb_jump == 0:
                    if df >= 1.: cb_labs.append(str(int(rr)))
                    if df <  1.: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./df)))+1) ))
                else:
                    cb_labs.append(' ')
                cpt = cpt + 1
        else:
            for rr in vc_fld: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./df)))+1) ))

        clb.ax.set_xticklabels(cb_labs, **cfont_clb)
        clb.set_label(cunit, **cfont_clb)
        clb.ax.yaxis.set_tick_params(color=color_top) ; # set colorbar tick color
        clb.outline.set_edgecolor(color_top) ; # set colorbar edgecolor
        clb.ax.tick_params(which = 'minor', length = 2, color = color_top )
        clb.ax.tick_params(which = 'major', length = 4, color = color_top )


    if l_show_clock:
        xl = float(x_clock)/rfact_zoom
        yl = float(y_clock)/rfact_zoom
        ax.annotate('Date: '+cday+' '+chour+':00', xy=(1, 4), xytext=(xl,yl), **cfont_clock)

    #ax.annotate('laurent.brodeau@ocean-next.fr', xy=(1, 4), xytext=(xl+150, 20), **cfont_mail)

    if l_annotate_name:
        xl = rnxr/20./rfact_zoom
        yl = rnyr/1.33/rfact_zoom
        ax.annotate(CNEMO, xy=(1, 4), xytext=(xl, yl), **cfont_titl)

    if l_add_logo_on:
        datafile = cbook.get_sample_data(cf_logo_on, asfileobj=False)
        im = image.imread(datafile)
        #im[:, :, -1] = 0.5  # set the alpha channel
        fig.figimage(im, x_logo, y_logo, zorder=9)
        del datafile, im
        #
    if l_add_logo_ige:
        datafile = cbook.get_sample_data(cf_logo_ige, asfileobj=False)
        im = image.imread(datafile)
        fig.figimage(im, x_logo+144, y_logo-150., zorder=9)
        del datafile, im
        #
    if l_add_logo_prace:
        datafile = cbook.get_sample_data(cf_logo_prace, asfileobj=False)
        im = image.imread(datafile)
        fig.figimage(im, x_logo-77, y_logo-140., zorder=9)
        del datafile, im

    plt.savefig(cfig, dpi=dpi, orientation='portrait', facecolor='w') ; # white outside ploting region
    print cfig+' created!\n'
    plt.close(1)

    #if l_show_lsm: del clsm
    del cf, fig, ax
    if l_show_cb: del clb

id_fld.close()






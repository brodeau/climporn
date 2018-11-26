#!/usr/bin/env python

#       B a r a K u d a
#
#  Prepare 2D maps (monthly) that will later become a GIF animation!
#  NEMO output and observations needed
#
#    L. Brodeau, May 2018

import sys
from os import path
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

import barakuda_tool as bt
import barakuda_ncio as bnc


vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

fig_type='png'
dpi = 110
color_top = 'white'
#color_top = 'k'

cv_out = 'unknown'

jt0 = 0


i2=0
j2=0
l_get_name_of_run = False
l_show_lsm = True
l_do_ice  = True
l_show_cb = True
l_log_field = False
l_pow_field = False
l_annotate_name = True
l_show_clock = True

l_add_logo = True
cf_logo = '/home/brodeau/Dropbox/OceanNext/Graphic_Identity/0LOGO/logo_trans_white_H14_20180917.png'

l_save_nc = False ; # save the field we built in a netcdf file !!!

romega = 2.*nmp.pi/86400.0



narg = len(sys.argv)
if narg < 10: print 'Usage: '+sys.argv[0]+' <NEMOCONF> <BOX> <WHAT (CURL,CURLOF,CSPEED)><fileX> <varX> <fileY> <varY> <LSM_file> <YYYYMMDD (start)>'; sys.exit(0)
CNEMO  = sys.argv[1]
CBOX   = sys.argv[2]
CWHAT  = sys.argv[3]
cfx_in = sys.argv[4] ; cvx_in = sys.argv[5]
cfy_in = sys.argv[6] ; cvy_in = sys.argv[7]
cf_lsm = sys.argv[8] ; cf_clock0=sys.argv[9]

x_logo  = 50 ; y_logo  = 50



l_do_crl=False; l_do_cof=False; l_do_cspd=False ; 

if   CWHAT == 'CURL':
    l_do_crl = True  ; # do curl (relative-vorticity) !!!
elif CWHAT == 'CURLOF':
    l_do_cof = True  ; # do curl/f
elif CWHAT == 'CSPEED':
    l_do_cspd = True  ; # do current speed
else:
    print ' ERROR: unknow diagnostic "'+CWHAT+'" !!!'
    sys.exit(0)
    
cv_out = CWHAT


    
if cvx_in=='sozocrtx' and cvy_in=='somecrty' and l_do_cof:
    cpal_fld = 'on2' ; tmin=-0.6 ;  tmax=-tmin ;  df = 0.05
    #cpal_fld = 'ncview_bw' ; tmin=-0.4 ;  tmax=0.4 ;  df = 0.05
    cunit = ''
    cb_jump = 1

elif cvx_in=='sozocrtx' and cvy_in=='somecrty' and l_do_crl:
    cpal_fld = 'on2' ; tmin=-0.035 ;  tmax=0.035 ;  df = 0.05
    #cpal_fld = 'ncview_bw' ; tmin=-0.03 ;  tmax=0.03 ;  df = 0.05
    cpal_fld = 'bone' ; tmin=-0.025 ;  tmax=0.025 ;  df = 0.05
    cunit = '';  cb_jump = 1
    l_show_clock = False
    l_add_logo   = False
    l_annotate_name = False

elif cvx_in=='sozocrtx' and cvy_in=='somecrty' and l_do_cspd:
    # Current speed !
    l_save_nc=False
    cpal_fld = 'on2' ; tmin=0. ;  tmax=1.4 ;  df = 0.1
    cunit = 'm/s' ; cb_jump = 1
    
elif cvx_in=='vozocrtx' and cvy_in=='vomecrty' and l_do_cof:
    l_3d_field = True
    #cpal_fld = 'on2' ; tmin=-1. ;  tmax=1. ;  df = 0.05
    cpal_fld = 'ncview_bw' ; tmin=-0.4 ;  tmax=0.4 ;  df = 0.05
    cunit = ''
    cb_jump = 1

elif cvx_in=='vozocrtx' and cvy_in=='vomecrty' and l_do_crl:
    l_3d_field = True
    #cpal_fld = 'on2' ; tmin=-0.025 ;  tmax=0.025 ;  df = 0.05
    #cpal_fld = 'ncview_bw' ; tmin=-0.025 ;  tmax=0.025 ;  df = 0.05
    #cpal_fld = 'gray' ; tmin=-0.025 ;  tmax=0.025 ;  df = 0.05
    cpal_fld = 'bone' ; tmin=-0.025 ;  tmax=0.025 ;  df = 0.05
    cunit = '';  cb_jump = 1
    l_show_clock = False
    l_add_logo   = False

else:
    print 'ERROR: we do not know cvx_in and cvy_in! ("'+cvx_in+'", "'+cvy_in+'")'
    sys.exit(0)

    
if CNEMO == 'eNATL60':

    # Defaults:
    Ni0 = 8354
    Nj0 = 4729
    l_do_ice  = False
    l_show_clock = True
    x_clock = 1600 ; y_clock = 200 ; x_logo = 2200 ; y_logo  = 50
    cdt = '1h'
    l_get_name_of_run = True

    # Boxes:
    if   CBOX == 'ALL':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1440./float(Nj0) ; vcb=[0.61, 0.1, 0.36, 0.018]  ; font_rat=8.*rfact_zoom
        x_clock = 1600 ; y_clock = 200 ; x_logo=2200 ; y_logo=1200

    elif   CBOX == 'ALLFR':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.39, 0.018]  ; font_rat=8.*rfact_zoom
        x_clock = 4000 ; y_clock = 200 ; x_logo = 6000 ; y_logo  = 50; l_show_clock=False ; l_annotate_name=False; l_add_logo=False

    elif   CBOX == 'SALL':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1080./float(Nj0) ; vcb=[0.59, 0.1, 0.39, 0.018]  ; font_rat=8.*rfact_zoom
        x_clock = 1600 ; y_clock = 200 ; x_logo = 2200 ; y_logo  = 50

    elif CBOX == 'FullMed':
        i1=5520; j1=1525; i2=i1+2560 ; j2=j1+1440 ; rfact_zoom=1. ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        l_annotate_name=False

    elif CBOX == 'Med+BS':
        i1=5400; j1=1530; i2=Ni0 ; j2=3310 ; rfact_zoom=1440./float(j2-j1)   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        l_annotate_name=False

    elif CBOX == 'BlackSea':
        i1=Ni0-1920; j1=3330-1080; i2=Ni0 ; j2=3330 ; rfact_zoom=1.   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
        x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16 ; l_annotate_name=False

    elif CBOX == 'Brittany':
        i1=5400; j1=2850; i2=5700 ; j2=3100 ; rfact_zoom=4.   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=1.*rfact_zoom
        x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16 ; l_annotate_name=False
        if cfield == 'SST': tmin=7. ;  tmax=13.   ;  df = 1. ; cpal_fld = 'ncview_nrl'

    elif CBOX == 'Portrait':
        i1=2760; j1=1000; i2=4870; j2=4000 ; rfact_zoom=1.     ; vcb=[0.59, 0.1, 0.38, 0.018]  ; font_rat=1.*rfact_zoom
        l_annotate_name=False; l_show_clock=False

    elif CBOX == 'EATL':
        i1=3100; j1=2290; i2=i1+1800; j2=j1+1080 ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.38, 0.018] ; font_rat = 2.           ; l_annotate_name=False
        x_clock = 1420 ; y_clock = 1030 ; x_logo = 1500 ; y_logo = 16

    elif CBOX == 'EATL2':
        i1=3040; j1=2140; i2=i1+1920; j2=j1+1376 ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.38, 0.018] ; font_rat = 2.
        l_annotate_name=False; l_show_clock=False ; l_add_logo = False
        
    elif CBOX == 'GrlIcl':
        i1=3370; j1=3941; i2=5062; j2=Nj0 ; rfact_zoom=1. ; vcb=[0.3, 0.1, 0.38, 0.018] ; font_rat = 2. ; l_annotate_name=False
        x_clock = 1350 ; y_clock = 750 ; x_logo = 1400 ; y_logo = 16 ; l_save_nc=True

    elif CBOX == 'Band':
        i1=5100-1920; j1=2200; i2=5100; j2=j1+1080 ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.38, 0.018] ; font_rat = 2.           ; l_annotate_name=False
        l_show_clock = False ; l_add_logo = False ; #x_clock = 1420 ; y_clock = 1030 ; x_logo = 1500 ; y_logo = 16

    elif CBOX == 'Balear':
        i1=5750; j1=1880; i2=6470; j2=2600 ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.38, 0.018] ; font_rat = 2. ; l_annotate_name=False
        x_clock = 1420 ; y_clock = 1030 ; x_logo = 1500 ; y_logo = 16

    else:
        print ' ERROR: unknow box "'+CBOX+'" for config "'+CNEMO+'" !!!'
        sys.exit(0)



elif CNEMO == 'NATL60':
    Ni0 = 5422
    Nj0 = 3454
    #l_pow_field = True ; pow_field = 1.5
    l_do_ice  = False
    l_show_clock = False
    cdt = '1h'
    #CBOX = 'zoom1' ; i1 = 1800 ; j1 = 950 ; i2 = i1+1920 ; j2 = j1+1080 ; rfact_zoom = 1. ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom ; l_show_lsm = False
    #CBOX = 'zoom1' ; i1 = 1800 ; j1 = 950 ; i2 = i1+2560 ; j2 = j1+1440 ; rfact_zoom = 1. ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom
    CBOX = 'ALL' ; i1=0 ; j1=0 ; i2=Ni0 ; j2=Nj0 ; rfact_zoom = 0.4 ; vcb=[0.59, 0.1, 0.38, 0.018] ; font_rat = 4.*rfact_zoom
    x_clock = 350 ; y_clock = 7 ; # where to put the date


elif CNEMO == 'NANUK025':
    l_do_ice = True
    cdt = '3h'; CBOX = 'ALL' ; i1 = 0 ; j1 = 0 ; i2 = 492 ; j2 = 614 ; rfact_zoom = 2. ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat = 8.*rfact_zoom
    x_clock = 350 ; y_clock = 7 ; # where to put the date




elif CNEMO == 'eNATL4':
    # Defaults:
    Ni0 = 559
    Nj0 = 313
    l_do_ice  = False
    l_annotate_name = False
    l_show_clock = False
    l_add_logo = False
    x_clock = 1600 ; y_clock = 200 ; x_logo = 2200 ; y_logo  = 50
    cdt = '1h'

    # Boxes:
    if   CBOX == 'ALL':
        i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=3. ; vcb=[0.59, 0.1, 0.39, 0.018]  ; font_rat=0.7*rfact_zoom
        l_show_cb = True

    
else:
    print 'ERROR: we do not know NEMO config "'+str(CNEMO)+'" !'
    sys.exit(0)


CRUN = ''
if l_get_name_of_run:
    # Name of RUN:
    vv = split('-|_', path.basename(cfx_in))
    if vv[0] != CNEMO:
        print 'ERROR: your file name is not consistent with "'+CNEMO+'" !!! ('+vv[0]+')' ; sys.exit(0)
    CRUN = vv[1]
    print '\n Run is called: "'+CRUN+'" !\n'

    


    

print '\n================================================================'
print '\n rfact_zoom = ', rfact_zoom
print ' font_rat = ', font_rat, '\n'


print ' i1,i2,j1,j2 =>', i1,i2,j1,j2

nx_res = i2-i1
ny_res = j2-j1

print ' *** nx_res, ny_res =', nx_res, ny_res

yx_ratio = float(nx_res+1)/float(ny_res+1)

rnxr = rfact_zoom*nx_res ; # widt image (in pixels)
rnyr = rfact_zoom*ny_res ; # height image (in pixels)

# Target resolution for figure:
rh_fig = round(rnyr/float(dpi),3) ; # width of figure
rw_fig = round(rh_fig*yx_ratio      ,3) ; # height of figure
rh_img = rh_fig*float(dpi)
rw_img = rw_fig*float(dpi)
while rw_img < round(rnxr,0):
    rw_fig = rw_fig + 0.01/float(dpi)
    rw_img = rw_fig*float(dpi)
while rh_img < round(rnyr,0):
    rh_fig = rh_fig + 0.01/float(dpi)
    rh_img = rh_fig*float(dpi)
    print ' *** size figure =>', rw_fig, rh_fig, '\n'
    print ' *** Forecasted dimension of image =>', rw_img, rh_img

print '\n================================================================\n\n\n'



cyr0=cf_clock0[0:4]
cmn0=cf_clock0[4:6]
cdd0=cf_clock0[6:8]


l_3d_field = False



# Ice:
if l_do_ice:
    cv_ice  = 'siconc'
    cf_ice = replace(cfx_in, 'grid_T', 'icemod')
    rmin_ice = 0.5
    cpal_ice = 'ncview_bw'
    vcont_ice = nmp.arange(rmin_ice, 1.05, 0.05)





if l_do_ice: bt.chck4f(cf_ice)

bt.chck4f(cf_lsm)
bt.chck4f(cfx_in)
bt.chck4f(cfy_in)

id_fx = Dataset(cfx_in)
vtime = id_fx.variables['time_counter'][:]
id_fx.close()

Nt = len(vtime)

if l_show_lsm or l_do_crl or l_do_cof or l_do_cspd:
    print "\nReading record metrics in "+cf_lsm
    id_lsm = Dataset(cf_lsm)
    nb_dim = len(id_lsm.variables['tmask'].dimensions)
    print ' The mesh_mask has '+str(nb_dim)+' dimmensions!'
    if l_show_lsm:
        cv_msk = 'tmask'
        if l_do_crl:  cv_msk = 'fmask'
        if l_do_cspd: cv_msk = 'tmask'
        print '\n Reading mask as "'+cv_msk+'" in:'; print cv_msk,'\n'
        if nb_dim==4: XMSK = id_lsm.variables[cv_msk][0,0,j1:j2,i1:i2]
        if nb_dim==3: XMSK = id_lsm.variables[cv_msk][0,j1:j2,i1:i2]
        if nb_dim==2: XMSK = id_lsm.variables[cv_msk][j1:j2,i1:i2]
        (nj,ni) = nmp.shape(XMSK)
    if l_do_crl or l_do_cof:
        # e2v, e1u, e1f, e2f
        e2v = id_lsm.variables['e2v'][0,j1:j2,i1:i2]
        e1u = id_lsm.variables['e1u'][0,j1:j2,i1:i2]
        e1f = id_lsm.variables['e1f'][0,j1:j2,i1:i2]
        e2f = id_lsm.variables['e2f'][0,j1:j2,i1:i2]
    if l_do_cof:
        ## Coriolis Parameter:
        ff  = id_lsm.variables['gphif'][0,j1:j2,i1:i2]
        ff[:,:] = 2.*romega*nmp.sin(ff[:,:]*nmp.pi/180.0)
        (nj,ni) = nmp.shape(XMSK)
    if l_save_nc:
        Xlon = id_lsm.variables['gphit'][j1:j2,i1:i2]
        Xlat = id_lsm.variables['glamt'][j1:j2,i1:i2]
    id_lsm.close()

    print 'Shape Arrays => ni,nj =', ni,nj

    print 'Done!\n'



params = { 'font.family':'Helvetica Neue',
           'font.weight':    'normal',
           'font.size':       int(9.*font_rat),
           'legend.fontsize': int(9.*font_rat),
           'xtick.labelsize': int(9.*font_rat),
           'ytick.labelsize': int(9.*font_rat),
           'axes.labelsize':  int(9.*font_rat) }
mpl.rcParams.update(params)
cfont_clb  =  { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(7.*font_rat), 'color':color_top}
cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(9.*font_rat), 'color':'w' }
cfont_mail =  { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*font_rat), 'color':'0.8'}
cfont_titl =  { 'fontname':'Helvetica Neue', 'fontweight':'light', 'fontsize':int(30.*font_rat), 'color':'w' }


# Colormaps for fields:
pal_fld = bcm.chose_colmap(cpal_fld)
if l_log_field:
    norm_fld = colors.LogNorm(  vmin = tmin, vmax = tmax, clip = False)
if l_pow_field:
    norm_fld = colors.PowerNorm(gamma=pow_field, vmin = tmin, vmax = tmax, clip = False)
else:
    norm_fld = colors.Normalize(vmin = tmin, vmax = tmax, clip = False)


if l_show_lsm:
    pal_lsm = bcm.chose_colmap('land_dark')
    norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

if l_do_ice:
    pal_ice = bcm.chose_colmap(cpal_ice)
    norm_ice = colors.Normalize(vmin = rmin_ice, vmax = 1, clip = False)



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

jd = int(cdd0) - 1
jm = int(cmn0)


Xplot = nmp.zeros((nj,ni))

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



    cfig = 'figs/'+cv_out+'_'+CNEMO+'-'+CRUN+'_'+CBOX+'_'+cday+'_'+chour+'_'+cpal_fld+'.'+fig_type

    ###### FIGURE ##############

    fig = plt.figure(num = 1, figsize=(rw_fig, rh_fig), dpi=None, facecolor='w', edgecolor='0.5')

    ax  = plt.axes([0., 0., 1., 1.], axisbg = '0.5')

    vc_fld = nmp.arange(tmin, tmax + df, df)


    print "Reading record #"+str(jt)+" of "+cvx_in+" in "+cfx_in
    id_fx = Dataset(cfx_in)
    if not l_3d_field:
        XFLD  = id_fx.variables[cvx_in][jt,j1:j2,i1:i2] ; # t, y, x
    else:
        print 'j1:j2 =', j1,j2
        print 'i1:i2 =', i1,i2
        XFLD  = id_fx.variables[cvx_in][jt,0,j1:j2,i1:i2] ; # t, y, x
        id_fx.close()
    print "Done!"

    print "Reading record #"+str(jt)+" of "+cvy_in+" in "+cfy_in
    id_fy = Dataset(cfy_in)
    if not l_3d_field:
        YFLD  = id_fy.variables[cvy_in][jt,j1:j2,i1:i2] ; # t, y, x
    else:
        YFLD  = id_fy.variables[cvy_in][jt,0,j1:j2,i1:i2] ; # t, y, x
        id_fy.close()
    print "Done!"


    if l_do_crl or l_do_cof:

        print '\nComputing curl...'
        lx = nmp.zeros((nj,ni))
        ly = nmp.zeros((nj,ni))

        lx[:,1:ni-1] =   e2v[:,2:ni]*YFLD[:,2:ni] - e2v[:,1:ni-1]*YFLD[:,1:ni-1]
        ly[1:nj-1,:] = - e1u[2:nj,:]*XFLD[2:nj,:] + e1u[1:nj-1,:]*XFLD[1:nj-1,:]

        if l_do_cof: Xplot[:,:] = ( lx[:,:] + ly[:,:] )*XMSK[:,:] / ( e1f[:,:]*e2f[:,:]*ff[:,:] ) # Relative Vorticity...
        if l_do_crl: Xplot[:,:] = ( lx[:,:] + ly[:,:] )*XMSK[:,:] / ( e1f[:,:]*e2f[:,:] ) * 1000. # Curl...

        del lx, ly

        print '... '+cv_out+' computed!\n'



    if l_do_cspd:
        print '\nComputing current speed at T-points ...'
        lx = nmp.zeros((nj,ni))
        ly = nmp.zeros((nj,ni))

        lx[:,2:ni] = 0.5*( XFLD[:,1:ni-1] + XFLD[:,2:ni] )
        ly[2:nj,:] = 0.5*( YFLD[1:nj-1,:] + YFLD[2:nj,:] )

        Xplot[:,:] = nmp.sqrt( lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] ) * XMSK[:,:]
        

    del XFLD,YFLD


    if l_save_nc:
        cf_out = 'nc/'+cv_out+'_NEMO_'+CNEMO+'_'+CBOX+'_'+cday+'_'+chour+'_'+cpal_fld+'.nc'
        print ' Saving in '+cf_out
        bnc.dump_2d_field(cf_out, Xplot, xlon=Xlon, xlat=Xlat, name=cv_out)
        print ''


    

    print "Ploting"

    plt.axis([ 0, ni, 0, nj])

    idx_miss = nmp.where( XMSK < 0.001)
    Xplot[idx_miss] = nmp.nan

    cf = plt.imshow(Xplot[:,:], cmap = pal_fld, norm = norm_fld, interpolation='none')

    # Ice
    if l_do_ice:
        print "Reading record #"+str(jt)+" of "+cv_ice+" in "+cf_ice
        id_ice = Dataset(cf_ice)
        XICE  = id_ice.variables[cv_ice][jt,:,:] ; # t, y, x
        id_ice.close()
        print "Done!"

        #XM[:,:] = XMSK[:,:]
        #bt.drown(XICE, XM, k_ew=2, nb_max_inc=10, nb_smooth=10)
        #ci = plt.contourf(XICE[:,:], vcont_ice, cmap = pal_ice, norm = norm_ice) #

        pice = nmp.ma.masked_where(XICE < rmin_ice, XICE)
        ci = plt.imshow(pice, cmap = pal_ice, norm = norm_ice, interpolation='none') ; del pice, ci
        del XICE

    #LOLO: rm ???
    if l_show_lsm:
        clsm = plt.imshow(nmp.ma.masked_where(XMSK>0.0001, XMSK), cmap = pal_lsm, norm = norm_lsm, interpolation='none')

    if l_show_cb:
        ax2 = plt.axes(vcb)
        clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend='both')
        cb_labs = []
        if cb_jump > 1:
            cpt = 0
            for rr in vc_fld:
                if cpt % cb_jump == 0:
                    if df >= 1.: cb_labs.append(str(int(rr)))
                    if df <  1.: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./df))))))
                else:
                    cb_labs.append(' ')
                cpt = cpt + 1
        else:
            for rr in vc_fld: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./df))))))

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

    if l_add_logo:
        datafile = cbook.get_sample_data(cf_logo, asfileobj=False)
        im = image.imread(datafile)
        #im[:, :, -1] = 0.5  # set the alpha channel
        fig.figimage(im, x_logo, y_logo, zorder=9)


    plt.savefig(cfig, dpi=dpi, orientation='portrait', facecolor='k')
    print cfig+' created!\n'
    plt.close(1)

    if l_show_lsm: del clsm
    del cf, fig, ax
    if l_show_cb: del clb

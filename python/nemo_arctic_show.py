#!/usr/bin/env python3
#
#     CLIMPORN
#
#  Show SST + sea-ice concentration in the Arctic on a polar stereographic projection!
#  NEMO output + mesh_mask needed.
#
#    L. Brodeau, November 2019

import sys
from os import path, getcwd, mkdir
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

import climporn as cp

ldrown = True ; #lolo
l_add_topo_land = False

l_show_ice_colbar = True

l_show_logos = True
f_logo_on    = 'ocean-next_trans_white_120x82.png'
f_logo_ifrmr = 'IFREMER_blanc_small.png'
f_logo_nersc = 'NERSC_white_120p.png'
f_logo_ige   = 'IGE_blanc_small.png'


vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

fig_type='png'
rDPI = 150

color_top = 'white'
color_top_cb = 'white'
if l_add_topo_land:
    color_top = 'k'
    color_top_cb = 'k'

#color_top = 'k'


jt0 = 0



cdt = '6h'
l_get_name_of_run = True

# Ice:
#rmin_ice=0.2
rmin_ice=0.
rmax_ice=1.
#cpal_ice = 'ncview_bw'
#cpal_ice = 'Blues_r'
#cpal_ice = 'bone'
#cpal_ice = 'ice_on'
#cpal_ice = 'ice2_on'
#cpal_ice = 'ice3_on'
cpal_ice = 'ice4_on'

# Continents:
rof_log = 150.
rof_dpt = 0.



vp =  ['nanuk1', 'stere', -60., 40., 122., 57.,    75.,  -12., 10., 'h' ]  # North Pole
#vp =  ['nanuk2', 'stere', -62., 54., 126., 60.,    90.,  -12., 10., 'h' ]  # North Pole


# Normally logos should be found there:
dir_logos = str.replace( path.dirname(path.realpath(__file__)) , 'climporn/python', 'climporn/misc/logos' )
print("\n --- logos found into : "+dir_logos+" !\n")


################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--fin' , required=True,                help='specify the NEMO netCDF file to read from...')
requiredNamed.add_argument('-w', '--what', required=True, default="sst", help='specify the field/diagnostic to plot (ex: sst)')

parser.add_argument('-C', '--conf', default="NANUK025",     help='specify NEMO config (ex: eNATL60)')
parser.add_argument('-N', '--name', default="X",            help='specify experiment name')
parser.add_argument('-m', '--fmm' , default="mesh_mask.nc", help='specify the NEMO mesh_mask file (ex: mesh_mask.nc)')
parser.add_argument('-s', '--sd0' , default="19950101",     help='specify initial date as <YYYYMMDD>')
parser.add_argument('-l', '--lev' , type=int, default=0,    help='specify the level to use if 3D field (default: 0 => 2D)')
parser.add_argument('-I', '--ice' , action='store_true',    help='draw sea-ice concentration layer onto the field')

args = parser.parse_args()

CNEMO = args.conf
CEXP  = args.name
#CBOX  = args.box
CWHAT = args.what
cf_in = args.fin
cf_mm = args.fmm
csd0  = args.sd0
jk    = args.lev
lshow_ice = args.ice
#cf_topo_land = args.zld

print('')
print(' *** CNEMO = ', CNEMO)
print(' *** CEXP  = ', CEXP)
print(' *** CWHAT = ', CWHAT)
print(' *** cf_in = ', cf_in)
print(' *** cf_mm = ', cf_mm)
print(' *** csd0 = ', csd0)
print(' ***   jk  = ', jk)
print(' *** Show ice? =>', lshow_ice)

#l_add_topo_land = False
#if args.zld != None:
#    print(' *** cf_topo_land = ', cf_topo_land)
#    l_add_topo_land = True

l3d = False
if jk > 0:
    l3d=True
else:
    jk=0
###############################################################################################################################################



if l_get_name_of_run:
    if CEXP == 'X':
        # Name of RUN:
        vv = split('-|_', path.basename(cf_in))
        if vv[0] != CNEMO:
            print('ERROR: your file name is not consistent with "'+CNEMO+'" !!! ('+vv[0]+')')
            print(' ==> specify experiment name with "-N"\n'); sys.exit(0)
        CEXP = vv[1]
    print('\n Run is called: "'+CEXP+'" !\n')

#---------------------------------------------------------------

if CNEMO in ['NANUK025', 'CREG025', 'CREG025.L75']:
    jk=0
    j1=0 ; j2=603
    i1=0 ; i2=528
else:
    print('ERRO: unknow conf '+CNEMO)
    ###############################

if CNEMO == 'NANUK025': cxtra_info1 = "OPA - neXtSIM" ; #cxtra_info2 = "   (CREG025)"
if CNEMO[:7] == 'CREG025':  cxtra_info1 = "OPA - LIM3"    ; #cxtra_info2 = "(CREG025)"

#cv_bg = ''
l_only_over_ice = False ; # only plot fields in regions with sea-ice
r_oi_thr = 0.01
rexp_ctrl = 0.

if  CWHAT == 'sst':
    # SST
    cv_in = 'sst'
    cv_if = 'siconc'
    cv_out = CWHAT
    tmin=-2. ;    tmax=20. ; df = 1. ; cb_jump = 2    
    cpal_fld = 'on3'  #cpal_fld = 'ncview_nrl'
    cunit = r'SST [$^{\circ}$C]'

elif CWHAT in [ 'sit', 'sivolu' ] :
    # Sea Ice Thickness
    cv_in = CWHAT
    if CWHAT=='sit'   : cv_if = 'sic'
    if CWHAT=='sivolu': cv_if = 'siconc'
    #
    cv_out = 'sit'
    tmin=0. ; tmax=5.; df = 1 ; cb_jump = 1
    #cpal_fld = 'on3'
    cpal_fld = 'viridis'
    #cpal_fld = 'cividis'
    #cpal_fld = 'cubehelix'
    #cpal_fld = 'cool'
    cunit = 'Sea-Ice thickness [m]'
    l_only_over_ice=True

elif CWHAT == 'damage':
    # Sea Ice Thickness
    cv_in = 'damage'
    cv_if = 'sic'
    #cv_bg = 'sst'
    cv_out = CWHAT
    #tmin=0.9 ; tmax=1.; df = 0.05 ; cb_jump = 1 ; rexp_ctrl = 2. ; #rexp_ctrl = 3.5
    tmin=0.8 ; tmax=1.; df = 0.05 ; cb_jump = 1 ; rexp_ctrl = 3.5
    #cpal_fld = 'ncview_oslo_r' ; l_only_over_ice=True
    cpal_fld = 'ncview_bone_r' ; l_only_over_ice=True
    cunit = 'Damage [-]'

elif CWHAT in [ 'Qns', 'nshfls' ]:
    cv_in = 'nshfls'
    cv_if = 'ice_cover'
    cv_out = 'Qns'
    #tmin=-1250 ;  tmax=250. ; df = 50. ; cb_jump = 5 ; cpal_fld = 'gist_stern_r'
    #tmin=-500. ;  tmax=500. ; df = 100. ; cb_jump = 2 ; cpal_fld = 'ncview_parula'
    tmin=-100. ;  tmax=0. ; df = 25. ; cb_jump = 1 ; cpal_fld = 'ncview_parula_r' ; rexp_ctrl=2.
    #tmin=-100. ;  tmax=0. ; df = 25. ; cb_jump = 1 ; cpal_fld = 'ncview_ssec_r'
    cunit = r'     Non-solar heat flux [$W/m^{2}$]'
    l_only_over_ice = True
    
elif CWHAT == 'Qnet':
    cv_in = 'qt'
    cv_if = 'ice_cover'
    cv_out = CWHAT
    tmin=-1000 ;  tmax=1000. ;  df = 200. ; cb_jump = 2
    #cpal_fld = 'Spectral_r'
    cpal_fld = 'RdBu_r'
    cunit = r'Net heat flux [$W/m^{2}$]'
    
else:
    print('ERROR: we do not know variable "'+str(cv_in)+'" !')
    sys.exit(0)
    

cp.chck4f(cf_mm)


if not path.exists("figs"): mkdir("figs")
cdir_figs = './figs/'+cv_out
if not path.exists(cdir_figs): mkdir(cdir_figs)



l_notime=False
cp.chck4f(cf_in)
id_in = Dataset(cf_in)
list_var = id_in.variables.keys()
if 'time_counter' in list_var:
    vtime = id_in.variables['time_counter'][:]
elif 'time' in list_var:
    vtime = id_in.variables['time'][:]
else:
    l_notime=True
    print('Did not find a time variable! Assuming no time and Nt=1')
    id_in.close()
    Nt = 1
if not l_notime: Nt = len(vtime)
print(' Nt =', Nt)



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

(nj,ni) = nmp.shape(XMSK)  ; print('Shape Arrays => ni,nj =', ni,nj)

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
cfont_clb  =  { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(8.5*fontr), 'color':color_top_cb}
cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(9.*fontr) , 'color':color_top }
#cfont_exp= { 'fontname':'Open Sans'  , 'fontweight':'light', 'fontsize':int(9.*fontr), 'color':color_top }
#cfont_mail  =  { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*fontr), 'color':'0.8'}
cfont_titl1 = { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(18.*fontr), 'color':color_top }
cfont_titl2 = { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(14.*fontr), 'color':color_top }


# for colorbar:
vc_fld = nmp.arange(tmin, tmax + df, df)

l_show_ice_colbar = l_show_ice_colbar and lshow_ice
if l_show_ice_colbar: vc_ice = nmp.arange(rmin_ice, rmax_ice+0.2, 0.2)





# Avec rebords:
#vfig_size = [ 10., 10.8 ]
#vsporg = [0.03, 0.1, 1., 0.8]

# For movie
vfig_size = [ 7.54, 7.2 ]
#vsporg = [0.001, 0.0011, 0.997, 0.999]
#vsporg = [0., 0.0001, 1., 1.001]
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
    print('ERROR: unknown dt!')


if rexp_ctrl > 0.:
    pal_fld = cp.chose_colmap(cpal_fld, exp_ctrl=rexp_ctrl)
else:
    pal_fld = cp.chose_colmap(cpal_fld)

nrm_fld = colors.Normalize(vmin=tmin, vmax=tmax, clip=False)

if lshow_ice:
    pal_ice = cp.chose_colmap(cpal_ice) ; #lolo, exp_ctrl=1.5)
    nrm_ice = colors.Normalize(vmin=rmin_ice, vmax=rmax_ice, clip = False)






    
ntpd = 24/dt

vm = vmn
if isleap(int(cyr0)): vm = vml

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
    print(' ct = ', ct)
    cday  = ct[:10]   ; print(' *** cday  :', cday)
    chour = ct[11:13] ; print(' *** chour :', chour)

    if l3d:
        cfig = cdir_figs+'/'+cv_out+'_'+CNEMO+'-'+CEXP+'_lev'+str(jk)+'_'+cday+'_'+chour+'.'+fig_type
    else:
        cfig = cdir_figs+'/'+cv_out+'_'+CNEMO+'-'+CEXP+'_'+cday+'_'+chour+'.'+fig_type

    # Getting field and sea-ice concentration at time record "jt":
    id_in = Dataset(cf_in)
    XFLD = id_in.variables[cv_in]   [jt,j1:j2,i1:i2]
    if lshow_ice or l_only_over_ice:
        XIFR = id_in.variables[cv_if][jt,j1:j2,i1:i2]
        XIFR = nmp.ma.masked_where(XMSK <= 0.0001, XIFR)
        XIFR_bkp = nmp.zeros(nmp.shape(XIFR)) ; XIFR_bkp[:,:] = XIFR[:,:]
        XIFR_bkp = nmp.ma.masked_where(XMSK <= 0.0001, XIFR_bkp)
    id_in.close()

    cjt = '%4.4i'%(jt)

    #######################################################################################################
    #col_bg = '#3b3b63'
    col_bg = '#041a4d'
    fig = plt.figure(num = 1, figsize=(vfig_size), dpi=None, facecolor=col_bg, edgecolor=col_bg)
    ax  = plt.axes(vsporg, facecolor = col_bg)

    carte = Basemap(llcrnrlon=vp[2], llcrnrlat=vp[3], urcrnrlon=vp[4], urcrnrlat=vp[5], \
                    resolution=vp[9], area_thresh=1000., projection='stere', \
                    lat_0=vp[6], lon_0=vp[7], epsg=None)

    x0,y0 = carte(Xlon,Xlat)

    #f = carte.pcolor(x0, y0, XIFR, cmap = pal_ice, norm = nrm_ice) #, interpolation='none')

    if ldrown:
        print(' Drowning...')
        if lshow_ice: cp.drown(XIFR, XMSK, k_ew=-1, nb_max_inc=10, nb_smooth=10)
        cp.drown(XFLD, XMSK, k_ew=-1, nb_max_inc=10, nb_smooth=10)

    if lshow_ice: XFLD = nmp.ma.masked_where(XIFR >= 0.2, XFLD)
    XFLD = nmp.ma.masked_where(XMSK <= 0.1, XFLD)

    if lshow_ice:
        XIFR = nmp.ma.masked_where(XIFR  < r_oi_thr, XIFR)
    
    if l_only_over_ice:
        XFLD = nmp.ma.masked_where(XIFR <  r_oi_thr, XFLD)


    if l_add_topo_land:
        #carte.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='ESRI_Imagery_World_2D', xpixels=400, ypixels=None, dpi=96, verbose=True)
        #carte.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 2000, verbose= True)
        #carte.etopo()
        carte.shadedrelief()
        carte.drawlsmask(ocean_color='0.5',land_color=(255,255,255,1), alpha=1)
        #carte.drawlsmask(land_color=(255,255,255,1))

    ft = carte.pcolormesh(x0, y0, XFLD, cmap = pal_fld, norm = nrm_fld )
    #ft = carte.pcolor(x0, y0, XFLD, cmap = pal_fld, norm = nrm_fld )
    if lshow_ice:
        fi = carte.pcolormesh(x0, y0, XIFR, cmap = pal_ice, norm = nrm_ice )
        #fc = carte.contour(   x0, y0, XIFR_bkp, [r_oi_thr], colors='k', linewidths=1. )

    #if l_only_over_ice:
    #    # A contour to make limit smoother...
    #    #fc = carte.contour(x0, y0, XIFR, [r_oi_thr/2.], colors='red', linewidths=2.1 )
    #    fc = carte.contour(x0, y0, XIFR_bkp, [0.001], colors='k', linewidths=3. )
    
    carte.drawcoastlines(linewidth=0.5)

    if not l_add_topo_land: carte.fillcontinents(color='grey') #, alpha=0)

    #carte.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    #carte.drawmapboundary()

    carte.drawmeridians(nmp.arange(-180,180,20), labels=[0,0,0,1], linewidth=0.3)
    carte.drawparallels(nmp.arange( -90, 90,10), labels=[1,0,0,0], linewidth=0.3)

    # ----------- Color bar for field -----------
    ax2 = plt.axes([0.64, 0.965, 0.344, 0.018])
    clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=nrm_fld, orientation='horizontal', extend='both')
    cb_labs = []
    cpt = 0
    for rr in vc_fld:
        if cpt % cb_jump == 0 or ( (tmin == -tmax) and (int(rr) == 0 ) ):
            if df >= 1.: cb_labs.append(str(int(rr)))
            if df <  1.: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./df)))+1) ))
        else:
            cb_labs.append(' ')
        cpt = cpt + 1
    clb.ax.set_xticklabels(cb_labs, **cfont_clb)
    clb.set_label(cunit, **cfont_clb)
    clb.ax.yaxis.set_tick_params(color=color_top_cb) ; # set colorbar tick color                                                                      
    clb.outline.set_edgecolor(color_top_cb) ; # set colorbar edgecolor                                                                                
    clb.ax.tick_params(which = 'minor', length = 2, color = color_top_cb )
    clb.ax.tick_params(which = 'major', length = 4, color = color_top_cb )
    del clb

    # ----------- Color bar for ICE fraction -----------
    if l_show_ice_colbar:
        ax3 = plt.axes([0.02, 0.965, 0.29, 0.018])
        clb = mpl.colorbar.ColorbarBase(ax3, ticks=vc_ice, cmap=pal_ice, norm=nrm_ice, orientation='horizontal') #, extend='min')
        cb_labs = []
        for rr in vc_ice:
            if rr == 1.:
                cb_labs.append(str(int(rr)))
            else:
                cb_labs.append(str(round(rr,1)))
        clb.ax.set_xticklabels(cb_labs, **cfont_clb)
        clb.set_label('Ice fraction', **cfont_clb)
        clb.ax.yaxis.set_tick_params(color=color_top_cb) ; # set colorbar tick color                                                                      
        clb.outline.set_edgecolor(color_top_cb) ; # set colorbar edgecolor                                                                                
        #clb.ax.tick_params(which = 'minor', length = 2, color = color_top_cb )
        clb.ax.tick_params(which = 'major', length = 4, color = color_top_cb )




    
    #ax.annotate('Date: '+cday+' '+chour+':00', xy=(0.785, 0.015), xycoords='figure fraction', **cfont_clock)
    ax.annotate('Date: '+cday+' '+chour+':00', xy=(0.76, 0.88), xycoords='figure fraction', **cfont_clock)
    
    ry0 = 0.78
    #ry0 = 0.9
    ax.annotate(cxtra_info1, xy=(0.02, ry0+0.05), xycoords='figure fraction', **cfont_titl1)
    #ax.annotate(cxtra_info2, xy=(0.05, ry0 ),     xycoords='figure fraction', **cfont_titl2)

    #
    if l_show_logos:

        datafile = cbook.get_sample_data(dir_logos+'/'+f_logo_ige, asfileobj=False)
        im = image.imread(datafile)
        fig.figimage(im, 1010, 855, zorder=9)
        del datafile, im

        datafile = cbook.get_sample_data(dir_logos+'/'+f_logo_nersc, asfileobj=False)
        im = image.imread(datafile)
        fig.figimage(im, 990, 775, zorder=9)
        del datafile, im
    
        datafile = cbook.get_sample_data(dir_logos+'/'+f_logo_ifrmr, asfileobj=False)
        im = image.imread(datafile)
        if vp[0] == 'nanuk1':
            fig.figimage(im, 999, 722, zorder=9)
        elif vp[0] == 'nanuk2':
            fig.figimage(im, 1010, 722, zorder=9)
        del datafile, im
    
        
        datafile = cbook.get_sample_data(dir_logos+'/'+f_logo_on, asfileobj=False)
        im = image.imread(datafile)
        if vp[0] == 'nanuk1':
            fig.figimage(im, 990, 9, zorder=9)
        elif  vp[0] == 'nanuk2':
            fig.figimage(im, 990, 620, zorder=9)
        del datafile, im
    
    print(' Saving figure: '+cfig+'\n\n')

    plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)
    plt.close(1)

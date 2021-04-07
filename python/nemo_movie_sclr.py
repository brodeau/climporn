#!/usr/bin/env python3
#
#     CLIMPORN
#
#  Prepare 2D maps (monthly) that will later become a movie!
#  NEMO output and observations needed
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

import warnings
warnings.filterwarnings("ignore")

from calendar import isleap
import datetime

from re import split

import climporn as cp



cwd = getcwd()

l_smooth = False

vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

fig_type='png'
rDPI = 110
color_top = 'white'
color_top_cb = 'white'
#color_top = 'k'

cv_out = 'unknown'

jt0 = 0


jk=0
l_show_lsm = True
l_do_ice  = False
l_log_field = False
l_pow_field = False


rof_log = 150.
rof_dpt = 0.

grav = 9.80665 # same as in NEMO 3.6

l_save_nc = False ; # save the field we built in a netcdf file !!!

romega = 7.292115083046062E-005 # same as in NEMO 3.6 / #romega = 2.*nmp.pi/86400.0

l_apply_lap   = False
l_apply_hgrad = False
l_apply_geov  = False

cb_extend = 'both' ;#colorbar extrema

# Normally logos should be found there:
dir_logos = str.replace( path.dirname(path.realpath(__file__)) , 'climporn/python', 'climporn/misc/logos' )
print("\n --- logos found into : "+dir_logos+" !\n")


################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--fin' , required=True,                help='specify the NEMO netCDF file to read from...')
requiredNamed.add_argument('-w', '--what', required=True, help='specify the field/diagnostic to plot (ex: CSPEED,CURLOF,ect.)')

parser.add_argument('-C', '--conf', default="none",           help='specify NEMO config (ex: eNATL60)')
parser.add_argument('-b', '--box' , default="ALL",            help='specify extraction box name (ex: ALL)')
parser.add_argument('-m', '--fmm' , default="mesh_mask.nc",   help='specify the NEMO mesh_mask file (ex: mesh_mask.nc)')
parser.add_argument('-s', '--sd0' , default="20090101",       help='specify initial date as <YYYYMMDD>')
parser.add_argument('-l', '--lev' , type=int, default=0,      help='specify the level to use if 3D field (default: 0 => 2D)')
parser.add_argument('-z', '--zld' ,                           help='specify the topography netCDF file to use (field="z")')
parser.add_argument('-t', '--tstep', type=int, default=1, help='specify the time step (hours) in input file')

args = parser.parse_args()

CNEMO = args.conf
CBOX  = args.box
CWHAT = args.what
cf_in = args.fin
cf_mm = args.fmm
csd0  = args.sd0
jk    = args.lev
cf_topo_land = args.zld
dt    = args.tstep  ; # time step in hours

print('')
print(' *** CNEMO = ', CNEMO)
print(' *** CBOX  = ', CBOX)
print(' *** CWHAT = ', CWHAT)
print(' *** cf_in = ', cf_in)
print(' *** cf_mm = ', cf_mm)
print(' *** csd0 = ', csd0)
print(' ***   jk  = ', jk)
l_add_topo_land = False
if args.zld != None:
    print(' *** cf_topo_land = ', cf_topo_land)
    l_add_topo_land = True
l3d = False
if jk > 0:
    l3d=True
else:
    jk=0
###############################################################################################################################################

if not path.exists('figs'): mkdir('figs')
cdir_figs = './figs/'+CWHAT
if not path.exists(cdir_figs): mkdir(cdir_figs)

if l_save_nc and not path.exists('nc'): mkdir('nc')

CRUN = ''
if CNEMO == 'none':
    # Name of RUN:
    vv = split('-|_', path.basename(cfx_in))
    if vv[0] != CNEMO:
        print('ERROR: your file name is not consistent with '+CNEMO+' !!! ('+vv[0]+')'); sys.exit(0)
    CRUN = vv[1]
    print('\n Run is called: ''+CRUN+'' !\n')

#---------------------------------------------------------------

nemo_box = cp.nemo_hbox(CNEMO,CBOX)

(Ni0,Nj0) = nemo_box.size()
print(' '+CNEMO+': Ni0,Nj0 => ', Ni0,Nj0)

(i1,j1, i2,j2) = nemo_box.idx()
print(' i1,j1, i2,j2 => ', i1,j1, i2,j2,'\n')

if nemo_box.l_show_clock: (x_clock,y_clock) = nemo_box.clock
#print(' x_clock,y_clock =', x_clock,y_clock

if nemo_box.l_show_exp: (x_exp,y_exp) = nemo_box.exp
#print(' x_exp,y_exp =', x_exp,y_exp

if nemo_box.l_add_logo: (x_logo,y_logo) = nemo_box.logo

#---------------------------------------------------------------



cv_out = CWHAT ; # default

if   CWHAT == 'MLD':
    cv_in = 'somxl010' ; cv_out = 'MLD'
    tmin=0. ;  tmax=1800.  ;  df = 50. ; cpal_fld = 'ncview_hotres' ;     cb_jump = 4
    cunit = r'MLD [m]'
    if CBOX == 'Nordic': tmin=0. ; tmax=1000. ; cb_jump = 2 ; cpal_fld = 'magma_r' ; color_top_cb='k'

elif CWHAT == 'SST':
    cv_in = 'sosstsst' ; cv_out = CWHAT ; #in ['sosstsst','tos']:    
    tmin=-2 ;  tmax=32.   ;  df = 1. ; cpal_fld = 'ncview_nrl' ;     cb_jump = 2
    cunit = r'SST ($^{\circ}$C)'
    if CBOX == 'ALL':      l_do_ice = True
    if CBOX == 'EUROPA':   tmin=0. ;  tmax=25.
    if CBOX == 'EUROPAs':  tmin=6. ;  tmax=18.
    if CBOX == 'Med':      tmin=10.;  tmax=30. ; cb_jump = 1
    if CBOX == 'Med+BS':   tmin=7. ;  tmax=25.
    if CBOX == 'LabSea':   tmin=-2.;  tmax=15.
    if CBOX == 'Brittany': tmin=7. ;  tmax=13.
    if CBOX == 'GrlIcl':   tmax = 12.
    if CBOX in [ 'AzoresP','AzoresL','AzoresS']:  tmin = 15. ; tmax = 25. ; df=0.5
    if CBOX == 'Bretagne': tmin = 10. ; tmax = 22. ; df=1.
    
elif CWHAT == 'T_1000':
    cv_in = 'votemper' ; cv_out = CWHAT ;
    tmin=0. ;  tmax=14.   ;  df = 1. ; cpal_fld = 'ncview_nrl' ; cb_jump = 1
    cunit = r'Potential temperature at 1000 m'

elif CWHAT == 'T_60':
    cv_in = 'votemper' ; cv_out = CWHAT ;
    tmin=0. ;  tmax=14.   ;  df = 1. ; cpal_fld = 'ncview_nrl' ; cb_jump = 1
    cunit = r'Potential temperature at 60 m'
    if CBOX == 'BlackSea' : tmin=0. ; tmax=20. ;  df = 1.   ; cb_jump = 1 ; #cpal_fld = 'gist_stern_r'

elif CWHAT == 'SSS':
    cv_in = 'sosaline' ; cv_out = CWHAT ; #in ['sosstsst','tos']:    
    tmin=20. ;  tmax=40.   ;  df = 2. ; cpal_fld = 'ncview_ssec' ; cb_jump = 2
    cunit = r'Sea surface salinity'
    if CBOX == 'Med' :    tmin=33. ; tmax=39.5 ;  df = 0.25 ; cpal_fld = 'magma'
    if CBOX == 'LabSea' : tmin=28. ; tmax=35.5 ;  df = 1.   ; cb_jump = 1 ; cpal_fld = 'gist_stern_r'

elif CWHAT == 'S_1000':
    cv_in = 'vosaline' ; cv_out = CWHAT ;
    tmin=33.5 ;  tmax=36.5   ;  df = 0.5 ; cpal_fld = 'ncview_helix2' ; cb_jump = 1
    cunit = r'Salinity at 1000 m'
    if CBOX == 'MeddiesW' : tmin=35. ;  tmax=36.6 ; df = 0.1

elif CWHAT == 'GRAD_SST':
    cv_in = 'sosstsst' ; cv_out = CWHAT ;
    l_apply_hgrad = True
    l_smooth = True ; nb_smooth  = 5
    tmin=0. ;  tmax=0.001 ;  df = 0.0001 ; cpal_fld = 'ncview_hotres' ; cb_jump = 1
    cunit = r'$\left|\vec{\nabla}SST\right|$ (K/m)'

elif CWHAT == 'SSU':
    cv_in = 'sozocrtx' ; cv_out = CWHAT ;
    tmin=-0.8 ;  tmax=-tmin  ;  df = 0.1 ; cpal_fld = 'bone' ; cb_jump = 1
    cunit = r'Zonal component of current speed [m/s]'
    cv_msk = 'umask'

elif CWHAT == 'SSV':
    cv_in = 'somecrty' ; cv_out = CWHAT ;
    tmin=-1. ;  tmax=-tmin  ;  df = 0.2 ; cpal_fld = 'bone' ; cb_jump = 1
    cunit = r'Meridional component of current speed [m/s]'
    cv_msk = 'vmask'

elif CWHAT == 'SSH':
    cv_in = 'sossheig' ; cv_out = CWHAT ;
    cpal_fld = 'RdBu_r' ; tmin=-3. ;  tmax=-tmin   ;  df = 0.5 ;
    cb_jump = 1
    cunit = r'SSH [m]'
    if CBOX == 'Med' or CBOX == 'Med+BS': tmin=-0.7; tmax=0.2   ; df = 0.1
    if CRUN[:4] == 'BLB0':                tmin=-1.2; tmax=-tmin ; df = 0.2
    if CBOX in [ 'Bretagne']:             tmin=-4.;  tmax=-tmin ; df = 0.5

elif CWHAT == 'CURLOF':
    cv_in = 'socurloverf' ; cv_out = CWHAT ;
    #tmin=-0.8 ;  tmax=-tmin  ;  df = 0.1  ; cb_jump = 2 ;
    tmin=-0.7 ;  tmax=-tmin  ;  df = 0.1  ; cb_jump = 1
    cpal_fld='RdBu_r' ; color_top_cb='k' ; # cpal_fld = 'on2' 
    cunit = r'$\zeta/f$'
    cv_msk = 'vmask'

    
elif CWHAT == 'GEOSSV':
    # Geostrophic velocity speed out of SSH
    cv_in = 'sossheig' ; cv_out = CWHAT ;
    l_apply_geov = True
    cpal_fld = 'on3' ; tmin=0. ;  tmax=1.2   ;  df = 0.2 ; cb_extend = 'max'
    cb_jump = 1
    cunit = r'Surface geostrophic velocity speed [m/s]'
    l_save_nc = True

elif CWHAT == 'LAP_SSH':
    cv_in = 'sossheig' ; cv_out = CWHAT ;
    l_apply_lap = True
    cpal_fld = 'on2' ; tmin=-1.2 ;  tmax=1.2   ;  df = 0.05 ; 

elif CWHAT == 'W_1000':
    cv_in = 'vovecrtz'  ; cv_out = CWHAT ;
    tmin=-0.01 ;  tmax=-tmin   ;  df = 0.005 ; cb_jump = 1
    cpal_fld='RdBu_r' ;    #cpal_fld='PiYG_r' ; #cpal_fld='BrBG_r'
    cunit = r'Vertical velocity at 1000 m [m/s]'
    if CBOX in [ 'AzoresP','AzoresL','AzoresS']: color_top = 'k'

elif CWHAT == 'Amplitude':
    cv_in = 'r'     ; cv_out = cv_in ;
    cpal_fld = 'RdBu_r' ; tmin=-0.5 ;  tmax=-tmin   ;  df = 0.1 ; cb_jump = 1
    cunit = r'Amplitude [m]'

elif CWHAT == 'Phase':
    cv_in = 'phi'     ; cv_out = cv_in ;
    cpal_fld = 'RdBu_r' ; tmin=-30. ;  tmax=-tmin   ;  df = 5. ; cb_jump = 1
    #
    cunit = r'Phase (deg.)'
    
else:
    print('ERROR: we do not know variable ''+str(cv_in)+'' !')
    sys.exit(0)



    
rfz   = nemo_box.rfact_zoom
fontr = nemo_box.font_rat
    
print('\n================================================================')
print('\n rfact_zoom = ', rfz)
print(' font_rat = ', fontr, '\n')


print(' i1,i2,j1,j2 =>', i1,i2,j1,j2)

nx_res = i2-i1
ny_res = j2-j1

print(' *** nx_res, ny_res =', nx_res, ny_res)

yx_ratio = float(nx_res+1)/float(ny_res+1)

rnxr = rfz*nx_res ; # widt image (in pixels)
rnyr = rfz*ny_res ; # height image (in pixels)

# Target resolution for figure:
rh_fig = round(rnyr/float(rDPI),3) ; # width of figure
rw_fig = round(rh_fig*yx_ratio      ,3) ; # height of figure
rh_img = rh_fig*float(rDPI)
rw_img = rw_fig*float(rDPI)
while rw_img < round(rnxr,0):
    rw_fig = rw_fig + 0.01/float(rDPI)
    rw_img = rw_fig*float(rDPI)
while rh_img < round(rnyr,0):
    rh_fig = rh_fig + 0.01/float(rDPI)
    rh_img = rh_fig*float(rDPI)
    print(' *** size figure =>', rw_fig, rh_fig, '\n')
    print(' *** Forecasted dimension of image =>', rw_img, rh_img)

print('\n================================================================\n\n\n')



l_3d_field = False



# Ice:
if l_do_ice:
    cv_ice  = 'ileadfrac'
    cf_ice = replace(cf_in, 'gridT-2D', 'icemod')
    rmin_ice = 0.25
    #cpal_ice = 'ncview_bw'
    #cpal_ice = 'Blues_r'
    cpal_ice = 'bone'
    vcont_ice = nmp.arange(rmin_ice, 1.05, 0.05)
    #
    pal_ice = cp.chose_colmap(cpal_ice)
    norm_ice = colors.Normalize(vmin = rmin_ice, vmax = 0.95, clip = False)





if l_do_ice: cp.chck4f(cf_ice)

cp.chck4f(cf_mm)
cp.chck4f(cf_in)

l_notime=False
id_f = Dataset(cf_in)
list_var = id_f.variables.keys()
if 'time_counter' in list_var:
    vtime = id_f.variables['time_counter'][:]
elif 'time' in list_var:
    vtime = id_f.variables['time'][:]
else:
    l_notime=True
    print('Did not find a time variable! Assuming no time and Nt=1')
id_f.close()

Nt = 1
if not l_notime: Nt = len(vtime)

if l_show_lsm:
    cv_msk = 'tmask'
    print('\nReading record metrics in '+cf_mm)
    id_lsm = Dataset(cf_mm)
    nb_dim = len(id_lsm.variables[cv_msk].dimensions)
    print(' The mesh_mask has '+str(nb_dim)+' dimmensions!')
    print(' *** Reading '+cv_msk+' !')
    if nb_dim==4: XMSK = id_lsm.variables[cv_msk][0,jk,j1:j2,i1:i2]
    if nb_dim==3: XMSK = id_lsm.variables[cv_msk][jk,j1:j2,i1:i2]
    if nb_dim==2: XMSK = id_lsm.variables[cv_msk][j1:j2,i1:i2]
    (nj,ni) = nmp.shape(XMSK)
    
    if l_apply_lap:
        print(' *** Reading e1t and e2t !')
        XE1T2 = id_lsm.variables['e1t'][0,j1:j2,i1:i2]
        XE2T2 = id_lsm.variables['e2t'][0,j1:j2,i1:i2]
        XE1T2 = XE1T2*XE1T2
        XE2T2 = XE2T2*XE2T2
    if l_apply_hgrad or l_apply_geov:
        print(' *** Reading e1u and e2v !')
        e1u = id_lsm.variables['e1u'][0,j1:j2,i1:i2]
        e2v = id_lsm.variables['e2v'][0,j1:j2,i1:i2]
        print(' *** Reading umask and vmask !')
        if nb_dim==4:
            UMSK = id_lsm.variables['umask'][0,jk,j1:j2,i1:i2]
            VMSK = id_lsm.variables['vmask'][0,jk,j1:j2,i1:i2]
        if nb_dim==3:
            UMSK = id_lsm.variables['umask'][jk,j1:j2,i1:i2]
            VMSK = id_lsm.variables['vmask'][jk,j1:j2,i1:i2]
        if nb_dim==2:
            UMSK = id_lsm.variables['umask'][j1:j2,i1:i2]
            VMSK = id_lsm.variables['vmask'][j1:j2,i1:i2]
    
    if l_apply_geov:        
        ## Coriolis Parameter:
        ff  = id_lsm.variables['gphif'][0,j1:j2,i1:i2]
        ff[:,:] = 2.*romega*nmp.sin(ff[:,:]*nmp.pi/180.0)

    if l_save_nc:
        Xlon = id_lsm.variables['glamt'][0,j1:j2,i1:i2]
        Xlat = id_lsm.variables['gphit'][0,j1:j2,i1:i2]
        
    if l3d:
        vdepth = id_lsm.variables['gdept_1d'][0,:]
        zdepth = vdepth[jk]
        if zdepth>990. and zdepth<1010.: rof_log = 1000.
        rof_dpt = zdepth
        cdepth = str(round(zdepth,1))+'m'
    id_lsm.close()

    print('Shape Arrays => ni,nj =', ni,nj)

    print('Done!\n')



idx_land = nmp.where( XMSK <  0.5 )
#idx_ocea = nmp.where( XMSK >= 0.5 )

#(vjj_land, vji_land) = idx_land
#print(idx_land
#print(vjj_land
#print(vji_land
#print(len(vjj_land), len(vji_land)
#sys.exit(0)


XLSM = nmp.zeros((nj,ni)) ; # will define real continents not NEMO mask...

if l_add_topo_land:
    cp.chck4f(cf_topo_land)
    id_top = Dataset(cf_topo_land)
    print(' *** Reading 'z' into:\n'+cf_topo_land)
    xtopo = id_top.variables['z'][0,j1:j2,i1:i2]
    id_top.close()
    if nmp.shape(xtopo) != (nj,ni):
        print('ERROR: topo and mask do not agree in shape!'); sys.exit(0)
    xtopo = xtopo*(1. - XMSK)
    #cp.dump_2d_field('topo_'+CBOX+'.nc', xtopo, name='z')    
    if l3d: xtopo = xtopo + rof_dpt
    xtopo[nmp.where( XMSK > 0.01)] = nmp.nan
    if nemo_box.l_fill_holes_k and not l3d:
        XLSM[nmp.where( xtopo < 0.0)] = 1
        xtopo[nmp.where( xtopo < 0.0)] = nmp.nan

XLSM[nmp.where( XMSK > 0.5)] = 1

        
params = { 'font.family':'Helvetica Neue',
           'font.weight':    'normal',
           'font.size':       int(9.*fontr),
           'legend.fontsize': int(9.*fontr),
           'xtick.labelsize': int(9.*fontr),
           'ytick.labelsize': int(9.*fontr),
           'axes.labelsize':  int(9.*fontr) }
mpl.rcParams.update(params)
cfont_clb_tcks = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(7.5*fontr), 'color':color_top_cb}
cfont_clb  =  { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(8.5*fontr), 'color':color_top_cb}
cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(9.*fontr), 'color':color_top }
cfont_exp= { 'fontname':'Open Sans'  , 'fontweight':'light', 'fontsize':int(9.*fontr), 'color':color_top }
cfont_mail =  { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*fontr), 'color':'0.8'}
cfont_titl =  { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(30.*fontr), 'color':color_top }


# Colormaps for fields:
pal_fld = cp.chose_colmap(cpal_fld)
if   l_log_field:
    norm_fld = colors.LogNorm(                   vmin=tmin, vmax=tmax, clip=False)
elif l_pow_field:
    norm_fld = colors.PowerNorm(gamma=pow_field, vmin=tmin, vmax=tmax, clip=False)
else:
    norm_fld = colors.Normalize(                 vmin=tmin, vmax=tmax, clip=False)


if l_show_lsm or l_add_topo_land:
    if l_add_topo_land:
        xtopo = nmp.log10(xtopo+rof_log)
        pal_lsm = cp.chose_colmap('gray_r')
        #norm_lsm = colors.Normalize(vmin = nmp.log10(min(-100.+rof_dpt/3.,0.) + rof_log), vmax = nmp.log10(4000.+rof_dpt + rof_log), clip = False)
        norm_lsm = colors.Normalize(vmin = nmp.log10(-100. + rof_log), vmax = nmp.log10(4000.+rof_dpt + rof_log), clip = False)
    else:
        pal_lsm = cp.chose_colmap('land_dark')
        norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)

cyr0=csd0[0:4]
cmn0=csd0[4:6]
cdd0=csd0[6:8]

# Time step as a string
if not dt in [ 24, 6, 3, 1 ]:
    print('ERROR: unknown dt! '+str(dt))
    sys.exit(0)
ntpd = 24/dt

vm = vmn
if isleap(int(cyr0)): vm = vml
#print(' year is ', vm, nmp.sum(vm)

jd = int(cdd0) - 1
jm = int(cmn0)

# Opening netCDF files:
id_f = Dataset(cf_in)

for jt in range(jt0,Nt):

    #---------------------- Calendar stuff --------------------------------------------    
    #jh  = (jt*dt)%24
    jh  = int( (float(jt)+0.5)*float(dt) ) % 24 ; # average is centered
    if jt%ntpd == 0: jd = jd + 1
    if jd == vm[jm-1]+1 and (jt)%ntpd == 0 :
        jd = 1
        jm = jm + 1
    ch = '%2.2i'%(jh)
    cd = '%3.3i'%(jd)
    cm = '%2.2i'%(jm)
    ct = str(datetime.datetime.strptime(cyr0+'-'+cm+'-'+cd+' '+ch, '%Y-%m-%j %H'))    
    ct=ct[:5]+cm+ct[7:] #lolo bug !!! need to do that to get the month and not '01    
    cday  = ct[:10]   ; #print(' *** cday  :', cday        
    if dt >= 24:
        cdate = cday
        cdats = cday
    else:
        chour = ct[11:13] ; #print(' *** chour :', chour
        cdate = cday+'_'+chour
        cdats = cday+' '+chour+':00'
    print('\n Current date = ', cdate+' !\n')
    #-----------------------------------------------------------------------------------

    if l3d:
        cfig = cdir_figs+'/'+cv_out+'_'+CNEMO+'-'+CRUN+'_lev'+str(jk)+'_'+CBOX+'_'+cdate+'_'+cpal_fld+'.'+fig_type
    else:
        cfig = cdir_figs+'/'+cv_out+'_'+CNEMO+'-'+CRUN+'_'+CBOX+'_'+cdate+'_'+cpal_fld+'.'+fig_type


    if not path.exists(cfig):
    ###### FIGURE ##############

        fig = plt.figure(num = 1, figsize=(rw_fig, rh_fig), dpi=None, facecolor='w', edgecolor='0.5')
    
        ax  = plt.axes([0., 0., 1., 1.], facecolor = '0.85') # missing seas will be in 'facecolor' !
    
        vc_fld = nmp.arange(tmin, tmax + df, df)
    
    
        print('Reading record #'+str(jt)+' of '+cvx_in+' in '+cfx_in)
        if l3d: print('            => at level #'+str(jk)+' ('+cdepth+')!')
    
        if l_notime:
            if l3d:
                Xplot  = id_f.variables[cv_in][jk,j1:j2,i1:i2]
            else:
                Xplot  = id_f.variables[cv_in][j1:j2,i1:i2]
        else:
            if l3d:
                Xplot  = id_f.variables[cv_in][jt,jk,j1:j2,i1:i2] ; # t, y, x        
            else:
                Xplot  = id_f.variables[cv_in][jt,j1:j2,i1:i2] ; # t, y, x        
    
        print('Done!\n')
    
    
        # Ice
        if l_do_ice:
            print('Reading record #'+str(jt)+' of '+cv_ice+' in '+cf_ice)
            id_ice = Dataset(cf_ice)
            XICE  = id_ice.variables[cv_ice][jt,j1:j2,i1:i2] ; # t, y, x
            id_ice.close()
            print('Done!\n')
    
    
    
    
    
        if l_apply_lap:
            print(' *** Computing Laplacian of "'+cv_in+'"!')
            lx = nmp.zeros((nj,ni))
            ly = nmp.zeros((nj,ni))
            lx[:,1:ni-1] = 1.E9*(Xplot[:,2:ni] -2.*Xplot[:,1:ni-1] + Xplot[:,0:ni-2])/XE1T2[:,1:ni-1]
            ly[1:nj-1,:] = 1.E9*(Xplot[2:nj,:] -2.*Xplot[1:nj-1,:] + Xplot[0:nj-2,:])/XE2T2[1:nj-1,:]
            Xplot[:,:] = lx[:,:] + ly[:,:]
            del lx, ly
    
        if l_apply_hgrad:
            print(' *** Computing gradient of "'+cv_in+'"!')
            lx = nmp.zeros((nj,ni))
            ly = nmp.zeros((nj,ni))
    
            if l_smooth: cp.smoother(Xplot, XMSK, nb_smooth=nb_smooth)
            
            # Zonal gradient on T-points:
            lx[:,1:ni-1] = (Xplot[:,2:ni] - Xplot[:,0:ni-2]) / (e1u[:,1:ni-1] + e1u[:,0:ni-2]) * UMSK[:,1:ni-1] * UMSK[:,0:ni-2]
            lx[:,:] = XMSK[:,:]*lx[:,:]
            #cp.dump_2d_field('dsst_dx_gridT.nc', lx, xlon=Xlon, xlat=Xlat, name='dsst_dx')
            # Meridional gradient on T-points:
            ly[1:nj-1,:] = (Xplot[2:nj,:] - Xplot[0:nj-2,:]) / (e2v[1:nj-1,:] + e2v[0:nj-2,:]) * VMSK[1:nj-1,:] * VMSK[0:nj-2,:]
            ly[:,:] = XMSK[:,:]*ly[:,:]
            #cp.dump_2d_field('dsst_dy_gridT.nc', ly, xlon=Xlon, xlat=Xlat, name='dsst_dy')
            Xplot[:,:] = 0.0
            # Modulus of vector gradient:        
            Xplot[:,:] = nmp.sqrt(  lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] )
            #cp.dump_2d_field('mod_grad_sst.nc', Xplot, xlon=Xlon, xlat=Xlat, name='dsst')
            del lx, ly
    
    
        if l_apply_geov:
            print(' *** Computing gradient of "'+cv_in+'"!')
            lx = nmp.zeros((nj,ni))
            ly = nmp.zeros((nj,ni))
    
            # Zonal gradient on T-points:
            lx[:,1:ni-1] = (Xplot[:,2:ni] - Xplot[:,0:ni-2]) / (e1u[:,1:ni-1] + e1u[:,0:ni-2]) * UMSK[:,1:ni-1] * UMSK[:,0:ni-2]
            lx[:,:] = XMSK[:,:]*lx[:,:]
            #cp.dump_2d_field('dsst_dx_gridT.nc', lx, xlon=Xlon, xlat=Xlat, name='dsst_dx')
            # Meridional gradient on T-points:
            ly[1:nj-1,:] = (Xplot[2:nj,:] - Xplot[0:nj-2,:]) / (e2v[1:nj-1,:] + e2v[0:nj-2,:]) * VMSK[1:nj-1,:] * VMSK[0:nj-2,:]
            ly[:,:] = XMSK[:,:]*ly[:,:]
            #cp.dump_2d_field('dsst_dy_gridT.nc', ly, xlon=Xlon, xlat=Xlat, name='dsst_dy')
            Xplot[:,:] = 0.0
            # Modulus of vector gradient:        
            Xplot[:,:] = grav/ff * nmp.sqrt( lx[:,:]*lx[:,:] + ly[:,:]*ly[:,:] )
            #cp.dump_2d_field('mod_grad_sst.nc', Xplot, xlon=Xlon, xlat=Xlat, name='dsst')
            del lx, ly
    
    
            
        print('')
        if not l_show_lsm and jt == jt0: ( nj , ni ) = nmp.shape(Xplot)
        print('  *** dimension of array => ', ni, nj, nmp.shape(Xplot))
        

        print('Ploting')
    
        plt.axis([ 0, ni, 0, nj])
        
        if nemo_box.c_imshow_interp == 'none':
            Xplot[idx_land] = nmp.nan
        else:
            cp.drown(Xplot, XMSK, k_ew=-1, nb_max_inc=10, nb_smooth=10)
    
        if l_save_nc:
            if l3d:
                cf_out = 'nc/'+CWHAT+'_NEMO_'+CNEMO+'-'+CRUN+'_lev'+str(jk)+'_'+CBOX+'_'+cdate+'_'+cpal_fld+'.nc'
            else:
                cf_out = 'nc/'+CWHAT+'_NEMO_'+CNEMO+'-'+CRUN+'_'+CBOX+'_'+cdate+'_'+cpal_fld+'.nc'
            print(' Saving in '+cf_out)
            cp.dump_2d_field(cf_out, Xplot, xlon=Xlon, xlat=Xlat, name=CWHAT)
            print('')
    
    
        cf = plt.imshow( Xplot[:,:], cmap=pal_fld, norm=norm_fld, interpolation=nemo_box.c_imshow_interp )
    
        # Ice
        if l_do_ice:
            #XM[:,:] = XMSK[:,:]
            #cp.drown(XICE, XM, k_ew=2, nb_max_inc=10, nb_smooth=10)
            #ci = plt.contourf(XICE[:,:], vcont_ice, cmap = pal_ice, norm = norm_ice) #
            pice = nmp.ma.masked_where(XICE < rmin_ice, XICE)
            ci = plt.imshow(pice, cmap = pal_ice, norm = norm_ice, interpolation='none') ; del pice, ci
            del XICE
    
        #LOLO: rm ???
        if l_show_lsm or l_add_topo_land:
            if l_add_topo_land:
                clsm = plt.imshow(nmp.ma.masked_where(XLSM>0.0001, xtopo), cmap = pal_lsm, norm = norm_lsm, interpolation='none')
                if nemo_box.c_imshow_interp == 'none':
                    plt.contour(XLSM, [0.9], colors='k', linewidths=0.5)
            else:
                clsm = plt.imshow(nmp.ma.masked_where(XLSM>0.0001, XLSM), cmap = pal_lsm, norm = norm_lsm, interpolation='none')
    
        ##### COLORBAR ######
        if nemo_box.l_show_cb:
            ax2 = plt.axes(nemo_box.vcb)
            clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=norm_fld, orientation='horizontal', extend=cb_extend)
            cb_labs = []
            #if cb_jump > 1:
            cpt = 0
            for rr in vc_fld:
                if cpt % cb_jump == 0:
                    if df >= 1.: cb_labs.append(str(int(rr)))
                    if df <  1.: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./df)))+1) ))
                else:
                    cb_labs.append(' ')
                cpt = cpt + 1
            #else:
            #    for rr in vc_fld: cb_labs.append(str(round(rr,int(nmp.ceil(nmp.log10(1./df)))+1) ))
    
            clb.ax.set_xticklabels(cb_labs, **cfont_clb_tcks)
            clb.set_label(cunit, **cfont_clb)
            clb.ax.yaxis.set_tick_params(color=color_top_cb) ; # set colorbar tick color
            clb.outline.set_edgecolor(color_top_cb) ; # set colorbar edgecolor
            clb.ax.tick_params(which = 'minor', length = 2, color = color_top_cb )
            clb.ax.tick_params(which = 'major', length = 4, color = color_top_cb )
    
    
        if nemo_box.l_show_clock:
            xl = float(x_clock)/rfz
            yl = float(y_clock)/rfz
            ax.annotate('Date: '+cdats, xy=(1, 4), xytext=(xl,yl), **cfont_clock)
    
        if nemo_box.l_show_exp:
            xl = float(x_exp)/rfz
            yl = float(y_exp)/rfz
            ax.annotate('Experiment: '+CNEMO+'-'+CRUN, xy=(1, 4), xytext=(xl,yl), **cfont_exp)
    
    
        #ax.annotate('laurent.brodeau@ocean-next.fr', xy=(1, 4), xytext=(xl+150, 20), **cfont_mail)
    
        if nemo_box.l_annotate_name:
            xl = rnxr/20./rfz
            yl = rnyr/1.33/rfz
            ax.annotate(CNEMO, xy=(1, 4), xytext=(xl, yl), **cfont_titl)
    
        if nemo_box.l_add_logo:
            datafile = cbook.get_sample_data(dir_logos+'/'+nemo_box.cf_logo_on, asfileobj=False)
            im = image.imread(datafile)
            #im[:, :, -1] = 0.5  # set the alpha channel
            fig.figimage(im, x_logo, y_logo, zorder=9)
            del datafile, im
            #
            if nemo_box.l_add_logo_ige:
                datafile = cbook.get_sample_data(dir_logos+'/'+nemo_box.cf_logo_ige, asfileobj=False)
                im = image.imread(datafile)
                fig.figimage(im, x_logo+144, y_logo-150., zorder=9)
                del datafile, im
                #
            if nemo_box.l_add_logo_prc:
                datafile = cbook.get_sample_data(dir_logos+'/'+nemo_box.cf_logo_prc, asfileobj=False)
                im = image.imread(datafile)
                fig.figimage(im, x_logo-77, y_logo-140., zorder=9)
                del datafile, im
    
        plt.savefig(cfig, dpi=rDPI, orientation='portrait', facecolor='k')
        print(cfig+' created!\n')
        plt.close(1)

        if l_show_lsm: del clsm
        del cf, fig, ax
        if nemo_box.l_show_cb: del clb

    else:
        print('\n Figure '+cfig+' already there!\n')

id_f.close()


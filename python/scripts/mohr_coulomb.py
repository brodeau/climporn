#!/usr/bin/env python3
#
#     CLIMPORN
#
#  Shows:
#  SI3 output + mesh_mask needed.
#
#    L. Brodeau, March 2022

import sys
from os import path, getcwd, mkdir
import argparse as ap
import numpy as np

from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from calendar import isleap
import datetime

from re import split
#import warnings
#warnings.filterwarnings("ignore")

import climporn as cp
from climporn import fig_style_mov as fsm

vmn = [ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]
vml = [ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 ]

fig_type='png'
rDPI = 100
rfig_fact=1.2

color_top = 'k'
#color_top_cb = 'white'

cpal_fld = 'viridis'

tmin = 0.      # Pa
tmax = 40000.  # Pa
df=1000. # Pa

xmin = -30000. ; xmax = 10000
ymin = 0.      ; ymax = 30000


cv_in1 = 'ice_sigI-t'
cv_in2 = 'ice_sigII-t'

#l_get_name_of_run = True


# Normally logos should be found there:
#dir_logos = str.replace( path.dirname(path.realpath(__file__)) , 'climporn/python', 'climporn/misc/logos' )
#print("\n --- logos found into : "+dir_logos+" !\n")


################## ARGUMENT PARSING / USAGE ################################################################################################
parser = ap.ArgumentParser(description='Shows the normal and shear invariants of the stress tensor w.r.t Mohr-Coulomb line')

requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--fin' , required=True, help='specify the NEMO netCDF file to read from...')
requiredNamed.add_argument('-m', '--fmm' , required=True, help='specify the NEMO mesh_mask file (ex: mesh_mask.nc)')

#requiredNamed.add_argument('-w', '--what', required=True, default="sst", help='specify the field/diagnostic to plot (ex: sst)')
#parser.add_argument('-C', '--conf', default="NANUK025",     help='specify NEMO config (ex: eNATL60)')
#parser.add_argument('-N', '--name', default="X",            help='specify experiment name')
parser.add_argument('-s', '--sd0' , default="20160101",     help='specify initial date as <YYYYMMDD>')
#parser.add_argument('-l', '--lev' , type=int, default=0,    help='specify the level to use if 3D field (default: 0 => 2D)')
#parser.add_argument('-I', '--ice' , action='store_true',    help='draw sea-ice concentration layer onto the field')
#parser.add_argument('-T', '--title', default="",            help='specify experiment title')
parser.add_argument('-f', '--freq',  default="1h",          help='specify output frequency in input file (1h,3h,6h,1d, etc.)')
parser.add_argument('-B', '--box',   default=None,             help='specify rectangular box to focus => "-B i1,j1,i2,j2"')

args = parser.parse_args()

cf_in = args.fin
cf_mm = args.fmm
cfreq = args.freq
csd0  = args.sd0
cbox  = args.box

lbox = False
if cbox:
    lbox = True
    [i1,j1,i2,j2]   = np.array( split(',',cbox), dtype=int )
    print(" *** Box to restric diag to: i1,j1,i2,j2 =",i1,j1,i2,j2)

print('')
print(' *** cf_in = ', cf_in)
print(' *** cf_mm = ', cf_mm)
###############################################################################################################################################

cp.chck4f(cf_mm)

if not path.exists("figs"): mkdir("figs")
cdir_figs = './figs'
if not path.exists(cdir_figs): mkdir(cdir_figs)

l_notime=False
cp.chck4f(cf_in)
id_in = Dataset(cf_in)
Nj, Ni = id_in.dimensions['y'].size, id_in.dimensions['x'].size
print(' *** Ni, Nj = ',Ni, Nj,'\n')
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


if not lbox:
    # Use the whole domain
    i1, i2 = 0, Ni
    j1, j2 = 0, Nj


id_lsm = Dataset(cf_mm)
n2, n1 = id_lsm.dimensions['y'].size, id_lsm.dimensions['x'].size
if n2!=Nj and n1!=Ni:
    print('ERROR: fields in meshmask file do not have the same shape as those in input file!'); exit(0)
nb_dim = len(id_lsm.variables['tmask'].dimensions)
print(' The mesh_mask has '+str(nb_dim)+' dimmensions!')
print(' *** Reading '+'tmask'+' !')
if nb_dim==4: XMSK = id_lsm.variables['tmask'][0,0,j1:j2,i1:i2]
if nb_dim==3: XMSK = id_lsm.variables['tmask'][0,j1:j2,i1:i2]
if nb_dim==2: XMSK = id_lsm.variables['tmask'][j1:j2,i1:i2]
Xlon = id_lsm.variables['glamt'][0,j1:j2,i1:i2]
Xlat = id_lsm.variables['gphit'][0,j1:j2,i1:i2]
id_lsm.close()




(Nj,Ni) = np.shape(XMSK) ; # update with shape of box we read...


VSIG1 = []
VSIG2 = []

id_in = Dataset(cf_in)

# A/ store everything in arrays for all the time steps and plot the result:
for jt in range(Nt):

    XS1 = id_in.variables[cv_in1][jt,j1:j2,i1:i2]
    XS2 = id_in.variables[cv_in2][jt,j1:j2,i1:i2]

    for jj in range(0,Nj):
        for ji in range(0,Ni):
            zS1, zS2 = XS1[jj,ji], XS2[jj,ji]
            if zS2<0:
                print('ERROR: WTF?'); exit(0)
            if XMSK[jj,ji]>0.5 and (zS2>0. and abs(zS1)>0):
                VSIG1.append(zS1)
                VSIG2.append(zS2)
    
id_in.close()

VSIG1 = np.array(VSIG1)
VSIG2 = np.array(VSIG2)

(Np,) = np.shape(VSIG1)


print(' *** Selected valid '+str(Np)+' points out of '+str(Nj*Ni*Nt)+'!')


# Font style:
font_ratio = 1.7*rfig_fact

kk = fsm( font_ratio );#r_top=fa.color_top, clr_top_cb=fa.color_top_cb )

col_bg = 'w'
vfig_size = [ 10.*rfig_fact, 7.2*rfig_fact ]
vsporg = [0.13, 0.13, 0.8, 0.8]

cfig = 'test.png'

fig = plt.figure(num = 1, figsize=(vfig_size), dpi=None, facecolor=col_bg, edgecolor=col_bg)
ax  = plt.axes(vsporg, facecolor = col_bg)

#XS1 = np.ma.masked_where(XMSK <= 0.1, XS1)
#XS2 = np.ma.masked_where(XMSK <= 0.1, XS2)

#ft = carte.pcolormesh(x0, y0, XS1, cmap = pal_fld, norm = nrm_fld )


plt.axis([ xmin, xmax, ymin, ymax])

sp = plt.scatter( VSIG1, VSIG2, s=1., c='k', alpha=0.1)

vx = np.arange(xmin,xmax+1000.,1000.)
vy = -0.7 * vx[:] + 6000.
#print(" vx =", vx[:] )
#sys.exit(0)

mc = plt.plot( vx, vy, color='#041a4d', linewidth=4. )

plt.xlabel(r'$\sigma_{I}$')
plt.ylabel(r'$\sigma_{II}$')

#ax.annotate('Date: '+cday+' '+chour+':00', xy=(0.785, 0.015), xycoords='figure fraction', **fsm.cfont_clock)
#ax.annotate('Date: '+cday+' '+chour+':00', xy=(0.5, 0.95), xycoords='figure fraction', **fsm.cfont_clock)

#ry0 = 0.78
#ax.annotate(cxtra_info1, xy=(0.02, ry0+0.05), xycoords='figure fraction', **fsm.cfont_titl1)

#if ctitle != "":
#    ax.annotate(ctitle, xy=(0.03, ry0-0.01), xycoords='figure fraction', **fsm.cfont_titl2)


print(' Saving figure: '+cfig+'\n\n')

plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)

plt.close(1)


exit(0)







#fontr=1.2*rfig_fact
#params = { 'font.family':'Open Sans',
#           'font.weight':    'normal',
#           'font.size':       int(12.*fontr),
#           'legend.fontsize': int(12.*fontr),
#           'xtick.labelsize': int(12.*fontr),
#           'ytick.labelsize': int(12.*fontr),
#           'axes.labelsize':  int(15.*fontr) }
#mpl.rcParams.update(params)
###cfont_clb  =  { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(10*fontr), 'color':color_top_cb}
#cfont_clock = { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':int(13*fontr) , 'color':color_top }
##cfont_exp= { 'fontname':'Open Sans'  , 'fontweight':'light', 'fontsize':int(9.*fontr), 'color':color_top }
##cfont_mail  =  { 'fontname':'Times New Roman', 'fontweight':'normal', 'fontstyle':'italic', 'fontsize':int(14.*fontr), 'color':'0.8'}
#cfont_titl1 = { 'fontname':'Open Sans', 'fontweight':'light', 'fontsize':int(18.*fontr), 'color':color_top }
#cfont_titl2 = { 'fontname':'Open Sans', 'fontweight':'normal','fontsize':int(12.*fontr), 'color':color_top }


# for colorbar:
vc_fld = np.arange(tmin, tmax + df, df)


# For movie
vfig_size = [ 10.*rfig_fact, 7.2*rfig_fact ]
vsporg = [0.13, 0.13, 0.8, 0.8]
#vsporg = [0., 0.0001, 1., 1.001]
#vsporg = [0., 0., 1., 1.]


vcbar = [0.05, 0.065, 0.92, 0.03]

cyr0=csd0[0:4]
cmn0=csd0[4:6]
cdd0=csd0[6:8]

if cfreq == '6h':
    dt = 6
elif cfreq == '3h':
    dt = 3
elif cfreq == '1h':
    dt = 1
else:
    print('ERROR: unknown frequency!'); sys.exit(0)

pal_fld = cp.chose_colmap(cpal_fld)
nrm_fld = colors.Normalize(vmin=tmin, vmax=tmax, clip=False)

ntpd = 24/dt

vm = vmn
if isleap(int(cyr0)): vm = vml

jd = int(cdd0) - 1
jm = int(cmn0)




for jt in range(Nt):

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

    cfig = cdir_figs+'/MC_SI3-BBM_'+cday+'_'+chour+'.'+fig_type

    # Getting field and sea-ice concentration at time record "jt":

    for jj in range(0,Nj):
        VSIG1.append(XS1[jj,:])
        VSIG2.append(XS2[jj,:])

    
    cjt = '%4.4i'%(jt)

    #######################################################################################################
    #col_bg = '#3b3b63'
    #col_bg = '#041a4d'
    col_bg = 'w'
    fig = plt.figure(num = 1, figsize=(vfig_size), dpi=None, facecolor=col_bg, edgecolor=col_bg)
    ax  = plt.axes(vsporg, facecolor = col_bg)

    XS1 = np.ma.masked_where(XMSK <= 0.1, XS1)
    XS2 = np.ma.masked_where(XMSK <= 0.1, XS2)

    #ft = carte.pcolormesh(x0, y0, XS1, cmap = pal_fld, norm = nrm_fld )


    plt.axis([ xmin, xmax, ymin, ymax])
    
    sp = plt.scatter( XS1, XS2, s=1., c='k', alpha=0.5)

    vx = np.arange(xmin,xmax+1000.,1000.)
    vy = -0.7 * vx[:] + 6000.
    #print(" vx =", vx[:] )
    #sys.exit(0)
    
    mc = plt.plot( vx, vy, color='#041a4d', linewidth=4. )

    plt.xlabel(r'$\sigma_{I}$')
    plt.ylabel(r'$\sigma_{II}$')

    
    # ----------- Color bar for field -----------
    #ax2 = plt.axes([0.64, 0.965, 0.344, 0.018])
    #clb = mpl.colorbar.ColorbarBase(ax2, ticks=vc_fld, cmap=pal_fld, norm=nrm_fld, orientation='horizontal', extend='both')
    #cb_labs = []
    #cpt = 0
    #for rr in vc_fld:
    #    if cpt % cb_jump == 0 or ( (tmin == -tmax) and (int(rr) == 0 ) ):
    #        if df >= 1.: cb_labs.append(str(int(rr)))
    #        if df <  1.: cb_labs.append(str(round(rr,int(np.ceil(np.log10(1./df)))+1) ))
    #    else:
    #        cb_labs.append(' ')
    #    cpt = cpt + 1
    #clb.ax.set_xticklabels(cb_labs, **fsm.cfont_clb)
    #clb.set_label(cunit, **fsm.cfont_clb)
    #clb.ax.yaxis.set_tick_params(color=color_top_cb) ; # set colorbar tick color
    #clb.outline.set_edgecolor(color_top_cb) ; # set colorbar edgecolor
    #clb.ax.tick_params(which = 'minor', length = 2, color = color_top_cb )
    #clb.ax.tick_params(which = 'major', length = 4, color = color_top_cb )
    #del clb


    #ax.annotate('Date: '+cday+' '+chour+':00', xy=(0.785, 0.015), xycoords='figure fraction', **fsm.cfont_clock)
    ax.annotate('Date: '+cday+' '+chour+':00', xy=(0.5, 0.95), xycoords='figure fraction', **fsm.cfont_clock)

    #ry0 = 0.78
    #ax.annotate(cxtra_info1, xy=(0.02, ry0+0.05), xycoords='figure fraction', **fsm.cfont_titl1)

    #if ctitle != "":
    #    ax.annotate(ctitle, xy=(0.03, ry0-0.01), xycoords='figure fraction', **fsm.cfont_titl2)


    print(' Saving figure: '+cfig+'\n\n')

    plt.savefig(cfig, dpi=rDPI, orientation='portrait', transparent=False)
    plt.close(1)

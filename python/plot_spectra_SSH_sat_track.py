#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

#    L. Brodeau, 2018

import sys
from os import path, getcwd, mkdir
import argparse as ap
import numpy as nmp
#
import scipy.signal as signal
from netCDF4 import Dataset
#
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import time

import clprn_plot as bp
import clprn_tool as bt

l_tapper      = True ; # apply tappering !
l_detrend_lin = True ; # apply a linear detrending on data segment before computing spectrum...
l_rm_i_noise  = False

# Plots for each segment (only to debug!)
l_plot_rawd = True
l_plot_maps = False
l_plot_trck = False ; # plots individual model and satellite tracks...
l_plot_spct = False

dir_figs='./figs'
fig_ext='svg'
#fig_ext='png'


rcut_by_time = 1.2 # specify in seconds the time gap between two obs to detect and discontinuity and therefore
#                  # should be considered as a cut!

rcut_by_dist = 7.8 # same as rcut_by_time, but in terms of distance (in km) between two consecutive points!
#                  # time criterion "rcut_by_time" would have been sufficient if SARAL data was not bugged!!!
#                  # => in SARAL data, spotted two consecutive points with the usual dt and a huge dL !!!
#                  #   => like if the satellite undergone an extremely short huge acceleration !!!
#                  #   => ex: 3rd of August 2016 around 07:53:43 !!!

#nlen_valid_seg = 120  # specify the minimum number of values a segment should contain to be considered and used!
nlen_valid_seg = 18  # specify the minimum number of values a segment should contain to be considered and used!

r_max_amp_ssh = 1.5 # in meters

r2Pi = 2.*nmp.pi


# Look and feel for the plot:

clr_red = '#AD0000'

# Color OceanNext style:
clr_sat = '#ffed00'
clr_mod = '#008ab8'

rDPI=100.


#=============== en of configurable part ================================

if not path.exists(dir_figs): mkdir(dir_figs)



################## ARGUMENT PARSING / USAGE ################################################################################################

# Defaults before reading command-line arguments:
cn_mod = "NEMO"
cn_sat = "Altimetry"
cn_box = "Unknown box"
pow10_min_y = -8
pow10_max_y =  2

parser = ap.ArgumentParser(description='Generate along-track spectral comparison of SSH, model versus satellite')
#
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-i', '--fin' , required=True, help='specify output "result" file from interp_to_ground_track.x@SOSIE')
requiredNamed.add_argument('-m', '--vmod', required=True, help='specify the name of the model variable for SSH' )
requiredNamed.add_argument('-s', '--vsat', required=True, help='specify the name of the satellite variable for SSH')
#
parser.add_argument('-B', '--name_box', default=cn_box,   help='name of the rectangular region (box) considered')
parser.add_argument('-M', '--name_mod', default=cn_mod,   help='name of the model (ex: NEMO-eNATL60)')
parser.add_argument('-S', '--name_sat', default=cn_sat,   help='name of the satellite (ex: SARAL-Altika)')
#parser.add_argument('-z', '--zld' ,                           help='specify the topography netCDF file to use (field="z")')
parser.add_argument('-a', '--pmin_y', type=int, default=pow10_min_y, help='minimum y-axis value to display in terms of 10^pmin_y')
parser.add_argument('-b', '--pmax_y', type=int, default=pow10_max_y, help='maximum y-axis value to display in terms of 10^pmax_y')
#
args = parser.parse_args()

cf_in  = args.fin
cv_mod = args.vmod
cv_sat = args.vsat
cn_box = args.name_box
cn_mod = args.name_mod
cn_sat = args.name_sat
pow10_min_y = args.pmin_y
pow10_max_y = args.pmax_y


cfs  = path.basename(cf_in)

cseas = ''
if ('JFM' in cfs) and not('JAS' in cfs) : cseas = 'JFM'; vseas = ['01','02','03']
if ('JAS' in cfs) and not('JFM' in cfs) : cseas = 'JAS'; vseas = ['07','08','09']
print('\n *** Season: '+cseas+'\n')



cextra = ''

print(' *** Opening file '+cf_in+'!')
id_in    = Dataset(cf_in)
vt_epoch = id_in.variables['time'][:]
vmodel   = id_in.variables[cv_mod][:]
vsatel   = id_in.variables[cv_sat][:]
vlon     = id_in.variables['longitude'][:]
vlat     = id_in.variables['latitude'][:]
vdist    = id_in.variables['distance'][:]
id_in.close()
print("  => Everything read!\n")


# Debug: wave-signal with a wave-length of 150 km and 25 km (vdist in km)
#vmodel = 0.25*nmp.cos(vdist*r2Pi/150.) + 0.75*nmp.cos(vdist*r2Pi/25.)


nbr = len(vt_epoch)

# Initial raw plot:
# Create Matplotlib time array:


if l_plot_rawd or l_plot_trck:
    import matplotlib.dates as mdates


cyear = time.strftime("%Y", time.localtime(vt_epoch[2]))




ii=nbr//300
ib=max(ii-ii%10,1)
xticks_d=30.*ib


if l_plot_rawd:
    vtime = nmp.zeros(nbr)
    for jt in range(nbr): vtime[jt] = mdates.epoch2num(vt_epoch[jt])
    #
    fig = plt.figure(num = 1, figsize=(12,7), facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.07, 0.24, 0.9, 0.75])
    ax1.set_xticks(vtime[::xticks_d])
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    plt.xticks(rotation='60')
    plt.plot(vtime, vsatel, '-', color=clr_sat, linewidth=2, label=cn_sat+' ("'+cv_sat+'")', zorder=10)
    plt.plot(vtime, vmodel, '-', color=clr_mod, linewidth=2,  label=cn_mod+' ("'+cv_mod+'")', zorder=15)
    ax1.set_ylim(-r_max_amp_ssh,r_max_amp_ssh) ; ax1.set_xlim(vtime[0],vtime[nbr-1])
    plt.xlabel('Time [seconds since 1970]')
    plt.ylabel('SSH [m]')
    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(bbox_to_anchor=(0.55, 0.98), ncol=1, shadow=True, fancybox=True)
    plt.savefig(dir_figs+'/fig_raw_data_'+cn_mod+'--'+cn_sat+'.'+fig_ext, dpi=120, transparent=True)
    plt.close(1)






vmask = vmodel.mask

(idx_ok,) = nmp.where(vmask==False) # indexes with valid values!

#print(vmodel[:])
#print(idx_ok)
#sys.exit(0)

nbr_v = len(idx_ok)

print(' *** '+str(nbr_v)+' valid points out of '+str(nbr)+' !')

# Will extract the N valid data segments:
nb_seg=0
idx_seg_start = [] ; # index of first valid point of the segment
idx_seg_stop  = [] ; # index of last  valid point of the segment

jr=0
while jr < nbr:
    # Ignoring masked values and zeros...
    if (not vmask[jr]) and (vmodel[jr]!=0.0) and (vmodel[jr]<100.):
        nb_seg = nb_seg + 1
        print('\n --- found seg #'+str(nb_seg)+' !')
        idx_seg_start.append(jr)
        print(' => starting at jt='+str(jr))
        #while (not vmask[jr+1]) and (vmodel[jr+1]!=0.0) and (vt_epoch[jr+1]-vt_epoch[jr] < rcut_by_time) and (vmodel[jr]<100.) :
        while (not vmask[jr+1]) and (vmodel[jr+1]!=0.0) and (vt_epoch[jr+1]-vt_epoch[jr] < rcut_by_time) and (vdist[jr+1]-vdist[jr] < rcut_by_dist) and (vmodel[jr]<100.) :
            jr = jr+1
            if jr==nbr-1: break
        idx_seg_stop.append(jr)
        print(' => and stoping at jt='+str(jr))
    jr = jr+1

if len(idx_seg_start) != nb_seg: print(' ERROR #1!'); sys.exit(1)



# Maximum number of poins in the segments:

isd = nmp.asarray(idx_seg_stop[:]) - nmp.asarray(idx_seg_start[:]) + 1
print('\n lengths =>', isd[:])
nbp_max = max(isd)
print('\n *** Longest segments has nbp_max='+str(nbp_max)+' points!\n')

# Treat only segments with at least nlen_valid_seg points:
(vtreat,) = nmp.where(isd >= nlen_valid_seg)
nb_v_seg = len(vtreat)

if nb_v_seg==0:
    print('PROBLEM: could not find any valid segment with nlen_valid_seg = '+str(nlen_valid_seg)+' !')
    sys.exit(0)        

rN = nmp.mean(isd[vtreat])
print('\nMean segment-length for the '+str(nb_v_seg)+' segments with at least '+str(nlen_valid_seg)+' points:', round(rN,1))

Nsp = int(rN/10.)*10
print('  ==> Nsp = '+str(Nsp))
(vtreat,) = nmp.where(isd >= Nsp)
nb_v_seg = len(vtreat)
print(' ==> will use '+str(nb_v_seg)+' segments with a fixed length of '+str(Nsp)+' point!\n')


#print(' *** will treat '+str(nb_v_seg)+' segments out of '+str(nb_seg)+' (need at least '+str(nlen_valid_seg)+' points)')


x_all_spectra_s = nmp.zeros((nb_v_seg,Nsp)) # array to store all spectra in...
x_all_spectra_m = nmp.zeros((nb_v_seg,Nsp)) # array to store all spectra in...

clabel_mod=cn_mod+' ("'+cv_mod+'")'
clabel_sat=cn_sat+' ("'+cv_sat+'")'
#clabel_mod=cn_mod+'  ("'+cv_mod+'"), 1D along-track'
#clabel_mod=cn_mod+' (1D along-track)'
#clabel3=cn_mod+' (2D box)'






jcpt = -1
for js in vtreat:
    
    jcpt= jcpt+1
    it1 = idx_seg_start[js]
    it2 = idx_seg_stop[js]
    nbp = it2-it1+1    
    cseg = '%2.2i'%(js+1)

    print('\n\n ###################################')
    print('  *** Seg #'+cseg+' of '+cn_box+':')
    print('  ***   => originally '+str(nbp)+' points in this segment (from '+str(it1)+' to '+str(it2)+')')

    # nb of points in excess / Nsp:
    nxcs = nbp - Nsp
    jmp_strt = nxcs//2
    jmp_stop = nxcs//2 + nxcs%2
    it1 = it1+jmp_strt
    it2 = it2-jmp_stop
    print('  ***   => we only retain '+str(it2-it1+1)+' points from '+str(it1)+' to '+str(it2)+'!')
    
    # Checking the typical distance (in km) between two measures:
    dmean = nmp.mean(vdist[it1+1:it2+1]-vdist[it1:it2])
    print('\n Mean distance between two consecutive points is '+str(dmean)+' km\n')
    # Sample spacing in [km] (inverse of the sampling rate):
    rdist_sample = round(dmean,3)
    print(' => will use a spatial sample spacing of '+str(rdist_sample)+' km\n')

    if l_plot_trck:
        # Create Matplotlib time array:
        vtime = nmp.zeros(Nsp)
        for jt in range(Nsp): vtime[jt] = mdates.epoch2num(vt_epoch[it1+jt])

    # First centering the two time-series about 0, and tappering at extremities (filling with zeros)
    vs_s = nmp.zeros(Nsp) ; vs_m = nmp.zeros(Nsp)
    print('             (length = '+str(len(vsatel[it1:it2+1]))+' / '+str(Nsp)+')')
    
    vs_s[:] = vsatel[it1:it2+1]
    vs_m[:] = vmodel[it1:it2+1]

    # Linear detrending
    if l_detrend_lin:
        vs_s[:] = signal.detrend(vs_s[:],type='linear')
        vs_m[:] = signal.detrend(vs_m[:],type='linear')

    # Centering about 0:
    vs_s = vs_s - nmp.mean(vs_s)
    vs_m = vs_m - nmp.mean(vs_m)

    # Tappering:
    if l_tapper:        
        wdw =  signal.tukey(Nsp,0.5)
        vs_s = vs_s*wdw
        vs_m = vs_m*wdw
        
    # Power Spectrum:
    vYf_s = 2.*(rdist_sample/float(Nsp)) * nmp.abs(nmp.fft.fft(vs_s))**2
    vYf_m = 2.*(rdist_sample/float(Nsp)) * nmp.abs(nmp.fft.fft(vs_m))**2

    # Wave numbers:
    if jcpt==0:
        vk  = nmp.fft.fftfreq(Nsp, rdist_sample)
        idx = nmp.argsort(vk)

    # Saving for current segment:
    x_all_spectra_s[jcpt,:] = vYf_s[idx]
    x_all_spectra_m[jcpt,:] = vYf_m[idx]



    
    # Maps:
    if l_plot_maps:
        from mpl_toolkits.basemap import Basemap
        fig = plt.figure(num = 1, figsize=(12,8.3), facecolor='w', edgecolor='k')
        ax1 = plt.axes([0.01, 0.01, 0.98, 0.98])
        m = Basemap(projection='lcc', llcrnrlat=18, urcrnrlat=60, llcrnrlon=-77, urcrnrlon=+25,lat_0=45.,lon_0=-40., resolution='c', area_thresh=1000.)
        m.bluemarble()
        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        m.drawstates(linewidth=0.5)
        x, y = m(vlon[it1:it2+1], vlat[it1:it2+1])
        plt.plot(x, y, 'r-')
        plt.savefig(dir_figs+'/'+cn_mod+'--'+cn_sat+'_seg'+cseg+'_map.'+fig_ext, dpi=120, transparent=True)
        plt.close(1)


    if l_plot_trck:
        cfigure = dir_figs+'/'+cn_box+'_'+cseas+'_'+cn_mod+'--'+cn_sat+'_seg'+cseg+'_track.'+fig_ext
        ii=Nsp//300
        ib=max(ii-ii%10,1)
        xticks_d=30*ib
        # Finding appropriate amplitude as a multiple of 0.25:    
        rmult = 0.2
        rmax    = max( nmp.max(vs_m[:]) , nmp.max(vs_s[:])  )
        r_ssh_p = bt.round_to_multiple_of(rmax, prec=1, base=rmult)
        if rmax > r_ssh_p: r_ssh_p = r_ssh_p + rmult/2.
    
        rmin    = min( nmp.min(vs_m[:]) , nmp.min(vs_s[:])  )
        r_ssh_m = bt.round_to_multiple_of(rmin, prec=1, base=rmult)
        if rmin < r_ssh_m: r_ssh_m = r_ssh_m - rmult/2.
    
        fig = plt.figure(num = 1, figsize=(12,7.4), facecolor='w', edgecolor='k')
        ax1 = plt.axes([0.07, 0.22, 0.88, 0.73])
        ax1.set_xticks(vtime[::xticks_d])
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
        plt.xticks(rotation='60')
        plt.plot(vtime, vtime*0., '-', color='k', linewidth=2, label=None)
        plt.plot(vtime, vs_s[:], '-o', color=clr_sat, linewidth=2, label=cn_sat+' ("'+cv_sat+'")', zorder=10)
        plt.plot(vtime, vs_m[:], '-o', color=clr_mod, linewidth=2,  label=cn_mod+' ("'+cv_mod+'")', zorder=15)
        plt.yticks( nmp.arange(r_ssh_m, r_ssh_p+0.1, 0.1) )
        ax1.set_ylim(r_ssh_m*1.05,r_ssh_p*1.05)
        plt.ylabel('SSH [m]')
        ax1.set_xlim(vtime[0],vtime[Nsp-1])
        ax1.grid(color='k', linestyle='-', linewidth=0.3)
        plt.legend(loc="best", ncol=1, shadow=True, fancybox=True)
        cstart = str(round(bt.long_to_m180_p180(vlon[it1]),2))+"$^{\circ}$E, "+str(round(vlat[it1],2))+"$^{\circ}$N"
        cstop  = str(round(bt.long_to_m180_p180(vlon[it2]),2))+"$^{\circ}$E, "+str(round(vlat[it2],2))+"$^{\circ}$N"
        plt.title(r""+cn_box+": "+cstart+"  $\longrightarrow$ "+cstop)
        plt.savefig(cfigure, dpi=120, transparent=False)
        plt.close(1)


    if l_plot_spct:
        # Figure that shows all the spectrum for each segment:
        cfigure = dir_figs+'/'+cn_box+'_'+cseas+'_'+cn_mod+'--'+cn_sat+'_seg'+cseg+'.'+fig_ext
        ii = bp.plot("pow_spectrum_ssh")(vk[idx], x_all_spectra_s[jcpt,:], cfig_name=cfigure, \
                                         clab1=clabel_sat, cinfo=str(Nsp)+' points ('+str(nbp)+')', \
                                         L_min=13.5, L_max=1400., P_min_y=pow10_min_y, P_max_y=pow10_max_y, \
                                         vk2=vk[idx], vps2=x_all_spectra_m[jcpt,:], clab2=clabel_mod)
        

# Plotting mean spectrum:
vps_mod = nmp.mean(x_all_spectra_m[:,:],axis=0)
vps_sat = nmp.mean(x_all_spectra_s[:,:],axis=0) 
#cinfrm = cn_box+', '+cseas+' '+cyear+'\n  '+str(nb_v_seg)+' tracks\n  '+str(Nsp)+' points\n  '+r'$\Delta$L: '+str(round(rdist_sample,1))+' km'
cinfrm = cseas+' '+cyear+'\n'+str(nb_v_seg)+' tracks\n'+str(Nsp)+' points\n'+r'$\Delta$L: '+str(round(rdist_sample,1))+' km'


# remove white noise at fine scale for satellite (instrument) [advice from Clement Ubelmann]
cxtr_noise=''
if l_rm_i_noise:    
    rwn = nmp.mean(vps_sat[Nsp-15:Nsp])
    vps_sat = vps_sat - rwn
    cxtr_noise='_denoised'

# Sample spacing rdist_sample
cpout   = dir_figs+'/'+cn_box+'_MEAN_'+cn_mod+'--'+cn_sat+'__'+cseas+cextra+'___pow-spectrum'+cxtr_noise
cfigure = cpout+'.'+fig_ext


print(' *** cn_sat =', cn_sat)

ii = bp.plot("pow_spectrum_ssh")(vk[idx], vps_mod, clab1=clabel_mod, clr1=clr_mod, lw1=5, \
                                 cfig_name=cfigure, cinfo=cinfrm, logo_on=False, \
                                 L_min=10., L_max=1200., P_min_y=pow10_min_y, P_max_y=pow10_max_y, \
                                 l_show_k4=False, l_show_k5=True, l_show_k11o3=False, l_show_k2=True, \
                                 vk2=vk[idx], vps2=vps_sat, clab2=clabel_sat, clr2=clr_sat, lw2=4)


nmp.savez( cpout+'_mod.npz', vk[idx], vps_mod )
nmp.savez( cpout+'_sat'+cxtr_noise+'.npz', vk[idx], vps_sat )

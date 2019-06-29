#!/usr/bin/env python
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

#    L. Brodeau, 2018

import sys
from os import path as path
from string import replace
import numpy as nmp
import scipy.signal as signal
from netCDF4 import Dataset
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from string import find
import warnings
warnings.filterwarnings("ignore")
import time

import barakuda_plot as bp
import barakuda_tool as bt

reload(sys)
sys.setdefaultencoding('utf8')


cmodel = "eNATL60"

l_add_AJ_d = False
l_add_AJ_h = False

l_tapper      = True
l_detrend_lin = True ; # apply a linear detrending on data segment before computing spectrum...
l_rm_inst_noise = False

# Plots for each segment (only to debug!)
l_plot_rawd = False
l_plot_maps = False
l_plot_trck = False ; # plots individual model and satellite tracks...
l_plot_spct = False

r2Pi = 2.*nmp.pi

#dir_figs='/home/laurent/tmp/figs'
#dir_figs='../tmp'
dir_figs='./figs'
fig_ext='svg'
#fig_ext='png'


rcut_by_time = 1.2 # specify in seconds the time gap between two obs that means a discontinuity and therefore
#                  # should be considered as a cut!

rcut_by_dist = 7.8 # same for distance (in km) between two consecutive points!
#                  # time criterion "rcut_by_time" would have been sufficient if SARAL data was not bugged!!!
#                  # => in SARAL data, spotted two consecutive points with the usual dt and a huge dL !!!
#                  #   => like if the satellite undergone an extremely short huge acceleration !!!
#                  #   => ex: 3rd of August 2016 around 07:53:43 !!!

nvalid_seg = 120  # specify the minimum number of values a segment should contain to be considered and used!

dir_npz = './npz_AJ'

r_max_amp_ssh = 1.5 # in meters

clr_red = '#AD0000'
#clr_sat = '#091459'
#clr_mod = '#ED7C4C'

# Before:
#clr_mod = '#2C558A'
#if cn_sat_short=='SARAL':     clr_sat = '#ED7C4C' # (orange)
#if cn_sat_short=='Sentinel3': clr_sat = '#3DE079' # (green)
# ON:
clr_sat = '#ffed00'
clr_mod = '#008ab8'

rDPI=100.



jt1=0 ; jt2=0

narg = len(sys.argv)

if narg != 4 and narg != 6 :
    print 'Usage: '+sys.argv[0]+' <file> <variable model> <variable satel> (<jt1> <jt2>)'; sys.exit(0)

cf_in = sys.argv[1] ; cv_mdl=sys.argv[2] ; cv_eph=sys.argv[3]
if narg == 6:
    jt1=int(sys.argv[4]) ; jt2=int(sys.argv[5])

bt.chck4f(cf_in)

cfs  = path.basename(cf_in) ; nfs = len(cfs)
idx1 = nfs - find(cfs[nfs+1:0:-1],'_')
idx2 = nfs - find(cfs[nfs+1:0:-1],'.')-1
cn_satellite = cfs[idx1:idx2]
print '\n *** Satellite: '+cn_satellite

cn_sat_short = cn_satellite
if cn_satellite == 'SARAL-Altika': cn_sat_short='SARAL'

idx1 = find(cfs,'__')+2
idx2 = find(cfs,cmodel)
cn_box = cfs[idx1:idx2-1]
print ' *** Box: '+cn_box
idx1 = idx2
idx2 = find(cfs[idx2:],'_')+idx2
cn_model = cfs[idx1:idx2]
print ' *** Model: '+cn_model

cseas = ''
if ('JFM' in cfs) and not('JAS' in cfs) : cseas = 'JFM'; vseas = ['01','02','03']
if ('JAS' in cfs) and not('JFM' in cfs) : cseas = 'JAS'; vseas = ['07','08','09']
print ' *** Season: '+cseas+'\n'




cextra = ''
if l_add_AJ_d:
    # reading AJ file for current box (Spectrum from DAILY SSH!):
    cbox = str(int(cn_box[4:6]))
    print ' Reading AJ spectra for cbox = ', cbox
    jc=0
    for cm in vseas:
        cf_AJ = dir_npz+'/WaveSpec_Box_'+cbox+'_ssh_2013_'+cm+'.npz'
        print ' *** loading file '+cf_AJ
        data = nmp.load(cf_AJ)
        vkspec = data['kspec']
        vpspec = data['pspec']
        if jc==0:
            vkspec_3 = nmp.zeros((3,len(vkspec)))
            vpspec_3 = nmp.zeros((3,len(vkspec)))
        vkspec_3[jc,:] = vkspec[:]*1000. #*600. # he uses [m], I use [km]
        vpspec_3[jc,:] = vpspec[:]/1000. #/100.  #lolo: why sqrt ????
        jc=jc+1
    vk_AJ = nmp.mean(vkspec_3,axis=0) ; # Seasonal mean
    vp_AJ = nmp.mean(vpspec_3,axis=0) ; #      "
    cextra = '__AJ_d'




if l_add_AJ_h:
    # reading AJ file for current box (Spectrum from HOURLY SSH!):
    print ' Reading AJ spectra for cbox = ', cn_box
    cf_AJ = dir_npz+'/'+cseas+'_SSH_Spectra_hourly_dataset_'+cn_box+'.npz'
    print ' *** loading file '+cf_AJ
    data = nmp.load(cf_AJ)
    vk_AJ = data['kspec']*1000.
    vp_AJ = data['pspec']/1000.
    cextra = '__AJ_h'
    

id_in    = Dataset(cf_in)
vt_epoch = id_in.variables['time'][:]
vmodel   = id_in.variables[cv_mdl][:]
vsatel   = id_in.variables[cv_eph][:]
vlon     = id_in.variables['longitude'][:]
vlat     = id_in.variables['latitude'][:]
vdist    = id_in.variables['distance'][:]
id_in.close()
print "  => READ!\n"


# Debug: wave-signal with a wave-length of 150 km and 25 km (vdist in km)
#vmodel = 0.25*nmp.cos(vdist*r2Pi/150.) + 0.75*nmp.cos(vdist*r2Pi/25.)


nbr = len(vt_epoch)

# Initial raw plot:
# Create Matplotlib time array:


if l_plot_rawd or l_plot_trck:
    import matplotlib.dates as mdates


cyear = time.strftime("%Y", time.localtime(vt_epoch[2]))




ii=nbr/300
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
    plt.plot(vtime, vsatel, '-', color=clr_sat, linewidth=2, label=cn_sat_short+' ("'+cv_eph+'")', zorder=10)
    plt.plot(vtime, vmodel, '-', color=clr_mod, linewidth=2,  label=cn_model+' ("'+cv_mdl+'")', zorder=15)
    ax1.set_ylim(-r_max_amp_ssh,r_max_amp_ssh) ; ax1.set_xlim(vtime[0],vtime[nbr-1])
    plt.xlabel('Time [seconds since 1970]')
    plt.ylabel('SSH [m]')
    ax1.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(bbox_to_anchor=(0.55, 0.98), ncol=1, shadow=True, fancybox=True)
    plt.savefig(dir_figs+'/fig_raw_data_'+cn_model+'--'+cn_sat_short+'.'+fig_ext, dpi=120, transparent=True)
    plt.close(1)






vmask = vmodel.mask

(idx_ok,) = nmp.where(vmask==False) # indexes with valid values!


nbr_v = len(idx_ok)

print ' *** '+str(nbr_v)+' valid points out of '+str(nbr)+' !'

# Will extract the N valid data segments:
nb_seg=0
idx_seg_start = [] ; # index of first valid point of the segment
idx_seg_stop  = [] ; # index of last  valid point of the segment

jr=0
while jr < nbr:
    # Ignoring masked values and zeros...
    if (not vmask[jr]) and (vmodel[jr]!=0.0) and (vmodel[jr]<100.):
        nb_seg = nb_seg + 1
        print '\n --- found seg #'+str(nb_seg)+' !'
        idx_seg_start.append(jr)
        print ' => starting at jt='+str(jr)
        #while (not vmask[jr+1]) and (vmodel[jr+1]!=0.0) and (vt_epoch[jr+1]-vt_epoch[jr] < rcut_by_time) and (vmodel[jr]<100.) :
        while (not vmask[jr+1]) and (vmodel[jr+1]!=0.0) and (vt_epoch[jr+1]-vt_epoch[jr] < rcut_by_time) and (vdist[jr+1]-vdist[jr] < rcut_by_dist) and (vmodel[jr]<100.) :
            jr = jr+1
            if jr==nbr-1: break
        idx_seg_stop.append(jr)
        print ' => and stoping at jt='+str(jr)
    jr = jr+1

if len(idx_seg_start) != nb_seg: print ' ERROR #1!'; sys.exit(1)



# Maximum number of poins in the segments:

isd = nmp.asarray(idx_seg_stop[:]) - nmp.asarray(idx_seg_start[:]) + 1
print '\n lengths =>', isd[:]
nbp_max = max(isd)
print '\n *** Longest segments has nbp_max='+str(nbp_max)+' points!\n'

# Treat only segments with at least nvalid_seg points:
(vtreat,) = nmp.where(isd >= nvalid_seg)
nb_v_seg = len(vtreat)
rN = nmp.mean(isd[vtreat])
print '\nMean segment-length for the '+str(nb_v_seg)+' segments with at least '+str(nvalid_seg)+' points:', round(rN,1)

Nsp = int(rN/10.)*10
print '  ==> Nsp = '+str(Nsp)
(vtreat,) = nmp.where(isd >= Nsp)
nb_v_seg = len(vtreat)
print ' ==> will use '+str(nb_v_seg)+' segments with a fixed length of '+str(Nsp)+' point!\n'


#print ' *** will treat '+str(nb_v_seg)+' segments out of '+str(nb_seg)+' (need at least '+str(nvalid_seg)+' points)'


x_all_spectra_s = nmp.zeros((nb_v_seg,Nsp)) # array to store all spectra in...
x_all_spectra_m = nmp.zeros((nb_v_seg,Nsp)) # array to store all spectra in...

clabel_sat=cn_sat_short+' ("'+cv_eph+'")'
#clabel_mod=cn_model+'  ("'+cv_mdl+'"), 1D along-track'
clabel_mod=cn_model+' (1D along-track)'
#clabel3=cn_model+' (2D box)'






jcpt = -1
for js in vtreat:
    jcpt= jcpt+1
    it1 = idx_seg_start[js]
    it2 = idx_seg_stop[js]
    nbp = it2-it1+1    
    cseg = '%2.2i'%(js+1)
    print '\n\n ###################################'
    print '  *** Seg #'+cseg+' of '+cn_box+':'
    print '  ***   => originally '+str(nbp)+' points in this segment (from '+str(it1)+' to '+str(it2)+')'

    # nb of points in excess / Nsp:
    nxcs = nbp - Nsp
    jmp_strt = nxcs/2
    jmp_stop = nxcs/2 + nxcs%2
    it1 = it1+jmp_strt
    it2 = it2-jmp_stop
    print '  ***   => we only retain '+str(it2-it1+1)+' points from '+str(it1)+' to '+str(it2)+'!'
    
    # Checking the typical distance (in km) between two measures:
    dmean = nmp.mean(vdist[it1+1:it2+1]-vdist[it1:it2])
    print '\n Mean distance between two consecutive points is '+str(dmean)+' km\n'
    # Sample spacing in [km] (inverse of the sampling rate):
    rdist_sample = round(dmean,3)
    print ' => will use a spatial sample spacing of '+str(rdist_sample)+' km\n'

    if l_plot_trck:
        # Create Matplotlib time array:
        vtime = nmp.zeros(Nsp)
        for jt in range(Nsp): vtime[jt] = mdates.epoch2num(vt_epoch[it1+jt])

    # First centering the two time-series about 0, and tappering at extremities (filling with zeros)
    vs_s = nmp.zeros(Nsp) ; vs_m = nmp.zeros(Nsp)
    print '             (length = '+str(len(vsatel[it1:it2+1]))+' / '+str(Nsp)+')'
    
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
        plt.savefig(dir_figs+'/'+cn_model+'--'+cn_sat_short+'_seg'+cseg+'_map.'+fig_ext, dpi=120, transparent=True)
        plt.close(1)


    if l_plot_trck:
        cfigure = dir_figs+'/'+cn_box+'_'+cseas+'_'+cn_model+'--'+cn_sat_short+'_seg'+cseg+'_track.'+fig_ext
        ii=Nsp/300
        ib=max(ii-ii%10,1)
        xticks_d=30.*ib
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
        plt.plot(vtime, vs_s[:], '-o', color=clr_sat, linewidth=2, label=cn_sat_short+' ("'+cv_eph+'")', zorder=10)
        plt.plot(vtime, vs_m[:], '-o', color=clr_mod, linewidth=2,  label=cn_model+' ("'+cv_mdl+'")', zorder=15)
        plt.yticks( nmp.arange(r_ssh_m, r_ssh_p+0.1, 0.1) )
        ax1.set_ylim(r_ssh_m*1.05,r_ssh_p*1.05)
        plt.ylabel('SSH [m]')
        ax1.set_xlim(vtime[0],vtime[Nsp-1])
        ax1.grid(color='k', linestyle='-', linewidth=0.3)
        plt.legend(loc="best", ncol=1, shadow=True, fancybox=True)
        cstart = str(round(bt.lon_180_180(vlon[it1]),2))+"$^{\circ}$E, "+str(round(vlat[it1],2))+"$^{\circ}$N"
        cstop  = str(round(bt.lon_180_180(vlon[it2]),2))+"$^{\circ}$E, "+str(round(vlat[it2],2))+"$^{\circ}$N"
        plt.title(r""+cn_box+": "+cstart+"  $\longrightarrow$ "+cstop)
        plt.savefig(cfigure, dpi=120, transparent=False)
        plt.close(1)


    if l_plot_spct:
        # Figure that shows all the spectrum for each segment:
        cfigure = dir_figs+'/'+cn_box+'_'+cseas+'_'+cn_model+'--'+cn_sat_short+'_seg'+cseg+'.'+fig_ext
        ii = bp.plot("pow_spectrum_ssh")(vk[idx], x_all_spectra_s[jcpt,:], cfig_name=cfigure, \
                                         clab1=clabel_sat, cinfo=str(Nsp)+' points ('+str(nbp)+')', \
                                         L_min=13.5, L_max=1400., P_min_y=-8, P_max_y=2,    \
                                         vk2=vk[idx], vps2=x_all_spectra_m[jcpt,:], clab2=clabel_mod)
        

# Plotting mean spectrum:
vps_mod = nmp.mean(x_all_spectra_m[:,:],axis=0)
vps_sat = nmp.mean(x_all_spectra_s[:,:],axis=0) 
#cinfrm = cn_box+', '+cseas+' '+cyear+'\n  '+str(nb_v_seg)+' tracks\n  '+str(Nsp)+' points\n  '+r'$\Delta$L: '+str(round(rdist_sample,1))+' km'
cinfrm = cseas+' '+cyear+'\n'+str(nb_v_seg)+' tracks\n'+str(Nsp)+' points\n'+r'$\Delta$L: '+str(round(rdist_sample,1))+' km'

# Sample spacing rdist_sample
cpout   = dir_figs+'/'+cn_box+'_MEAN_'+cn_model+'--'+cn_sat_short+'__'+cseas+cextra+'___pow-spectrum'
cfigure = cpout+'.'+fig_ext


# remove white noise at fine scale for satellite (instrument) [advice from Clement Ubelmann]
if l_rm_inst_noise:    
    rwn = nmp.mean(vps_sat[Nsp-15:Nsp])
    vps_sat = vps_sat - rwn



print ' *** cn_sat_short =', cn_sat_short    

ii = bp.plot("pow_spectrum_ssh")(vk[idx], vps_mod, clab1=clabel_mod, clr1=clr_mod, lw1=5, \
                                 cfig_name=cfigure, cinfo=cinfrm, logo_on=False, \
                                 L_min=10., L_max=1200., P_min_y=-6, P_max_y=1, \
                                 l_show_k4=False, l_show_k5=True, l_show_k11o3=False, l_show_k2=True, \
                                 vk2=vk[idx], vps2=vps_sat, clab2=clabel_sat, clr2=clr_sat, lw2=4)


nmp.savez( cpout+'_mod.npz', vk[idx], vps_mod )
nmp.savez( cpout+'_sat.npz', vk[idx], vps_sat )

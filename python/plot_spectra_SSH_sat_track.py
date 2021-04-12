#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

#    L. Brodeau, 2018

from sys import exit
from os import path, mkdir
import argparse as ap
import numpy as nmp
from netCDF4 import Dataset
import time
import matplotlib.dates as mdates
#
import climporn as cp
import gonzag   as gzg ; # Gonzag package => https://github.com/brodeau/gonzag

ivrb = 1

l_rm_i_noise = False

l_plot_rawd = True ; # plot the whole SSH times series of model and satellite 
# Plots for each segment (only to debug!):
l_plot_trck = True ; # plots model and satellite tracks for each segment
l_plot_spct = False ; # plot spectra for each segment

dir_figs='./figs'
fig_ext='svg'
#fig_ext='png'


# Because track files downloaded are usually full of inconsistencies / time and distance done
# between 2 consecutive points, we need to spot these inconsitencies with the 2 following criteria:
# ex: normally SARAL takes one measure every 1 s, => ds~=1 s, time during which it should move by
#     roughly dx ~= 7 km on the ground...
#     However, the netCDF files downloaded are full of points with dt=2s and the corresponding dx is
#     still dx ~= 7 km, this is clearly a bug !!!

rcut_by_time = 1.2 # specify in seconds the time gap between two obs to detect and discontinuity and therefore
#                  # should be considered as a cut!

rcut_by_dist = 7.8 # same as rcut_by_time, but in terms of distance (in km) between two consecutive points!
#                  # time criterion "rcut_by_time" would have been sufficient if SARAL data was not bugged!!!
#                  # => in SARAL data, spotted two consecutive points with the usual dt and a huge dL !!!
#                  #   => like if the satellite undergone an extremely short huge acceleration !!!
#                  #   => ex: 3rd of August 2016 around 07:53:43 !!!

nlen_valid_seg_default = 120  # specify the minimum number of points a segment should contain to be considered and used!


# Look and feel for the plot:
clr_sat = '#AD0000'
clr_mod = '#008ab8'

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
parser.add_argument('-n', '--npseg', type=int, default=nlen_valid_seg_default, help='minimum number of points for a valid segment (default: '+str(nlen_valid_seg_default)+')')
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
nlen_valid_seg = args.npseg
cn_box = args.name_box
cn_mod = args.name_mod
cn_sat = args.name_sat
pow10_min_y = args.pmin_y
pow10_max_y = args.pmax_y


cfs  = path.basename(cf_in)

cseas = ''
if ('JFM' in cfs) and not('JAS' in cfs) : cseas = 'JFM'; vseas = ['01','02','03']
if ('JAS' in cfs) and not('JFM' in cfs) : cseas = 'JAS'; vseas = ['07','08','09']
if cseas != '': print('\n *** Season: '+cseas+'\n')



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


nbr = len(vt_epoch)

# Once for all, convertion from Epoch Unix time to "day since 0001":
VT = mdates.epoch2num(vt_epoch)


cyear = time.strftime("%Y", time.localtime(vt_epoch[2]))

print(' *** Year = '+cyear+'\n')


# Plot time-series of SSH as found in input file:
if l_plot_rawd:
    cfigure = dir_figs+'/fig_raw_data_'+cn_mod+'--'+cn_sat+'.'+fig_ext
    ii = cp.PlotInputSeries(VT, vsatel, vmodel, cfigure, \
                            clabS=cn_sat+' ("'+cv_sat+'")', clabM=cn_mod+' ("'+cv_mod+'")')


vmask = vmodel.mask
(idx_ok,) = nmp.where(vmask==False) # indexes with valid values!
nbr_v = len(idx_ok)
print(' *** Input data contains '+str(nbr_v)+' non masked data points among the '+str(nbr)+' points...\n')


# Extract the Ns continuous data segments:
ISeg_start, ISeg_stop = gzg.FindUnbrokenSegments( vt_epoch, vdist, vmodel, vmask, rcut_time=rcut_by_time, rcut_dist=rcut_by_dist )

# Selecting proper segments:
NbSeg, Nsl, IDEDSeg = gzg.SegmentSelection(ISeg_start, ISeg_stop, np_valid_seg=nlen_valid_seg)
# validity check:
#for js in range(NbSeg):
#    print(' * Seg # ',js+1,' => it1, it2 =', IDEDSeg[js,:], ' ==> len = ', IDEDSeg[js,1]-IDEDSeg[js,0]+1)
#print(' Nsl = ',Nsl)

# Process data on segment so ready for FFT:
XPs, XPm, rdist_sample = gzg.Process4FFT( IDEDSeg, vdist, vsatel, vmodel )

# Apply FFT !
Kwn, PwSpc_s, PwSpc_m = gzg.ApplyFFT( IDEDSeg, XPs, XPm, rdist_sample )

# Everything is computed, time for control plots:
clbl_sat = cn_sat+' ("'+cv_sat+'")'
clbl_mod = cn_mod+' ("'+cv_mod+'")'
#
for js in range(NbSeg):
    it1 = IDEDSeg[js,0]
    it2 = IDEDSeg[js,1]
    if l_plot_trck:
        # Figure that shows processed time-series of SSH for model and satellite on each segment:
        cfigure = dir_figs+'/'+cn_box+'_'+cseas+'_'+cn_mod+'--'+cn_sat+'_seg'+'%2.2i'%(js+1)+'_track.'+fig_ext
        cstart = str(round(cp.degE_to_degWE(vlon[it1]),2))+"$^{\circ}$E, "+str(round(vlat[it1],2))+"$^{\circ}$N"
        cstop  = str(round(cp.degE_to_degWE(vlon[it2]),2))+"$^{\circ}$E, "+str(round(vlat[it2],2))+"$^{\circ}$N"
        #
        ii = cp.PlotSegmentTrack(VT[it1:it2+1], XPs[js,:], XPm[js,:], cfig=cfigure, \
                              ctitle=r""+cn_box+": "+cstart+"  $\longrightarrow$ "+cstop, \
                              clabS=clbl_sat, clabM=clbl_mod)
    #
    if l_plot_spct:
        # Figure that shows all the spectrum for each segment:
        cfigure = dir_figs+'/'+cn_box+'_'+cseas+'_'+cn_mod+'--'+cn_sat+'_seg'+'%2.2i'%(js+1)+'.'+fig_ext
        ii = cp.plot("pow_spectrum_ssh")(Kwn, PwSpc_s[js,:], cfig_name=cfigure, \
                                         clab1=clbl_sat, cinfo=str(Nsl)+' points ('+str(it2-it1+1)+')', \
                                         L_min=13.5, L_max=1400., P_min_y=pow10_min_y, P_max_y=pow10_max_y, \
                                         vk2=Kwn, vps2=PwSpc_m[js,:], clab2=clbl_mod)
        

# Plotting mean spectrum:
vps_mod = nmp.mean(PwSpc_m[:,:],axis=0)
vps_sat = nmp.mean(PwSpc_s[:,:],axis=0) 

ctime = ''
if cseas != '': ctime=cseas+' '+cyear+'\n'
cinfrm = ctime+str(NbSeg)+' segments\n'+str(Nsl)+' points/segment\n'+r'$\Delta$d sat.: '+str(round(rdist_sample,1))+' km'


# remove white noise at fine scale for satellite (instrument) [advice from Clement Ubelmann]
cxtr_noise=''
if l_rm_i_noise:    
    rwn = nmp.mean(vps_sat[Nsl-15:Nsl])
    vps_sat = vps_sat - rwn
    cxtr_noise='_denoised'

# Sample spacing rdist_sample
cpout   = dir_figs+'/'+cn_box+'_MEAN_'+cn_mod+'--'+cn_sat+'__'+cseas+'___pow-spectrum'+cxtr_noise
cfigure = cpout+'.'+fig_ext

if ivrb>1: print(' *** cn_sat =', cn_sat)

ii = cp.plot("pow_spectrum_ssh")(Kwn, vps_mod, clab1=clbl_mod, clr1=clr_mod, lw1=5, \
                                 cfig_name=cfigure, cinfo=cinfrm, logo_on=False, \
                                 L_min=10., L_max=1200., P_min_y=pow10_min_y, P_max_y=pow10_max_y, \
                                 l_show_k4=False, l_show_k5=True, l_show_k11o3=False, l_show_k2=True, \
                                 vk2=Kwn, vps2=vps_sat, clab2=clbl_sat, clr2=clr_sat, lw2=4)

# Saved results:
nmp.savez( cpout+'_mod.npz', Kwn, vps_mod )
nmp.savez( cpout+'_sat'+cxtr_noise+'.npz', Kwn, vps_sat )

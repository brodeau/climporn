#!/usr/bin/env python3

# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

# L. Brodeau, 2022

import sys
import os
import numpy as nmp
from netCDF4 import Dataset
#import string

# Extra
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import climporn as cp

DPIsvg=200

def init_fig( font_rat=1, color_top='k' ):
    rr = font_rat
    params = { 'font.family':'Open Sans', 'font.weight':    'normal',
               'font.size':       int(12.*rr),
               'legend.fontsize': int(22.*rr),
               'xtick.labelsize': int(13.*rr), 'ytick.labelsize': int(12.*rr),
               'axes.labelsize':  int(14.*rr), # label of colorbar
               'text.color': color_top, 'axes.labelcolor': color_top, 'xtick.color':color_top, 'xtick.color':color_top }
    mpl.rcParams.update(params)
    return 0


################################################################################################

if __name__ == '__main__':

    print('\n')

    narg = len(sys.argv)
    if narg != 4:
        print('Usage: '+sys.argv[0]+' <nemo_coordinates_file> <CONF_NAME> <Nsubsamp>'); sys.exit(0)

    cf_mm    =     sys.argv[1]
    conf     =     sys.argv[2]
    isubsamp = int(sys.argv[3])

    #
    cp.chck4f(cf_mm)

    id_mm = Dataset(cf_mm)

    list_dim = id_mm.dimensions.keys()
    l_td = 'time_counter' in list_dim or 'time' in list_dim
    l_zd = 'z' in list_dim

    nbdim = 2
    if l_td: nbdim = nbdim + 1
    if l_zd: nbdim = nbdim + 1
    
    print(list_dim, l_td, ltz)
    
    Ni = id_mm.dimensions['x'].size
    Nj = id_mm.dimensions['y'].size

    print('\n *** Horizontal domain => '+str(Ni)+' x '+str(Nj))

    list_var = id_mm.variables.keys()
    print('\n *** Variables in '+cf_mm+':\n', list(list_var), '\n')
    list_treat = []
    list_ndims = []
    list_point = []
    for cv in list_var:
        if not cv in ['nav_lon','nav_lat','deptht','time_counter','mask']:
            nbd = len(id_mm.variables[cv].dimensions)
            if nbd > 1:
                list_treat.append(cv)
                list_ndims.append(nbd)
                cc=cv[len(cv)-1]
                if not cc in ['t','u','v','f']: cc = 't'
                list_point.append(cc)
    nbv = len(list_treat)
    xcv_n = nmp.asarray(list_treat)
    xcv_d = nmp.asarray(list_ndims)       
    xcv_p = nmp.asarray(list_point)


    XVF = nmp.zeros((nbv,Nj,Ni))
    jv  = 0
    for cv in xcv_n:

        print('\n  ==> reading variable # '+str(jv+1)+' :'+cv+' ('+str(xcv_d[jv])+'D), grid point = '+xcv_p[jv])

        if cv == 'glamf': id_glamf=jv
        if cv == 'gphif': id_gphif=jv
        if cv == 'glamt': id_glamt=jv
        if cv == 'gphit': id_gphit=jv
        if cv == 'e1t'  : id_e1t=jv
        if cv == 'e2t'  : id_e2t=jv

        # There is always time_counter as a dimmension
        if   xcv_d[jv] == nbdim:
            if l_td:
                XVF[jv,:,:] = id_mm.variables[cv][0,:,:]
            else:
                XVF[jv,:,:] = id_mm.variables[cv][:,:]
        else:
            print(' PROBLEM: variable '+cv+' has '+str(xcv_d[jv])+' dimensions!!!'); sys.exit(0)
           
        jv=jv+1

    id_mm.close()
    
    

    cf_out = 'grid_nextsimDG_2x'+conf+'.nc'
        


    # Location of north pole:
    

    NNN = nmp.argmax( XVF[id_gphit,:,:] )
    NjNP = NNN // Ni
    NiNP = NNN % Ni    
    print( "\n *** North Pole found at : Ni, Nj = ", NiNP, NjNP  )    
    #cp.dump_2d_field( 'LAT.nc',  XVF[id_gphit,:,:] ) ; #debug

    #########################################
    #             FIGURES
    #########################################

    # Zoom around North Pole:
    i1=NiNP-15 ; i2=NiNP+15
    j1=NjNP-15 ; j2=NjNP+15
    

    # Small box including the North Pole to spot possible inconsistencies in the grid:
    ii = cp.PlotGridGlobe( XVF[id_glamf,j1:j2,i1:i2], XVF[id_gphif,j1:j2,i1:i2], chemi='N', lon0=-35., cfig_name='zoom_NH_35W_f_OUT_ortho_WHITE.svg',
                           nsubsamp=1, rdpi=DPIsvg, ldark=False, nzoom=5 )

    # White:
    ii = cp.PlotGridGlobe( XVF[id_glamf,:,:], XVF[id_gphif,:,:], chemi='N', lon0=-35., cfig_name='NH_35W_f_OUT_ortho_WHITE.svg',  nsubsamp=isubsamp, rdpi=DPIsvg, ldark=False )
    ii = cp.PlotGridGlobe( XVF[id_glamf,:,:], XVF[id_gphif,:,:], chemi='N', lon0=130., cfig_name='NH_130E_f_OUT_ortho_WHITE.svg', nsubsamp=isubsamp, rdpi=DPIsvg, ldark=False )



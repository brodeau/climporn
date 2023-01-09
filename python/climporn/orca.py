#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
#   ///// https://github.com/brodeau/climporn \\\\\
#
#       L. Brodeau, 2021
#
############################################################################

import sys
import numpy as np

#TRG_NAME_CONF = 'HEMI'

#DPIsvg=200
#DPIpng=200



def get_basin_info( cf_bm ):    
    from netCDF4 import Dataset
    l_b_names = [] ; l_b_lgnms = []
    l_b_names.append(u'GLO') ; l_b_lgnms.append(u'Global Ocean')
    id_bm = Dataset(cf_bm)
    list_var = id_bm.variables.keys()
    for cv in list_var:
        if cv[:5] == 'tmask':
            l_b_names.append(cv[5:])
            l_b_lgnms.append(id_bm.variables[cv].long_name)
    id_bm.close()
    return l_b_names, l_b_lgnms


def lon_reorg_orca(ZZ, Xlong, ilon_ext=0, v_jx_jump_p=170., v_jx_jump_m=-170.):
    #
    #
    # IN:
    # ===
    # ZZ       : 1D, 2D, or 3D array to treat
    # Xlong    : 1D or 2D array containing the longitude, must be consistent with the shape of ZZ !
    # ilon_ext : longitude extention of the map you want (in degrees)
    #
    # OUT:
    # ====
    # ZZx     : re-organized array, mind the possibility of modified x dimension !
    #
    import climporn.utils as cpu
    #
    idim_lon = len(np.shape(Xlong))
    if idim_lon not in [ 1 , 2 ]:
        print('util_orca.lon_reorg_orca: ERROR => longitude array "Xlong" must be 1D or 2D !'); sys.exit(0)
    #
    if idim_lon == 2: (nj,ni) = np.shape(Xlong)
    if idim_lon == 1:      ni = len(Xlong)
    #
    vlon = np.zeros(ni)
    #
    if idim_lon == 2: vlon[:] = Xlong[nj/3,:]
    if idim_lon == 1: vlon[:] = Xlong[:]
    #
    lfound_jx_jump = False
    ji=0
    while ( not lfound_jx_jump and ji < ni-1):
        if vlon[ji] > v_jx_jump_p and vlon[ji+1] < v_jx_jump_m:
            jx_jump = ji + 1
            lfound_jx_jump = True
        ji = ji + 1
    print('  *** clprn_orca.lon_reorg_orca >> Longitude jump at ji = ', jx_jump)
    #
    lfound_jx_zero = False
    ji=0
    while ( not lfound_jx_zero and ji < ni-1):
        if vlon[ji] < 0. and vlon[ji+1] > 0.:
            jx_zero = ji + 1
            lfound_jx_zero = True
        ji = ji + 1
    print('  *** clprn_orca.lon_reorg_orca >> Longitude zero at ji = ', jx_zero)
    #
    del vlon
    
    jx_oo = 2  # orca longitude overlap...
    vdim = ZZ.shape
    ndim = len(vdim)

    if ndim < 1 or ndim > 4:
        print('util_orca.lon_reorg_orca: ERROR we only treat 1D, 2D, 3D or 4D arrays...')

    if ndim == 4:
        #
        [ nr, nz , ny , nx ] = vdim ;     nx0 = nx - jx_oo
        ZZx  = np.zeros((nr, nz, ny, nx0))
        ZZx_ext  = np.zeros((nr, nz, ny, (nx0+ilon_ext)))
        #
        for jx in range(jx_zero,nx):
            ZZx[:,:,:,jx-jx_zero] = ZZ[:,:,:,jx]
        for jx in range(jx_oo,jx_zero):
            ZZx[:,:,:,jx+(nx-jx_zero)-jx_oo] = ZZ[:,:,:,jx]
        #
        if ilon_ext == 0: ZZx_ext[:,:,:,:] = ZZx[:,:,:,:]
    #
    #
    if ndim == 3:
        #
        [ nz , ny , nx ] = vdim ;     nx0 = nx - jx_oo
        #print('nx, ny, nz = ', nx, ny, nz
        #
        ZZx  = np.zeros(nx0*ny*nz) ;  ZZx.shape = [nz, ny, nx0]
        ZZx_ext  = np.zeros((nx0+ilon_ext)*ny*nz) ;  ZZx_ext.shape = [nz, ny, (nx0+ilon_ext)]
        #
        for jx in range(jx_zero,nx):
            ZZx[:,:,jx-jx_zero] = ZZ[:,:,jx]
        for jx in range(jx_oo,jx_zero):
            ZZx[:,:,jx+(nx-jx_zero)-jx_oo] = ZZ[:,:,jx]
        #
        if ilon_ext == 0: ZZx_ext[:,:,:] = ZZx[:,:,:]
    #
    #
    if ndim == 2:
        #
        [ ny , nx ] = vdim ;     nx0 = nx - jx_oo
        #print('nx, ny = ', nx, ny
        #
        ZZx  = np.zeros(nx0*ny) ;  ZZx.shape = [ny, nx0]
        ZZx_ext  = np.zeros((nx0+ilon_ext)*ny) ;  ZZx_ext.shape = [ny, (nx0+ilon_ext)]
        #
        for jx in range(jx_zero,nx):
            ZZx[:,jx-jx_zero] = ZZ[:,jx]
        for jx in range(jx_oo,jx_zero):
            ZZx[:,jx+(nx-jx_zero)-jx_oo] = ZZ[:,jx]
        #
        if ilon_ext == 0: ZZx_ext[:,:] = ZZx[:,:]
    #
    #
    if ndim == 1:
        #
        [ nx ] = vdim ;     nx0 = nx - jx_oo
        #print('nx = ', nx
        #
        ZZx  = np.zeros(nx0) ;  ZZx.shape = [nx0]
        ZZx_ext  = np.zeros(nx0+ilon_ext) ;  ZZx_ext.shape = [nx0+ilon_ext]
        #
        for jx in range(jx_zero,nx):
            ZZx[jx-jx_zero] = ZZ[jx]
            #print(jx-jx_zero, 'prend', jx, '    ', vlon[jx]
            #
        #print(''
        for jx in range(jx_oo,jx_zero):
            ZZx[jx+(nx-jx_zero)-jx_oo] = ZZ[jx]
            #print(jx+(nx-jx_zero)-jx_oo, 'prend', jx, '    ', vlon[jx]
        #
        if ilon_ext == 0: ZZx_ext[:] = ZZx[:]
        #iwa = np.where(vlon0 < 0.) ; vlon0[iwa] = vlon0[iwa] + 360.
        #
        #
        #
        # Now longitude extenstion:
    if ilon_ext > 0: ZZx_ext = cpu.extend_domain(ZZx, ilon_ext)
    del ZZx

    return ZZx_ext








def conf_exp(ccexp):
    #
    # Find the CONF from CONF-EXP and exit if CONF does not exist!
    #
    i = 0 ; conf = ''
    while i < len(ccexp) and ccexp[i] != '-' : conf = conf+ccexp[i]; i=i+1
    #print('conf =', conf, '\n'
    return conf



def mean_3d(XD, LSM, XVOL):
    #
    # XD             : 3D+T array containing data
    # LSM            : 3D land sea mask
    # XVOL           : 3D E1T*E2T*E3T  : 3D mesh volume
    #
    # RETURN vmean: vector containing 3d-averaged values at each time record

    ( lt, lz, ly, lx ) = np.shape(XD)

    if np.shape(LSM) != ( lz, ly, lx ):
        print('ERROR: mean_3d.clprn_orca.py => XD and LSM do not agree in shape!')
        sys.exit(0)
    if np.shape(XVOL) != ( lz, ly, lx ):
        print('ERROR: mean_3d.clprn_orca.py => XD and XVOL do not agree in shape!')
        sys.exit(0)

    vmean = np.zeros(lt)
    
    XX = LSM[:,:,:]*XVOL[:,:,:]
    rd = np.sum( XX )
    XX = XX/rd
    if rd > 0.:
        for jt in range(lt):
            vmean[jt] = np.sum( XD[jt,:,:,:]*XX )
    else:
        vmean[:] = np.nan

    return vmean


def mean_2d(XD, LSM, XAREA):
    #
    # XD             : 2D+T array containing data
    # LSM            : 2D land sea mask
    # XAREA          : 2D E1T*E2T  : mesh area
    #
    # RETURN vmean: vector containing 2d-averaged values at each time record

    ( lt, ly, lx ) = np.shape(XD)

    if np.shape(LSM) != ( ly, lx ):
        print('ERROR: mean_2d.clprn_orca.py => XD and LSM do not agree in shape!')
        sys.exit(0)
    if np.shape(XAREA) != ( ly, lx ):
        print('ERROR: mean_2d.clprn_orca.py => XD and XAREA do not agree in shape!')
        sys.exit(0)

    vmean = np.zeros(lt)
    XX = LSM[:,:]*XAREA[:,:]
    rd = np.sum( XX )

    # Sometimes LSM can be 0 everywhere! => rd == 0. !
    if rd > 0.:
        XX = XX/rd
        for jt in range(lt):
            vmean[jt] = np.sum( XD[jt,:,:]*XX )
    else:
        vmean[:] = np.nan

    return vmean






def ij_from_xy(xx, yy, xlon, xlat):
    #
    #=============================================================
    # Input:
    #       xx : longitude of searched point (float)
    #       yy : latitude  of searched point (float)
    #       xlon  : 2D array of longitudes of the ORCA grid
    #       xlat  : 2D array of latitudes  of the ORCA grid
    # Output:
    #       ( ji, jj ) : coresponding i and j indices on the ORCA grid    
    #=============================================================    
    #
    ji=0 ; jj=0
    #
    if xx < 0.: xx = xx + 360.
    #
    (nj , ni) = xlon.shape
    iwa = np.where(xlon < 0.) ; xlon[iwa] = xlon[iwa] + 360. # Want only positive values in longitude:
    #
    # Southernmost latitude of the ORCA domain:
    ji0 = np.argmin(xlat[0,:])
    lat_min = xlat[0,ji0]
    ji = ji0
    yy = max( yy, lat_min )
    #
    # Northernmost latitude of the ORCA domain:    
    ji0 = np.argmax(xlat[nj-1,:])
    lat_max = xlat[nj-1,ji0]
    yy = min( yy, lat_max )
    #
    A = np.abs( xlat[:-2,:-2]-yy ) + np.abs( xlon[:-2,:-2]-xx )
    (jj,ji) = np.unravel_index(A.argmin(), A.shape)
    #
    return ( ji, jj )



def transect_zon_or_med(x1, x2, y1, y2, xlon, xlat): #, rmin, rmax, dc):
    #
    #=============================================================
    # Input:
    #       x1,x2 : longitudes of point P1 and P2 defining the section (zonal OR meridional)
    #       y1,y2 : latitudes of point P1 and P2 defining the section (zonal OR meridional)
    #            => so yes! either x1==x2 or y1==y2 !
    #       xlon  : 2D array of longitudes of the ORCA grid
    #       xlat  : 2D array of latitudes  of the ORCA grid
    # Output:
    #       ( ji1, ji2, jj1, jj2 ) : coresponding i and j indices on the ORCA grid
    #=============================================================    
    #
    ji1=0 ; ji2=0 ; jj1=0 ; jj2=0
    lhori = False ; lvert = False
    #
    if y1 == y2: lhori = True
    if x1 == x2: lvert = True
    #
    if not (lhori or lvert) :
        print('transect_zon_or_med only supports horizontal or vertical sections!')
        sys.exit(0)
    #
    (ji1,jj1) = ij_from_xy(x1, y1, xlon, xlat)
    (ji2,jj2) = ij_from_xy(x2, y2, xlon, xlat)
    #
    if lhori and (jj1 != jj2): jj2=jj1
    if lvert and (ji1 != ji2): ji2=ji1
    #
    return ( ji1, ji2, jj1, jj2 )


def shrink_domain(LSM):
    # Decrasing the domain size to only retain the (rectangular region) with
    # ocean points (LSM == 1)
    #
    # Input:
    #     LSM : 2D land sea mask array
    # Output:
    #  (i1,j1,i2,j2): coordinates of the 2 points defining the rectangular region 
    #
    ( ly, lx ) = np.shape(LSM)
    #
    #if np.shape(LSM) != ( lz, ly, lx ):
    #    print('ERROR: shrink_domain.clprn_orca.py => XD and LSM do not agree in shape!'
    #    sys.exit(0)
    (vjj , vji)  = np.where(LSM[:,:]>0.5)
    j1 = max( np.min(vjj)-2 , 0    )
    i1 = max( np.min(vji)-2 , 0    )
    j2 = min( np.max(vjj)+2 , ly-1 ) + 1
    i2 = min( np.max(vji)+2 , lx-1 ) + 1
    #
    if (i1,j1,i2,j2) != (0,0,lx,ly):
        print('       ===> zooming on i1,j1 -> i2,j2:', i1,j1,'->',i2,j2)
    if (i1,i2) == (0,lx): i2 = i2-2 ; # Mind east-west periodicity overlap of 2 points...
    #
    return (i1,j1,i2,j2)




##################


#def init_fig( font_rat=1, color_top='k' ):
#    rr = font_rat
#    params = { 'font.family':'Open Sans', 'font.weight':    'normal',
#               'font.size':       int(12.*rr),
#               'legend.fontsize': int(22.*rr),
##               'xtick.labelsize': int(13.*rr), 'ytick.labelsize': int(12.*rr),
#               'axes.labelsize':  int(14.*rr), # label of colorbar
#               'text.color': color_top, 'axes.labelcolor': color_top, 'xtick.color':color_top, 'xtick.color':color_top }
#    mpl.rcParams.update(params)
#    return 0


def PlotGridGlobe( Xglamf, Xgphif, Xglamt=[], Xgphit=[], chemi='N', lon0=-35., lat0=45., cfig_name='mesh_globe_ortho.svg',
                   nsubsamp=5, rdpi=200, nxcld_n=0, hres=0.25, ldark=False, nzoom=1, linew=0.2,
                   lNPzoom=False ):
    """
            Shows the actual grid meshes on an orthographic projection of the globe
            * chemi     => which hemisphere to look at (N/S)
            * lon0      => facing longitude [deg.E]
            * cfig_name => name of figure to create, extension is important as it will tell what format to use!
            * nsubsamp  => subsampling level, useful for high resolution grids
            * nxcld_n   => number of raws to exclude at the north fold
            * ldark     => hell yeah, make it dark babe!
    """
    #
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from mpl_toolkits.basemap import Basemap
    from mpl_toolkits.basemap import shiftgrid
    #
    ishape   = np.shape(Xglamf)
    (ny,nx,) = ishape
    if len(ishape) != 2:
        print('ERROR: `PlotGridGlobe()` => coordinate arrays should have 2 dimensions!'); sys.exit(0)
    if np.shape(Xgphif) != ishape:
        print('ERROR: `PlotGridGlobe()` => coordinate arrays should have the same shape!'); sys.exit(0)


    # Providing glamt and gphit means we want to plot one of them in the background:
    l_show_lat_bg = (np.shape(Xglamt)==np.shape(Xglamf)) and (np.shape(Xgphit)==np.shape(Xgphif))



        
    vsporg = [0., 0., 1., 1.]
    col_bg = 'w' ; col_fg = 'k' ; col_gr = 'b' ; col_cl = 'k' ; col_fc = '0.85'
    if ldark:
        col_bg = 'k' ; col_fg = 'w' ; col_gr = 'w' ; col_cl = '0.5' ; col_fc = '0.15'

    fig = plt.figure(num = 1, figsize=(nzoom*7.,nzoom*7.), dpi=rdpi, facecolor=col_fg, edgecolor=col_fg)
    ax  = plt.axes(vsporg, facecolor=col_bg)

    if   chemi== 'N':
        if lNPzoom:
            rr = hres/0.027
            carte = Basemap(projection='ortho', lon_0=lon0, lat_0=89.9999, llcrnrx=-rr*20000, llcrnry=-rr*20000, urcrnrx=rr*20000, urcrnry=rr*20000, resolution='c')
        else:
            carte = Basemap(projection='ortho', lon_0=lon0, lat_0= lat0, resolution='c')
    elif chemi== 'S':
        carte = Basemap(projection='ortho', lon_0=lon0, lat_0=-lat0, resolution='c')
    else:
        print('ERROR: `PlotGridGlobe()` => only "N" and "S" hemispheres are known!'); sys.exit(0)


    if l_show_lat_bg:
        rl0 = 90. - 2*hres
        # Add field of latitude:
        x0,y0 = carte(Xglamf[:,:], Xgphif[:,:])
        XP = np.ma.masked_where( Xgphit[:,:]<rl0, Xgphit[:,:] )
        carte.pcolor( x0, y0, XP, cmap=plt.get_cmap('cubehelix_r'), norm=colors.Normalize(vmin=rl0,vmax=90.,clip=False), zorder=1 )
        
    # Vertical lines connecting F-points:
    for ji in range(0,nx,nsubsamp):
        x0,y0 = carte(Xglamf[::nsubsamp,ji], Xgphif[::nsubsamp,ji])
        ftv = carte.plot( x0, y0, color=col_gr, linestyle='-', linewidth=linew, marker=None )

    # Horizontal lines connecting F-points:
    for jj in range(0,ny-nxcld_n,nsubsamp):
        x0,y0 = carte(Xglamf[jj,::nsubsamp], Xgphif[jj,::nsubsamp])
        fth = carte.plot( x0, y0, color=col_gr, linestyle='-', linewidth=linew, marker=None )
        
    carte.drawcoastlines(linewidth=1.,  color=col_cl)
    carte.fillcontinents( color=col_fc )


    plt.savefig(cfig_name, dpi=rdpi, orientation='portrait', transparent=(not l_show_lat_bg))
    print(cfig_name+' created!\n')
    plt.close(1)
    return 0



def CheckNBfolding( pX, jperio, cpnt='T', cwhat='longitude' ):
    '''
       Assuming we are dealing with an ORCA grid, look at the periodic northern
       boundary condition (bi-polar: Siberia and Canadian Archipelago).

    * pX:      2D array of latitudes or longitudes
    * jperio: just as in the NEMO world => what type of periodicity/folding to expect
              - jperio == 4 => North fold T-point pivot (ex: ORCA2, ORCA025, ORCA12, etc.)
              - jperio == 6 => North fold F-point pivot (ex: ORCA1, ORCA05)
    * cwhat:  string 'longitude' or 'latitude'
    
    '''

    print('\n *** `CheckNBfolding`: jperio='+str(jperio)+' / Testing at points "'+str(cpnt)+'"')
    
    if not cpnt in ['T','U','V','F']:
        print('util_orca.CheckNBfolding: ERROR => grid point type unknown: '+cpnt)
        sys.exit(0)
    
    (ny,nx) = np.shape(pX)

    
    if   jperio==4:

        if cpnt=='T':
            # For T-point fields:
            zv1 = pX[ ny-1 ,    1:nx//2-1      ]   ; # Fortran: pX(2:nx/2,ny)
            zv2 = pX[ ny-3 , nx-1:nx-nx//2+1:-1]   ; # Fortran: pX(nx:nx-nx/2+2:-1,ny-2)
            #
        if cpnt=='U':
            # For U-point fields:
            zv1 = pX[ ny-1 ,    1:nx//2-1      ]   ; # Fortran: Xtest(2:nx/2,ny) 
            zv2 = pX[ ny-3 , nx-2:nx-nx//2  :-1]   ; # Fortran: Xtest(nx-1:nx-nx/2+1:-1,ny-2)
            #
        if cpnt=='V':
            # For V-point fields:
            zv1 = pX[ ny-1 ,    1:nx//2-1      ]   ; # Fortran: Xtest(2:nx/2,ny)
            zv2 = pX[ ny-4 , nx-1:nx-nx//2+1:-1]   ; # Fortran: Xtest(nx:nx-nx/2+2:-1,ny-3)
            
        else:
            print('DO ME: point type =',cpnt)
        
    elif jperio==6:

        if cpnt=='T':
            # For T-point fields:
            zv1 = pX[ny-1 ,    1:nx//2-1    ]      ; # Fortran: pX(2:nx/2,ny)
            zv2 = pX[ny-2 , nx-2:nx-nx//2:-1]      ; # Fortran: pX(nx-1:nx-nx/2+1:-1,ny-1)
        else:
            print('DO ME: point type =',cpnt)

    # ....
    if len(zv1) != len(zv2):
        print('ERROR: len(v1) != len(v2) ! ', len(zv1), len(zv2))
        exit(0)


    diff = np.sum(np.abs(zv1-zv2))

    print('\n *** DIFF sum(|v1-v2|) = ', diff)

    if diff==0.:
        print('\n ==> SUCCESSFULY PASSED !  :D\n')
        irtrn = 0
    else:
        print('\n ==> NOT PASSED !  :(\n')
        irtrn = -1
        print(zv1[::10])
        print(zv2[::10])
    
    

    return irtrn

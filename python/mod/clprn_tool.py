
import sys
import os
import numpy as nmp


#    #Lfinite = nmp.isfinite(XF)
#    #Lshit   = nmp.logical_not(Lfinite)
#    #idx_bad = nmp.where(Lshit)


ris2 = 1./nmp.sqrt(2.)


def chck4f(cf, script_name=''):

    cmesg = 'ERROR: File '+cf+' does not exist !!!'
    if script_name != '': cmesg = 'ERROR in script '+script_name+': File '+cf+' does not exist !!!'

    if not os.path.exists(cf):
        print(cmesg) ; sys.exit(0)
    else:
        print(' *** will open file '+cf)



def check_env_var(cnm, list):

    # * cnm  : string for the name of the script calling this function
    # * list : a list of string supposed to be environement variables
    #
    # Returns a dictionary containing the variable names and their respective content

    print("\n In "+cnm+" :")

    env_var_dic = {}  ; # empty dictionary

    for cv in list:
        cenv = os.getenv(cv)
        if cenv is None:
            print(" ERROR in "+cnm+":")
            print("  => the {} environement is not set".format(cv))
            sys.exit(0)
        env_var_dic[cv] = cenv
        print(" *** "+cv+" => "+cenv)

    print("")
    return env_var_dic





def round_to_multiple_of(x, prec=2, base=0.5):
    ## Round to non-integer values, such as 0.5 !
    return round(base * round(float(x)/base),prec)

def int_as_multiple_of(x, base=5):
    # Closest integer multiple of base:
    return int(base * round(float(x)/base))


def lon_180_180(x):
    rsign = 1.
    if 180.-x < 0.: rsign=-1.
    return rsign*min(x,abs(x-360.))

def long_to_m180_p180(xx):
    ## Forces longitude to be in the -180:180 frame...
    ## xx: longitude
    xx   = xx % 360.
    rlon = copysign(1.,180.-xx)*min(xx, abs(xx-360.)) ;
    return rlon
#
def long_to_m180_p180_vct(xx):
    ## Forces longitude to be in the -180:180 frame...
    ## xx: longitude
    xx   = xx % 360.
    xlon = nmp.copysign(1.,180.-xx)*nmp.minimum(xx, nmp.abs(xx-360.)) ;
    return xlon
                        

def get_sections_from_file(cfile):
    list_sections = []
    f = open(cfile, 'r') ; cread_lines = f.readlines() ; f.close()
    jl=0 ; l_stop=False
    while not l_stop:        
        ll = cread_lines[jl]
        ls = ll.split()
        cc = ls[0]
        if cc == 'EOF':
            l_stop=True            
        elif cc[0] != '#':
            if cc != 'EOF' and cc[0:13] != 'ref_temp_sali':
                print(ls)
                list_sections.append(cc)
                jl=jl+1 ; # to skip coordinates line
        else:
            print('  ....  ')
        jl=jl+1    
    return list_sections


def iaxe_tick(ny):
    # I want 20 ticks on the absciss axe and multiple of 5
    itick = int( max( 1 , min(ny/20 , max(ny/20,5)/5*5) ) )
    if itick == 4 or itick == 3: itick = 5
    if ny >=  16 and itick == 1: itick = 2
    if ny >=  45 and ny < 100: itick = 5
    if ny >= 100 and ny < 250: itick = 10
    if ny >= 250 and ny < 750: itick = 25
    if ny >= 750 and ny <2000: itick = 50
    if ny >= 2000:             itick = 100
    return itick


def monthly_2_annual(vtm, XDm):

    # Transform a montly time series into an annual time series

    nbm = len(vtm)
    if len(nmp.shape(XDm)) == 1:
        # XDm is a vector
        nbcol = 1
        nt    = len(XDm)
    else:
        # XDm is an array
        [ nbcol, nt ] = nmp.shape(XDm)

    #print(nt


    if nt < nbm:
        print('ERROR: vmonthly_2_vannual.clprn_tool.py => vt and data disagree in size!')
        print('      => size vt = '+str(nbm)+' / size data = '+str(nt))
        sys.exit(0)
    if nt > nbm:
        print('WARNING: vmonthly_2_vannual.clprn_tool.py => vt and data disagree in size!')
        print('      => size vt = '+str(nbm)+' / size data = '+str(nt))
        print('      => deleting the last '+str(nt-nbm)+' of data array data...')
        #if nbcol == 1:
        XDm = nmp.delete(XDm, nmp.arange(nbm,nt))
        print('      => new shape of data =', nmp.shape(XDm),'\n')


    if nbm%12 != 0: print('ERROR: vmonthly_2_vannual.clprn_tool.py => not a multiple of 12!'); sys.exit(0)

    nby = nbm/12
    vty = nmp.zeros(nby)
    XDy = nmp.zeros((nbcol,nby))

    #print('DEBUG: monthly_2_annual.clprn_tool.py => nbm, nby, nbcol:', nbm, nby, nbcol)

    for jy in range(nby):
        jt_jan = jy*12
        vty[jy] = nmp.trunc(vtm[jt_jan]) + 0.5 ; #  1992.5, not 1992

        if nbcol == 1:
            XDy[0,jy] = nmp.mean(XDm[jt_jan:jt_jan+12])
        else:
            XDy[:,jy] = nmp.mean(XDm[:,jt_jan:jt_jan+12], axis=1)

    if nbcol == 1:
        return vty, XDy[0,:]
    else:
        #print('DEBUG: monthly_2_annual.clprn_tool.py => shape(vty):', nmp.shape(vty))
        return vty, XDy



def find_ij_region_box(vbox4, VX, VY):

    [x_min, y_min, x_max, y_max ] = vbox4

    print('')
    print('clprn_tool.find_ij_region_box : x_min, y_min, x_max, y_max => ', x_min, y_min, x_max, y_max)


    # fixing longitude:
    # ~~~~~~~~~~~~~~~~~
    if x_min < 0. : x_min = x_min + 360.
    if x_max < 0. : x_max = x_max + 360.

    VXtmp = nmp.zeros(len(VX)) ; VXtmp[:] = VX[:]
    idx = nmp.where(VX[:] < 0.0) ; VXtmp[idx] = VX[idx] + 360.


    # fixing latitude:
    # ~~~~~~~~~~~~~~~~~

    # Is latitude increasing with j ?
    jy_inc = 1
    if VY[1] < VY[0]: jy_inc = -1

    #print(jy_inc

    #VYtmp = nmp.zeros(len(VY)) ; VYtmp[:] = VY[:]

    j_y_min = find_index_from_value( y_min, VY )
    j_y_max = find_index_from_value( y_max, VY )
    i_x_min = find_index_from_value( x_min, VXtmp )
    i_x_max = find_index_from_value( x_max, VXtmp )

    if i_x_min == -1 or i_x_max == -1 or j_y_min == -1 or j_y_max == -1:
        print('ERROR: clprn_tool.find_ij_region_box, indiex not found')
        sys.exit(0)

    if jy_inc == -1: jdum = j_y_min; j_y_min = j_y_max; j_y_max = jdum

    #print('  * i_x_min = ', i_x_min, ' => ', VX[i_x_min])
    #print('  * j_y_min = ', j_y_min, ' => ', VY[j_y_min])
    #print('  * i_x_max = ', i_x_max, ' => ', VX[i_x_max])
    #print('  * j_y_max = ', j_y_max, ' => ', VY[j_y_max])
    #print('\n')

    return [ i_x_min, j_y_min, i_x_max, j_y_max ]

#-----------------------------------


def read_ascii_column(cfile, ivcol2read):
    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # INPUT
    #       cfile       : ASCII file
    #       ivcol2read  : vector containing indices of colum to be read (ex: [0, 1, 4])
    #
    # OUTPUT
    #      Xout         : array containg the extracted data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    chck4f(cfile)
    f = open(cfile, 'r')
    cread_lines = f.readlines()
    f.close()
    #
    nbcol = len(ivcol2read) ; #print("nbcol = ", nbcol
    #
    # Need to know how many "non-comment" lines:
    jl = 0
    for ll in cread_lines:
        ls = ll.split()
        if ls[0] != '#': jl = jl + 1
    nbl = jl
    #print('number of lines = ', nbl ; sys.exit)
    #
    Xout  = nmp.zeros((nbcol,nbl))
    #
    jl = -1
    for ll in cread_lines:
        ls = ll.split()
        if ls[0] != '#':
            jl = jl+1
            jc = -1
            for icol in ivcol2read:
                jc = jc+1
                Xout[jc,jl] = float(ls[icol])
    #
    return Xout





def get_min_max_df(ZZ, ndf):

    import math

    # RETURNS rounded Min and Max of array ZZ as well as the contour interval for "ndf" contours...

    # Testing where the array has finite values (non nan, non infinite)
    Lfinite = nmp.isfinite(ZZ)
    idx_good = nmp.where(Lfinite)

    zmin = nmp.amin(ZZ[idx_good]) ; zmax = nmp.amax(ZZ[idx_good])

    if abs(zmin) >= abs(zmax):
        zmax = -zmin
    else:
        zmin = -zmax


    zmin0 = zmin ; zmax0 = zmax

    #print(' Before in clprn_tool.get_min_max_df: zmin, zmax =', zmin, zmax)

    rmagn = 10.**(int(math.log10(zmax)))

    zmin = round(zmin/rmagn - 0.5)
    zmax = round(zmax/rmagn + 0.5)

    zdf0 = (zmax - zmin)/ndf

    zdf = 0.0 ; idec = 0
    while zdf == 0.0:
        idec = idec + 1
        zdf = round(zdf0,idec)

    if idec >= 1: zmin = round(zmin,idec-1) ; zmax = round(zmax,idec-1)

    zmin = zmin*rmagn
    zmax = zmax*rmagn

    zdf = (zmax - zmin)/ndf


    while zmax - zdf >= zmax0 and zmin + zdf <= zmin0:
        zmax = zmax - zdf
        zmin = zmin + zdf

    if abs(zmin) < zmax: zmax = abs(zmin)
    if abs(zmin) > zmax: zmin = -zmax

    # Might divide zdf by 2 of zmax and zmin really decreased...
    rn1 = (zmax-zmin)/zdf
    nn = int(round(float(ndf)/rn1,0))
    zdf = zdf/nn

    fact = 10**(-(int(math.log10(zdf))-1))
    zdf = round(zdf*fact,0)
    zdf = zdf/fact

    # we want zmax and rmin to be mutilples of zdf :
    zmax = zdf*int(round(zmax/zdf))
    zmin = zdf*int(round(zmin/zdf))
    
    return [ zmin, zmax, zdf ]





def find_index_from_value( val, VX ):
    if val > nmp.max(VX) or val < nmp.min(VX):
        print('ERROR: find_index_from_value.clprn_tool => value "'+str(val)+'"outside range of Vector!')
        print(VX[:]) ; print(' => value =', val)
        sys.exit(0)
    jval = -1; jj = 0 ; lfound = False
    while not lfound:
        if VX[jj] <= val and VX[jj + 1] > val:
            jval = jj+1; lfound = True
        jj = jj+1
    return jval


def drown(X, mask, k_ew=-1, nb_max_inc=5, nb_smooth=5):
    '''
    PURPOSE : fills continental areas of field X (defined by mask==0)
    -------   using nearest surrounding sea points (defined by mask==1)
              field X is absoluletly unchanged where mask==1
        
    Input/Output :
    --------------
         * X    :  treated array (float)   / modified!            [2D array]
         * mask :  land-sea mask (integer) / unchanged            [2D array]
    
     Optional :
     ----------
         * k_ew :  east-west periodicity on the input file/grid
                   k_ew = -1  --> no periodicity
                   k_ew >= 0  --> periodicity with overlap of k_ew points
    
         * nb_max_inc : distance (in grid points) of incursion into the continents for the drowning...

         * nb_smooth : number of times the smoother is applied on masked region (mask=0)
                       => default: nb_smooth = 50
    
    Author : Laurent BRODEAU, 2007, as part of SOSIE
             ported to python November 2013
    '''

    cmesg = 'ERROR, clprn_tool.py => drown :'

    nbdim = len(nmp.shape(X))

    if nbdim > 3 or nbdim <2:
        print(cmesg+' size of data array is wrong!!!'); sys.exit(0)


    nt = 1
    l_record = False
    if nbdim == 3: l_record = True

    if l_record:
        if nmp.shape(X[0,:,:]) != nmp.shape(mask):
            print(cmesg+' size of data and mask do not match!!!'); sys.exit(0)
        (nt,nj,ni) = nmp.shape(X)
    else:
        if nmp.shape(X) != nmp.shape(mask):
            print(cmesg+' size of data and mask do not match!!!'); sys.exit(0)
        (nj,ni) = nmp.shape(X)

    if nmp.sum(mask) == 0 :
        print('The mask does not have sea points! Skipping drown!')
        return

    Xtemp = nmp.zeros((nj,ni))

    for jt in range(nt):

        if l_record:
            Xtemp[:,:] = X[jt,:,:]
        else:
            Xtemp[:,:] = X[:,:]

        maskv = nmp.zeros((nj,ni), dtype=nmp.int)
        dold = nmp.zeros((nj,ni))
        xtmp = nmp.zeros((nj,ni))
        mask_coast = nmp.zeros((nj,ni))

        jinc = 0

        maskv[:,:] = mask[:,:]

        for jinc in range(1,nb_max_inc+1):

            dold[:,:] = Xtemp[:,:]

            # Building mask of the coast-line (belonging to land points)
            mask_coast[:,:] = 0

            mask_coast[1:-1,1:-1] = (maskv[1:-1,2:] + maskv[2:,1:-1] + maskv[1:-1,:-2] + maskv[:-2,1:-1])*(-(maskv[1:-1,1:-1]-1))

            if k_ew >= 0:
                # Left LBC:
                mask_coast[1:-1,0]    = (maskv[1:-1,1]    + maskv[2:,0]    + maskv[1:-1,ni-1-k_ew] + maskv[:-2,0]   )*(-(maskv[1:-1,0]   -1))
                # Right LBC:
                mask_coast[1:-1,ni-1] = (maskv[1:-1,k_ew] + maskv[2:,ni-1] + maskv[1:-1,ni-2]      + maskv[:-2,ni-1])*(-(maskv[1:-1,ni-1]-1))

            idx_coast = nmp.where(mask_coast[:,:] > 0)
            #mask_coast[:,:] = 0
            #mask_coast[idx_coast] = 1

            # Extrapolating sea values on that coast line:
            (idx_j_land,idx_i_land) = idx_coast

            ic = 0
            for jj in idx_j_land:
                ji = idx_i_land[ic]

                if ji == 0 and k_ew >= 0:
                    Xtemp[jj,0] = 1./(maskv[jj,1]+maskv[jj+1,0]+maskv[jj,ni-1-k_ew]+maskv[jj-1,0]+
                                   ris2*maskv[jj+1,1]+ris2*maskv[jj+1,ni-1-k_ew]+ris2*maskv[jj-1,ni-1-k_ew]+ris2*maskv[jj-1,1])*(
                        maskv[jj,1]*dold[jj,1] + maskv[jj+1,0]*dold[jj+1,0] +
                        maskv[jj,ni-1-k_ew]*dold[jj,ni-1-k_ew] + maskv[jj-1,0]*dold[jj-1,0] +
                        ris2*maskv[jj+1,1]*dold[jj+1,1] + ris2*maskv[jj+1,ni-1-k_ew]*dold[jj+1,ni-1-k_ew] +
                        ris2*maskv[jj-1,ni-1-k_ew]*dold[jj-1,ni-1-k_ew] + ris2*maskv[jj-1,1]*dold[jj-1,1]  )

                elif ji == ni-1 and k_ew >= 0:
                    Xtemp[jj,ni-1] = 1./(maskv[jj,k_ew]+maskv[jj+1,ni-1]+maskv[jj,ni-2]+maskv[jj-1,ni-1]+
                                   ris2*maskv[jj+1,k_ew]+ris2*maskv[jj+1,ni-2]+ris2*maskv[jj-1,ni-2]+ris2*maskv[jj-1,k_ew])*(
                        maskv[jj,k_ew]*dold[jj,k_ew] + maskv[jj+1,ni-1]*dold[jj+1,ni-1] +
                        maskv[jj,ni-2]*dold[jj,ni-2] + maskv[jj-1,ni-1]*dold[jj-1,ni-1] +
                        ris2*maskv[jj+1,k_ew]*dold[jj+1,k_ew] + ris2*maskv[jj+1,ni-2]*dold[jj+1,ni-2] +
                        ris2*maskv[jj-1,ni-2]*dold[jj-1,ni-2] + ris2*maskv[jj-1,k_ew]*dold[jj-1,k_ew]  )

                else:
                    Xtemp[jj,ji] = 1./(maskv[jj,ji+1]+maskv[jj+1,ji]+maskv[jj,ji-1]+maskv[jj-1,ji]+
                                   ris2*maskv[jj+1,ji+1]+ris2*maskv[jj+1,ji-1]+ris2*maskv[jj-1,ji-1]+ris2*maskv[jj-1,ji+1])*(
                        maskv[jj,ji+1]*dold[jj,ji+1] + maskv[jj+1,ji]*dold[jj+1,ji] +
                        maskv[jj,ji-1]*dold[jj,ji-1] + maskv[jj-1,ji]*dold[jj-1,ji] +
                        ris2*maskv[jj+1,ji+1]*dold[jj+1,ji+1] + ris2*maskv[jj+1,ji-1]*dold[jj+1,ji-1] +
                        ris2*maskv[jj-1,ji-1]*dold[jj-1,ji-1] + ris2*maskv[jj-1,ji+1]*dold[jj-1,ji+1]  )

                ic = ic+1

            # Loosing land for next iteration:
            maskv[idx_coast] = 1

        # Smoothing the what's been done on land:
        if nb_smooth >= 1:

            dold[:,:] = Xtemp[:,:]

            for kk in range(nb_smooth):

                xtmp[:,:] = Xtemp[:,:]

                Xtemp[1:-1,1:-1] = 0.35*xtmp[1:-1,1:-1] + 0.65*0.25*( xtmp[1:-1,2:] + xtmp[2:,1:-1] + xtmp[1:-1,:-2] + xtmp[:-2,1:-1] )

                if k_ew != -1:   # we can use east-west periodicity
                    Xtemp[1:-1,0] = 0.35*xtmp[1:-1,0] + 0.65*0.25*( xtmp[1:-1,1] + xtmp[2:,1] + xtmp[1:-1,ni-1-k_ew] + xtmp[:-2,1] )

                    Xtemp[1:-1,ni-1] = 0.35*xtmp[1:-1,ni-1] + 0.65*0.25*( xtmp[1:-1,k_ew] + xtmp[2:,ni-1] + xtmp[1:-1,ni-2] + xtmp[:-2,ni-1] )


            Xtemp[1:-1,:] = mask[1:-1,:]*dold[1:-1,:] - (mask[1:-1,:]-1)*Xtemp[1:-1,:]


        del maskv, dold, mask_coast, xtmp

        if l_record:
            X[jt,:,:] = Xtemp[:,:]
        else:
            X[:,:]    = Xtemp[:,:]

        Xtemp[:,:] = 0.

    # loop on nt over

    del Xtemp

    return 0





def extend_domain(ZZ, ext_east_deg, skp_west_deg=0):
    #
    # IN:
    # ===
    # ZZ      : array to extend in longitude, 2D field or 1D longitude vector
    # ext_east_deg : eastward extension in degrees...
    # skp_west_deg : what to skip at the west in degrees...
    #
    # OUT:
    # ====
    # ZZx     : zonally-extended array
    #
    #
    #
    vdim = ZZ.shape
    ndim = len(vdim)

    if ndim < 1 or ndim > 3: print('extend_conf.py: ERROR we only treat 1D or 2D arrays...'); sys.exit(0)
    if ndim == 3: [ nz , ny , nx ] = vdim
    if ndim == 2:      [ ny , nx ] = vdim
    if ndim == 1:           [ nx ] = vdim

    nb_skp = int(nx/360.*skp_west_deg)
    nb_ext = int(nx/360.*ext_east_deg) - nb_skp
    nx_ext = nx + nb_ext

    if ndim == 3:
        ZZx  = nmp.zeros((nz, ny, nx_ext))
        for jx in range(nx-skp_west_deg):         ZZx[:,:,jx] = ZZ[:,:,jx+skp_west_deg]
        for jx in range(nx-skp_west_deg, nx_ext): ZZx[:,:,jx] = ZZ[:,:,jx-nx+skp_west_deg]

    if ndim == 2:
        ZZx  = nmp.zeros((ny, nx_ext))
        for jx in range(nx-skp_west_deg):         ZZx[:,jx] = ZZ[:,jx+skp_west_deg]
        for jx in range(nx-skp_west_deg, nx_ext): ZZx[:,jx] = ZZ[:,jx-nx+skp_west_deg]

    if ndim == 1:
        ZZx  = nmp.zeros(nx_ext)
        for jx in range(nx-skp_west_deg):        ZZx[jx] = ZZ[jx+skp_west_deg]
        for jx in range(nx-skp_west_deg,nx_ext): ZZx[jx] = ZZ[jx-nx+skp_west_deg] + 360.


    return ZZx



def mk_zonal(XF, XMSK=[0.], r_mask_from_val=-9999.):
    #**************************************************************************************
    # Computes the zonal average of field XF, ignoring points where XMSK==0.
    #
    # INPUT:
    #        * XF:   2D [ny,nx] or 2D+time [ny,nx,Nt] array of input field
    #  "   opt:
    #        * XMSK: 2D [ny,nx] array, 1 on points to consider, 0 on points to exclude
    #        * r_mask_from_val: instead of providing the 2D mask provide the flag value
    #                           where mask should be 0
    # RETURNS:
    #        * VZ:   1D [ny] array of zonally-averaged XF
    #**************************************************************************************
    #
    vshp = nmp.shape(XF)
    ndim = len(vshp)
    if ndim == 3:
        ( Nt, ny, nx ) = vshp
    elif ndim == 2:
        (     ny, nx ) = vshp
        Nt = 1
    else:
        print(' ERROR (mk_zonal of clprn_tool.py): dimension of your field is weird!')
        sys.exit(0)
    #
    if len(nmp.shape(XMSK)) == 2:
        (n2,n1) = XMSK.shape
        if n2 != ny or n1 != nx:
            print('ERROR: mk_zonal.clprn_tool.py => XF and XMSK do not agree in size!')
            sys.exit(0)
    else:
        # Need to build the mask
        xtmp = nmp.zeros((ny,nx))
        if ndim == 3: xtmp = XF[0,:,:]
        if ndim == 2: xtmp = XF[  :,:]
        XMSK = nmp.zeros((ny,nx))
        idx1 = nmp.where(xtmp > r_mask_from_val + 1.E-6)
        XMSK[idx1] = 1.
        idx1 = nmp.where(xtmp < r_mask_from_val - 1.E-6)
        XMSK[idx1] = 1.
        del xtmp
    #
    VZ = nmp.zeros((Nt,ny))
    #
    vweights = nmp.zeros(ny)
    for jy in range(ny): vweights[jy] = nmp.sum(XMSK[jy,:])
    idx0 = nmp.where(vweights == 0.)
    vweights[idx0] = 1.E12
    #
    for jt in range(Nt):        
        for jy in range(ny):
            if ndim == 3: rmean = nmp.sum(XF[jt,jy,:]*XMSK[jy,:])/vweights[jy]
            if ndim == 2: rmean = nmp.sum(XF[   jy,:]*XMSK[jy,:])/vweights[jy]
            VZ[jt,jy] = rmean
            VZ[jt,idx0] = nmp.nan
    if ndim == 3: return VZ
    if ndim == 2: return VZ[0,:]



def read_coor(cf, ctype='int', lTS_bounds=False):

    # cf : file containing sections or coordinates (name lon1 lat1 lon2 lat2)
    # ctype: 'int' or 'float'
    # lTS_bounds: read and return the bounds (min and max) for T and S !
    
    chck4f(cf)

    f = open(cf, 'r')
    cread_lines = f.readlines()
    f.close()

    vboxes = [] ; vi1 = [] ; vj1 = [] ; vi2 = [] ; vj2 = []
    vTs = [] ; vTl = [] ; vSs = [] ; vSl = [] ; 
    leof = False
    jl   = -1

    while not leof:
        jl = jl + 1
        ll = cread_lines[jl]
        ls = ll.split()
        #
        c1 = ls[0] ; c1 = c1[0]
        #
        if c1 != '#' and ls[0] != 'EOF':
            vboxes.append(ls[0])
            if ctype == 'int':
                vi1.append(int(ls[1]))
                vj1.append(int(ls[2]))
                vi2.append(int(ls[3]))
                vj2.append(int(ls[4]))
            elif ctype == 'float':
                vi1.append(float(ls[1]))
                vj1.append(float(ls[2]))
                vi2.append(float(ls[3]))
                vj2.append(float(ls[4]))
            else:
                print('ERROR: read_coor => ctype must be "int" or "float"')
                sys.exit(0)
            #
            if lTS_bounds:
                # Min and max values for temperature and salinity
                if len(ls) == 9:
                    vTs.append(float(ls[5]))
                    vTl.append(float(ls[6]))
                    vSs.append(float(ls[7]))
                    vSl.append(float(ls[8]))
                elif len(ls) == 5:
                    # no temp. sal. bounds given, chosing some defaults for you...
                    vTs.append(float(-2.)) ; # deg.C
                    vTl.append(float(30.)) ; # deg.C
                    vSs.append(float(33.)) ; # PSU
                    vSl.append(float(37.)) ; # PSU
                else:
                    print('ERROR: read_coor => something wrong in line "'+ls[0]+'" !')
                    sys.exit(0)
        #
        if ls[0] == 'EOF': leof = True

    if lTS_bounds:
        return vboxes, vi1, vj1, vi2, vj2,  vTs, vTl, vSs, vSl
    else:
        return vboxes, vi1, vj1, vi2, vj2


def test_nb_years(vt, cd):
    nb_rec = len(vt)
    # Monthly or Annual time-series?
    idt = int( vt[1] - vt[0] )
    if (idt == 0) and (nb_rec%12 == 0):
        # Montly time-series
        nb_m = nb_rec
        nb_y = nb_rec/12
    elif idt == 1:
        # Annual time-series
        nb_m = -1
        nb_y = nb_rec
    else:
        print('ERROR: '+csn+' for diag='+cd)
        print('       => the time vector seems to be neither monthly nor annual!')
        print('       => Nb. rec. = '+str(nb_rec))
        sys.exit(0)
    ttick = iaxe_tick(nb_y)
    return (nb_y, nb_m, nb_rec, ttick)



def var_and_signs( csin ):
    # Ex: if csin = 'Qsol+Qr-Qlat-Qsens' (or '+Qsol+Qr-Qlat-Qsens')
    #     this function will return:
    #     ['Qsol', 'Qr', 'Qlat', 'Qsens'] , [1.0, 1.0, -1.0, -1.0]
    from re import split as splt
    #
    sgn1 = '+'
    cvar = splt('\-|\+',csin)
    if cvar[0] == '':
        if csin[0] == '-': sgn1 = '-'
        cvar = cvar[1:]
    ccum = ''
    for cc in cvar: ccum=ccum+'|'+cc
    ccum = ccum+'|'
    ctt=splt("'"+ccum+"'",csin)
    csgn = ctt[:-1]
    isgn = []
    for cs in csgn: isgn.append(float(cs+'1'))
    return cvar, isgn
    



def smoother(X, msk, nb_smooth=5):

    ### Do boundaries!!!
    
    cmesg = 'ERROR, clprn_tool.py => smoother :'

    nbdim = len(nmp.shape(X))

    if nbdim != 2:
        print(cmesg+' size of data array is wrong!!!'); sys.exit(0)
    
    (nj,ni) = nmp.shape(X)

    xtmp = nmp.zeros((nj,ni))

    for ii in range(nb_smooth):

        xtmp[:,:] = X[:,:]*msk[:,:]

        X[1:-1,1:-1] = 0.35*xtmp[1:-1,1:-1] + ( 0.65*( xtmp[1:-1,2:] + xtmp[2:,1:-1] + xtmp[1:-1,:-2] + xtmp[:-2,1:-1] \
                                    + ris2*( xtmp[2:,2:]   + xtmp[:-2,2:]  + xtmp[:-2,:-2] + xtmp[2:,:-2]  )  ) ) \
                                    / nmp.maximum( msk[1:-1,2:] + msk[2:,1:-1] + msk[1:-1,:-2] + msk[:-2,1:-1] \
                                  + ris2*( msk[2:,2:]   + msk[:-2,2:]  + msk[:-2,:-2] + msk[2:,:-2]  ) \
                                           , 1.E-6 )
                           
        X[:,:] = X[:,:]*msk[:,:]

    del xtmp

    return



def symetric_range( pmin, pmax ):
    # Returns a symetric f-range that makes sense for the anomaly of "f" we're looking at...
    from math import floor, copysign, log, ceil
    zmax = max( abs(pmax) , abs(pmin) )
    romagn = floor(log(zmax, 10)) ; # order of magnitude of the anomaly  we're dealing with
    rmlt = 10.**(int(romagn)) / 2.
    frng = copysign( ceil(abs(zmax)/rmlt)*rmlt , zmax)
    return frng



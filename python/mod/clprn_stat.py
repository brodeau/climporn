
# BaraKuda statistics misc. functions...

import sys
import numpy as nmp
import random

def least_sqr_line(ZX,ZY):
    Nx = len(ZX) ; Ny = len(ZY)
    if Nx != Ny: print 'clprn_stat.least_sqr_line: ERROR in size!'; sys.exit(0)
    Sx  = nmp.sum(ZX[:])       ; Sy  = nmp.sum(ZY[:])
    Sxx = nmp.sum(ZX[:]*ZX[:]) ; Sxy = nmp.sum(ZX[:]*ZY[:])
    rA = ( Sx*Sy - Nx*Sxy ) / ( Sx*Sx - Nx*Sxx )
    rB = ( Sy - rA*Sx ) / Nx
    return ( rA, rB )


def stat_serie(ZX,ZY):

    nt = len(ZX)

    # Least-square curve:
    [ A, B ] = least_sqr_line(ZX,ZY)

    ZYn = nmp.zeros(nt)

    # Removing the trend:
    ZYn = ZY - A*ZX[:] + B

    # Mean
    moy_y = nmp.mean(ZYn[:])

    # Ecart-type des Y:
    sigma = nmp.sqrt(nmp.sum((ZYn[:] - moy_y)*(ZYn[:] - moy_y))/nt)

    print 'A, B, moy_y, sigma_y =', A, B, moy_y, sigma

    ZYn = (ZYn - moy_y)/sigma

    return ZYn





def std_dev(VX):
    Np = len(VX)
    mn = nmp.mean(VX)
    sgm = nmp.sqrt( nmp.sum((VX[:] - mn)*(VX[:] - mn)) / Np )
    return sgm



def correl(VX,VY, Np_force=0):

    # Compute the Pearson correlation coefficient
    Nx = len(VX)      ; Ny = len(VY)
    if Nx != Ny: print 'clprn_stat.correl: ERROR in size! => ', Nx, Ny; sys.exit(0)

    Np = Nx
    if Np_force > 1: Np = Np_force

    mx = nmp.mean(VX[:])  ; my = nmp.mean(VY[:])
    sx = std_dev(VX)  ; sy = std_dev(VY)
    Rxy = nmp.sum((VX[:] - mx)*(VY[:] - my)) / ( Np*sx*sy )

    return Rxy


def cross_correl(VX, VY, klag):

    # Cross-correlation for a lag of klag points...

    Np = len(VX)
    if Np != len(VY): print 'clprn_stat.cross_correl: ERROR in size! => ', Np, len(VY); sys.exit(0)

    if klag >= 0:
        Rxyl = correl(VX[:Np-klag], VY[klag:])
    else:
        Rxyl = correl(VX[-klag:]  , VY[:Np+klag])

    return Rxyl



def auto_correl(VX, klag):

    # Auto-correlation for a lag of klag points...

    Rxx = cross_correl(VX, VX, klag)

    return Rxx






def running_mean_3(VV, l_fill_bounds=False):

    nv  = len(VV) ;     Vrm = nmp.zeros(nv)

    for jv in nmp.arange(1,nv-1):
        Vrm[jv] = (VV[jv-1] + VV[jv] + VV[jv+1]) / 3.

    Vrm[0] = nmp.nan ; Vrm[nv-1:nv] = nmp.nan

    if l_fill_bounds:
        jv = 0
        Vrm[jv] = ( VV[jv] + VV[jv+1] + VV[jv+2] ) / 3.   # Asymetric !!!
        jv = nv-1
        Vrm[jv] = ( VV[jv] + VV[jv-1] + VV[jv-2] ) / 3.   # Asymetric !!!

    return Vrm






def running_mean_5(VV, l_fill_bounds=False):

    nv  = len(VV) ;     Vrm = nmp.zeros(nv)

    for jv in nmp.arange(2,nv-2):
        Vrm[jv] = (VV[jv-2] + VV[jv-1] + VV[jv] + VV[jv+1] + VV[jv+2]) / 5.

    Vrm[0:2] = nmp.nan ; Vrm[nv-2:nv] = nmp.nan

    if l_fill_bounds:
        jv = 1
        Vrm[jv] = ( VV[jv-1] + VV[jv] + VV[jv+1] ) / 3.   # Asymetric !!!
        jv = nv-jv-1
        Vrm[jv] = ( VV[jv-1] + VV[jv] + VV[jv+1] ) / 3.   # Asymetric !!!

        jv = 0
        Vrm[jv] = ( VV[jv] + VV[jv+1] + VV[jv+2] ) / 3.   # Asymetric !!!
        jv = nv-1
        Vrm[jv] = ( VV[jv] + VV[jv-1] + VV[jv-2] ) / 3.   # Asymetric !!!

    return Vrm




def running_mean_11(VV, l_fill_bounds=False):

    # Input:
    #          VV: inut vector
    #          l_fill_bounds: fill or leave as 'nan' the 4 first and 4 last value of the vector
    #
    # Output:
    #          Vrm : shame size as VV

    nv  = len(VV)
    Vrm = nmp.zeros(nv)
    #
    for jv in nmp.arange(5,nv-5):
        Vrm[jv] = ( VV[jv-5] + VV[jv-4] + VV[jv-3] + VV[jv-2] + VV[jv-1] + VV[jv]
                    + VV[jv+1] + VV[jv+2] + VV[jv+3] + VV[jv+4] + VV[jv+5] ) / 11.

    # About begining and end:
    Vrm[0:5] = nmp.nan ; Vrm[nv-5:nv] = nmp.nan

    if l_fill_bounds:

        jv = 4
        Vrm[jv] = ( VV[jv-4]+VV[jv-3]+VV[jv-2]+VV[jv-1]+VV[jv]+VV[jv+1]+VV[jv+2]+VV[jv+3]+VV[jv+4] ) / 9.
        jv = nv-jv-1
        Vrm[jv] = ( VV[jv-4]+VV[jv-3]+VV[jv-2]+VV[jv-1]+VV[jv]+VV[jv+1]+VV[jv+2]+VV[jv+3]+VV[jv+4] ) / 9.

        jv = 3
        Vrm[jv] = ( VV[jv-3] + VV[jv-2] + VV[jv-1] + VV[jv] + VV[jv+1] + VV[jv+2] + VV[jv+3] ) / 7.
        jv = nv-jv-1
        Vrm[jv] = ( VV[jv-3] + VV[jv-2] + VV[jv-1] + VV[jv] + VV[jv+1] + VV[jv+2] + VV[jv+3] ) / 7.

        jv = 2
        Vrm[jv] = ( VV[jv-2] + VV[jv-1] + VV[jv] + VV[jv+1] + VV[jv+2] ) / 5.
        jv = nv-jv-1
        Vrm[jv] = ( VV[jv-2] + VV[jv-1] + VV[jv] + VV[jv+1] + VV[jv+2] ) / 5.

        jv = 1
        Vrm[jv] = ( VV[jv-1] + VV[jv] + VV[jv+1] + VV[jv+2]+VV[jv+3] ) / 5.   # Asymetric !!!
        jv = nv-jv-1
        Vrm[jv] = ( VV[jv-1] + VV[jv-2]+VV[jv-3] + VV[jv] + VV[jv+1] ) / 5.   # Asymetric !!!

        jv = 0
        Vrm[jv] = ( VV[jv] + VV[jv+1] + VV[jv+2] + VV[jv+3] + VV[jv+4]) / 5.   # Asymetric !!!
        jv = nv-1
        Vrm[jv] = ( VV[jv] + VV[jv-1] + VV[jv-2] + VV[jv-3] + VV[jv-4]) / 5.   # Asymetric !!!


    return Vrm[:]





def running_mean_21(VV):

    # Input:
    #          VV: inut vector
    #
    # Output:
    #          Vrm : shame size as VV

    nv  = len(VV)
    Vrm = nmp.zeros(nv)

    for jv in nmp.arange(10,nv-10):
        Vrm[jv] = (   VV[jv-10]+ VV[jv-9] + VV[jv-8] + VV[jv-7] + VV[jv-6] \
                    + VV[jv-5] + VV[jv-4] + VV[jv-3] + VV[jv-2] + VV[jv-1] \
                    + VV[jv] \
                    + VV[jv+1] + VV[jv+2] + VV[jv+3] + VV[jv+4] + VV[jv+5] \
                    + VV[jv+6] + VV[jv+7] + VV[jv+8] + VV[jv+9] + VV[jv+10] ) / 21.

    # About begining and end:
    Vrm[0:10] = nmp.nan ; Vrm[nv-10:nv] = nmp.nan

    return Vrm[:]



def running_mean_31(VV):

    # Input:
    #          VV: inut vector
    #
    # Output:
    #          Vrm : shame size as VV

    nv  = len(VV)
    Vrm = nmp.zeros(nv)

    for jv in nmp.arange(15,nv-15):
        Vrm[jv] = (   VV[jv-15]+ VV[jv-14] +VV[jv-13] +VV[jv-12] +VV[jv-11] \
                    + VV[jv-10]+ VV[jv-9] + VV[jv-8] + VV[jv-7] + VV[jv-6] \
                    + VV[jv-5] + VV[jv-4] + VV[jv-3] + VV[jv-2] + VV[jv-1] \
                    + VV[jv] \
                    + VV[jv+1] + VV[jv+2] + VV[jv+3] + VV[jv+4] + VV[jv+5] \
                    + VV[jv+6] + VV[jv+7] + VV[jv+8] + VV[jv+9] + VV[jv+10] \
                    + VV[jv+11]+ VV[jv+12]+ VV[jv+13]+ VV[jv+14]+ VV[jv+15] ) / 31.

    return Vrm[:]



def shuffle(x):
    lx = list(x)
    random.shuffle(lx)
    return nmp.asarray(lx)


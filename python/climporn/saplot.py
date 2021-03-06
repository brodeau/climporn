#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
#   ///// https://github.com/brodeau/gonzag \\\\\
#
#       L. Brodeau, 2021
#
############################################################################

import numpy as nmp
#
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from .utils import sym_round_bounds

# Look and feel for the plot:
clr_sat = '#AD0000'
clr_mod = '#008ab8'


def PlotSLASeries(vt, VS, VM, cfig, clabS='Satellite', clabM='Model' ):
    ras = nmp.mean(VS) ; ram = nmp.mean(VM)
    ymin, ymax, dy = sym_round_bounds(min(nmp.min(VM-ram),nmp.min(VS-ras)), max(nmp.max(VM-ram), nmp.max(VS-ras)), base=0.1 )
    fig = plt.figure(num = 1, figsize=(12,7), facecolor='w', edgecolor='k')
    ax  = plt.axes([0.07, 0.22, 0.9, 0.75])
    #ax.set_xticks(vt[::xticks_d])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    plt.xticks(rotation='60')
    plt.plot(vt, vt*0.0         , '-', color='k',                   label=None,  zorder=5)
    plt.plot(vt, VS-ras, '.', color=clr_sat, markersize=6, alpha=0.5, label=clabS, zorder=10)
    plt.plot(vt, VM-ram, '.', color=clr_mod, markersize=6, alpha=0.5, label=clabM, zorder=15)
    plt.yticks( nmp.arange(ymin, ymax+dy, dy) )
    ax.set_ylim(ymin,ymax)
    #ax.set_xlim(vt[0],vt[-1])
    plt.ylabel('SSH [m]')
    ax.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(bbox_to_anchor=(0.55, 1.), ncol=1, shadow=True, fancybox=True)
    plt.savefig(cfig, dpi=120, transparent=False)
    plt.close(1)
    return 0



def PlotSegmentTrack( vt, VS, VM, cfig, ctitle='', clabS='Satellite', clabM='Model' ):
    '''
    '''
    ii=len(vt)//300 ; xticks_d=5*max(ii-ii%10,1)
    #
    ymin, ymax, dy = sym_round_bounds( min(nmp.min(VM[:]),nmp.min(VS[:])), max(nmp.max(VM[:]), nmp.max(VS[:])), base=0.1 )
    #
    fig = plt.figure(num = 1, figsize=(12,7.4), facecolor='w', edgecolor='k')
    ax = plt.axes([0.07, 0.22, 0.88, 0.73])
    ax.set_xticks(vt[::xticks_d])
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    plt.xticks(rotation='60')
    plt.plot(vt, vt*0., '-', color='k', linewidth=2, label=None)
    plt.plot(vt, VS[:], '-o', color=clr_sat, linewidth=2, label=clabS, zorder=10)
    plt.plot(vt, VM[:], '-o', color=clr_mod, linewidth=2, label=clabM, zorder=15)
    plt.yticks( nmp.arange(ymin, ymax+dy, dy) )
    ax.set_ylim(ymin,ymax)
    plt.ylabel('SSH [m]')
    ax.set_xlim(vt[0],vt[-1])
    ax.grid(color='k', linestyle='-', linewidth=0.3)
    plt.legend(loc="best", ncol=1, shadow=True, fancybox=True)
    plt.title(ctitle)
    plt.savefig(cfig, dpi=120, transparent=False)
    plt.close(1)
    return 0


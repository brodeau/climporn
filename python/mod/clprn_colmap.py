# Misc. colormaps...

# Check http://matplotlib.org/examples/color/colormaps_reference.html !!!

# Colors:  http://www.pitt.edu/~nisg/cis/web/cgi/rgb.html

import sys
import numpy as nmp

# List of Climporn home-made colormaps:
list_climporn = [ 'blk', 'land', 'land_dark', 'terre', 'cb1', 'eke', 'bathy', 'mld', 'tap1', 'tap2', 'jetblanc', 'amoc',
                  'sst1', 'sst2', 'sst3', 'ice', 'ice_on', 'ice2_on', 'ice3_on', 'ice4_on', 'blanc', 'rms',
                  'sigtr', 'bbr', 'bbr2', 'bbr0', 'bbr_cold', 'bbr_warm',
                  'cold0', 'warm0', 'graylb', 'graylb2', 'sigma', 'sigma0', 'mask', 'on0', 'on1', 'on2', 'on3' ]

# There are also NCVIEW colormaps, and the default Matplotlib colormaps...

l_debug = False

def chose_colmap( cname, log_ctrl=0, exp_ctrl=0 ):

    # 1st is it a ncview colormap ?
    if cname[:7] == 'ncview_':
        lr = ( cname[-2:] == '_r' )
        if lr: cname = cname[:-2]
        M = ncview_ncmap_to_array( cname, l_rev=lr )
        ColorMap = __build_colormap__(M, log_ctrl=log_ctrl, exp_ctrl=exp_ctrl)

        # Maybe a climporn colormap ?
    elif cname in list_climporn or ( cname[-2:] == '_r' and cname[:-2] in list_climporn):
        if l_debug: print('\n *** Getting Climporn colormap "'+cname+'" !')
        x = brkd_cmap(cname)
        ColorMap = x.clrmp(log_ctrl=log_ctrl, exp_ctrl=exp_ctrl)
    else:
        # Then it must be a Matplotlib colormap:
        if log_ctrl or exp_ctrl > 0: print('WARNING: cannot use LOG or EXP colormap with Matplotlib colormaps...')
        from matplotlib.pylab import cm
        import matplotlib.pyplot as mp
        list = mp.colormaps()
        if cname in list:
            # Yes it is!
            if l_debug: print('\n *** Getting Matplotlib colormap "'+cname+'" !')
            fToCall = getattr(cm, cname)
            ColorMap = fToCall
        else:
            print('ERROR: (chose_colmap of clprn_colmap.py) do not know where to get colormap "'+cname+'" !')
            sys.exit(0)

    return ColorMap



def ncview_h_to_array( cname ):
    #
    #########################################################################################
    #
    # Get the NCVIEW colormap in the C header file and return an array ready for a colormap
    #         Author: L. Brodeau, 2017
    #
    # cname: "ncview_" +  name of the NCVIEW colormap as in the NCVIEW header colormap files  [string]
    #   example : 'ncview_rainbow'
    #
    #
    # The environment variable 'DIR_NCVIEW_CMAP' must be set!
    #    => path to the directory containing the NCVIEW header colormap files
    #    => ex: colormap 'rainbow' is defined in header file 'colormaps_rainbow.h''
    #
    ##########################################################################################

    import os
    import re

    dir_ncview_cmap = os.getenv('DIR_NCVIEW_CMAP')
    if dir_ncview_cmap is None:
        print(" ERROR => the {} environement variable is not set".format('DIR_NCVIEW_CMAP'))
        sys.exit(0)

    if cname[:7] != 'ncview_' : print(' ERROR: a ncview colormap should begin with "ncview_" !'); sys.exit(0)
    ncview_name = cname[7:]

    cf_ncview_cmap = dir_ncview_cmap+'/colormaps_'+ncview_name+'.h'
    if not os.path.exists(cf_ncview_cmap):
        print('ERROR: NCVIEW colormap '+cf_ncview_cmap+' not found!'); sys.exit(0)
    if l_debug: print('\n *** Getting NCVIEW colormap "'+ncview_name+'" from file "'+cf_ncview_cmap+'"')

    f = open(cf_ncview_cmap, 'r')
    cread_lines = f.readlines()
    f.close()

    lstarted = False
    vec = []
    for ll in cread_lines:

        ll = re.sub(r'\s', '', ll)
        ll = re.sub(r'{', ',', ll)
        ll = re.sub(r'}', ',', ll)
        ls = re.split(',',ll)

        if lstarted:
            for ve in ls[:-1]: vec.append(float(ve)) ; # [:-1] is to ommit the ',' at the end

        if ls[0] == 'staticintcmap_'+ncview_name+'[]=':
            lstarted = True
            for ve in ls[1:-1]: vec.append(float(ve)) ; # [:-1] is to ommit the ',' at the end

    ctmp = []
    ii = 0
    while ii < len(vec):
        ctmp.append([vec[ii], vec[ii+1], vec[ii+2]])
        ii += 3

    MM = nmp.array(ctmp)/255.

    return MM


def ncview_ncmap_to_array( cname, l_rev=False ):
    #
    # Reads the *.ncmap files !
    #
    from os import path
    import re

    dir_scrpt = path.dirname(path.realpath(__file__))
    dir_ncview_cmap = str.replace( dir_scrpt ,  'python/mod', 'misc/ncview_colormaps')

    if cname[:7] != 'ncview_' : print(' ERROR: a ncview colormap should begin with "ncview_" !'); sys.exit(0)
    ncview_name = cname[7:]
    #
    cf_ncview_cmap = dir_ncview_cmap+'/'+ncview_name+'.ncmap'
    if not path.exists(cf_ncview_cmap):
        print('ERROR: NCVIEW colormap '+cf_ncview_cmap+' not found!'); sys.exit(0)
    if l_debug: print('\n *** Getting NCVIEW colormap "'+ncview_name+'" from file "'+cf_ncview_cmap+'"')
    #
    f = open(cf_ncview_cmap, 'r')
    cread_lines = f.readlines()
    f.close()
    #
    if l_rev:
        vrlines = cread_lines[::-1]
    else:
        vrlines = cread_lines[:]
    #
    vec = []
    for ll in vrlines:
        ls = re.split(' ',ll)
        for ve in ls[:]: vec.append(float(ve))
    #
    ctmp = []
    ii = 0
    while ii < len(vec):
        ctmp.append([vec[ii], vec[ii+1], vec[ii+2]])
        ii += 3
    MM = nmp.array(ctmp)/255.
    #
    return MM








# ===== local ======


def __build_colormap__(MC, log_ctrl=0, exp_ctrl=0):

    import matplotlib.colors as mplc

    [ nc, n3 ] = nmp.shape(MC)

    # Make x vector :
    x =[]
    for i in range(nc): x.append(255.*float(i)/((nc-1)*255.0))
    x = nmp.array(x)
    if log_ctrl > 0: x = nmp.log(x + log_ctrl)
    if exp_ctrl > 0: x = nmp.exp(x * exp_ctrl)
    rr = x[nc-1] ; x  = x/rr

    y =nmp.zeros(nc)
    for i in range(nc): y[i] = x[nc-1-i]

    x = 1 - y ; rr = x[nc-1] ; x  = x/rr

    vred  = [] ; vblue = [] ; vgreen = []

    for i in range(nc):
        vred.append  ([x[i],MC[i,0],MC[i,0]])
        vgreen.append([x[i],MC[i,1],MC[i,1]])
        vblue.append ([x[i],MC[i,2],MC[i,2]])

    cdict = {'red':vred, 'green':vgreen, 'blue':vblue}

    my_cm = mplc.LinearSegmentedColormap('my_colormap',cdict,256)

    return my_cm



#=======================================================================



class brkd_cmap:

    def __init__(self, name):
        self.name = name

    def clrmp(self, log_ctrl=0, exp_ctrl=0):

        cname = self.name

        lrev = False
        if cname[-2:] == '_r':
            lrev = True
            cname = cname[:-2]

        if cname == 'blk':
            M = nmp.array( [
                [ 0. , 0., 0. ], # black
                [ 0. , 0., 0. ]  # black
            ] )

        elif cname == 'land':
            M = nmp.array( [
                [ 0.75 , 0.75, 0.75 ],
                [ 0.75 , 0.75, 0.75 ]
            ] )

        elif cname == 'land_dark':
            M = nmp.array( [
                [ 0.35 , 0.35, 0.35 ],
                [ 0.35 , 0.35, 0.35 ]
            ] )

        elif cname == 'terre':
            M = nmp.array( [
                [ 156./255.,85./255.,54./255. ],
                [ 156./255.,85./255.,54./255. ],
            ] )

        elif cname == 'on0':
            M = nmp.array( [
                [ 1.,1.,1. ],               # blanc
                [ 0.,138./255.,184./255. ], # bleu
                [ 0.,0.,0. ],               # noir
            ] )

        elif cname == 'on1':
            M = nmp.array( [
                [ 1.,1.,1. ], # blanc
                [ 0.,138./255.,184./255. ], # bleu
            ] )

        elif cname == 'on2':
            M = nmp.array( [
                [ 0.,0.,0. ],               # noir
                [ 0.,138./255.,184./255. ], # bleu
                [ 1.,1.,1. ],               # blanc
            ] )

        elif cname == 'on3':
            M = nmp.array( [
                [ 0.,0.,0. ],               # noir                
                [ 0.,138./255.,184./255. ], # bleu                
                [ 1.,1.,1. ],               # blanc                
                [ 1.,237./255.,0 ],         # jaune
            ] )

        elif cname == 'cb1':
            M = nmp.array( [
                [ 255,255,204  ],
                [ 161,218,180  ],
                [ 65,182,196  ],
                [ 44,127,184  ],
                [ 37,52,148  ]
            ] ) / 255.

        elif cname == 'eke':
            M = nmp.array( [
                [ 0.  , 0.0 , 0.2  ], # black
                [ 0.1 , 0.5 , 1.0  ], # blue
                [ 0.2 , 1.0 , 0.0  ], # green
                [ 1.  , 1.0 , 0.0  ], # yellow
                [ 1.  , 0.0 , 0.0  ], # red
                [0.2  , 0.27, 0.07 ] # brown
            ] )

        elif cname == 'bathy':
            M = nmp.array( [
                [ 0.0 , 0.0 , 0.4 ], # dark blue
                [ 0.1 , 0.5 , 1.0  ], # blue
                [ 0.2 , 1.0 , 0.0  ], # green
                [ 1.  , 1.0 , 0.0  ], # yellow
                [ 1.  , 0.0 , 0.0  ], # red
                [0.2  , 0.27, 0.07 ] # brown
            ] )

        elif cname == 'mld':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 0.1 , 0.5 , 1.0 ], # light blue
                [ 0.13, 0.54, 0.13], # dark green
                [ 0.2 , 1.0 , 0.0 ], # light green
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 60./255.,0.,0. ] # uber dark redish brown
            ] )

        elif cname == 'tap1':
            M = nmp.array( [
                [232./255.,254./255.,255./255.], # light blue 1
                [ 0.1 , 0.5 , 1.0 ], #  blue 1
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 60./255.,0.,0. ]   # uber dark redish brown
            ] )

        elif cname == 'tap2':
            M = nmp.array( [
                [232./255.,254./255.,255./255.], # light blue 1
                #[23./255.,170./255.,1.], #light blue 2
                #[132./255.,207./255.,197./255.], # light blue 3
                ##[166./255.,204./255.,255./255.], # light blue 4
                #[ 0.1 , 0.5 , 1.0 ], #  blue 1
                [32./255.,55./255.,145./255.], # blue 2
                #[255./255.,166./255.,198./255.], # pink
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ],  # red
                [ 112./255.,4./255.,4./255. ] # dark redish brown
            ] )

        elif cname == 'jetblanc':
            M = nmp.array( [
                [ 0.6 , 0.0 , 0.8 ], # violet
                [ 0.0 , 0.0 , 0.4 ], # dark blue
                [ 0.1 , 0.5 , 1.0 ], # light blue
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ]  # dark redish brown
            ] )

        elif cname == 'amoc':
            M = nmp.array( [
                [ 0.4 , 0.0 , 0.6 ], # violet
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 1.0 , 1.0 , 1.0 ], # white
                [0.68 , 0.98, 0.98], # light blue
                [ 0.0 , 0.0 , 0.95], # dark blue
                [ 0.2 , 1.0 , 0.0 ], # green
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ] # dark read
            ] )

        elif cname == 'sst1':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 0.4 , 0.0 , 0.6 ], # violet
                [ 0. , 0.2 , 0.99], # dark blue
                [0.68 , 0.98, 0.98], # light blue
                [ 0.13, 0.54, 0.13], # dark green
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ] # dark read
            ] )

        elif cname == 'sst2':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 0.4 , 0.0 , 0.6 ], # violet
                [ 0.0 , 0.0 , 0.95], # dark blue
                [0.68 , 0.98, 0.98], # light blue
                [46./255., 203./255., 35./255.], # green
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ] # dark read
            ] )

        elif cname == 'sst3':
            M = nmp.array( [
                [ 0.4 , 0.0 , 0.6 ], # violet
                [ 0. , 0.2 , 0.99], # dark blue
                [0.68 , 0.98, 0.98], # light blue
                [ 0.13, 0.54, 0.13], # dark green
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ] # dark read
            ] )

        elif cname == 'ice':
            M = nmp.array( [
                [  0. , 0.  , 0.3 ], # dark blue
                [ 0.6 , 0.6 , 0.8 ], # light grey
                [ 0.95 , 0.95 , 0.95 ],  # white
                [ 1.0 , 1.0 , 1.0 ]  # white
            ] )

        elif cname == 'ice_on':
            M = nmp.array( [
                [ 0.,0.,0. ],        # noir   (has to match with coldest color of "on3" !
                [ 25./255. , 102./255. , 114./255. ],
                [ 1.0 , 1.0 , 1.0 ]  # white
            ] )
        elif cname == 'ice2_on':
            M = nmp.array( [
                [ 14./255.,16./255.,16./255. ],        # Dark gray
                [ 42./255.,42./255.,46./255. ],        # Dark gray
                #[ 40./255. , 58./255. , 58./255. ],
                [ 51./255. , 85./255. , 85./255. ],
                [ 70./255. ,104./255. ,104./255. ],
                [ 174./255. , 196./255. , 212./255. ], # light grey blueish
                #[ 215./255. , 232./255. , 246./255. ], # light grey blueish
                [ 90./255. , 90./255. , 90./255. ], # light grey
                [ 1.0 , 1.0 , 1.0 ]  # white
            ] )
            #                [ 0.6 , 0.6 , 0.8 ], # light grey
        elif cname == 'ice3_on':
            M = nmp.array( [
                [  15./255., 15./255., 15./255. ],        # Dark gray
                [ 70./255.,70./255.,70./255. ],        # Dark gray                
                [ 160./255.,160./255.,160./255. ],     # Light gray
                [  95./255.,156./255.,156./255. ],
                [ 51./255. , 85./255. , 85./255. ],
                [ 1.0 , 1.0 , 1.0 ]                   # white
            ] )
            #
        elif cname == 'ice4_on':
            M = nmp.array( [
                [  0./255., 0./255., 0./255. ],        # Dark gray
                [ 70./255.,70./255.,70./255. ],        # Dark gray                
                [ 160./255.,160./255.,160./255. ],        # Light gray
                [ 127./255.,144./255.,144./255. ],  # ?              
                [  95./255.,156./255.,156./255. ],
                [ 72./255. ,118./255.,118./255. ],  # ?
                [ 51./255. , 85./255. , 85./255. ],
                [ 40./255. , 75./255. , 75./255. ], # ? 
                [ 30./255. , 60./255. , 60./255. ], # ?
                [ 1.0 , 1.0 , 1.0 ]  # white
            ] )
            #
        elif cname == 'blanc':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ],  # white
                [ 1.0 , 1.0 , 1.0 ]  # white
            ] )

        elif cname == 'rms':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ],
                [ 0.1 , 0.5 , 1.0 ],
                [ 0.2 , 1.0 , 0.0 ],
                [ 1.0 , 1.0 , 0.0 ],
                [ 1.0 , 0.0 , 0.0 ],
                [ 0.2 , 0.3 , 0.1 ]
            ] )

        elif cname == 'sigtr':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 0.0 , 0.8 , 1.0 ], #light blue
                [ 0.1 , 0.5 , 1.0 ], #light blue
                [ 0.0 , 0.0 , 0.4 ], # blue
                [ 0.0 , 0.4 , 0.0 ], # dark green
                [ 0.1 , 1.0 , 0.0 ], # green
                [ 0.4 , 1.0 , 0.0 ], # vert pomme
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 1.0 , 0.4 , 0.0 ], # orange
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 0.6 , 0.0 , 0.0 ], # red
                [ 0.2 , 0.3 , 0.1 ]  # dark red
            ] )

        elif cname == 'bbr':
            M = nmp.array( [
                [ 0.  , 0. , 0.2 ],
                [ 0.  , 0. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 0. , 0.  ],
                [ 0.6 , 0. , 0.  ]
            ] )

        elif cname == 'bbr2':
            M = nmp.array( [
                [ 0.  , 1. , 1.  ], # cyan
                [ 0.  , 0. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 0. , 0.  ],
                [ 1.  , 1. , 0.  ]  # jaune
            ] )

        elif cname == 'bbr0':
            M = nmp.array( [
                [ 0.  , 1. , 1.  ], # cyan
                [ 0.  , 0. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 0. , 0.  ],
                [ 1.  , 1. , 0.  ]  # jaune
            ] )

        elif cname == 'bbr_cold':
            M = nmp.array( [
                [ 0.  , 1. , 1.  ], # cyan
                [ 0.  , 0. , 1.  ],
                [ 19./255.  , 7./255. , 129./255  ], # dark blue
                [ .1  , .1 , .9  ],
                [ 1.  , 1. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 0.7  , 0. , 0.  ] # Dark red
            ] )

        elif cname == 'bbr_warm':
            M = nmp.array( [
                [ 19./255.  , 7./255. , 129./255  ], # dark blue
                [ 1.  , 1. , 1.  ],
                [ 1.  , 1. , 1.  ],
                [ 0.9  , 0.1 , 0.1  ],
                [ 0.7  , 0. , 0.  ], # Dark red
                [ 1.  , 1. , 0.  ],  # jaune
            ] )

        elif cname == 'cold0':
            M = nmp.array( [
                [ 177./255.  , 250./255. , 122./255. ],   # greenish
                [ 0.  , 1. , 1.  ], # cyan
                [ 7./255.  , 11./255. , 122./255. ], # dark blue
                [ 0.  , 0. , 1.  ], # true blue
                [ 177./255.  , 189./255. , 250./255. ], # light blue
                [ 1.  , 1. , 1.  ],
            ] )

        elif cname == 'warm0':
            M = nmp.array( [
                [ 1.  , 1. , 1.  ],
                [ 255./255.  , 254./255. , 198./255.  ], # very light yellow
                [ 1.  , 1. , 0.  ],  # yellow
                [ 244./255.  , 78./255. , 255./255.  ], # pink
                [ 1.  , 0. , 0.  ], # true red
                [ 139./255.  , 5./255. , 5./255.  ] # dark red
            ] )

        elif cname == 'graylb':
            M = nmp.array( [
                [ 1.  , 1. , 1. ],
                [ 0.1  , 0.1 , 0.1 ]
            ] )

        elif cname == 'graylb2':
            M = nmp.array( [
                [ 0.6  , 0.6 , 0.6 ],
                [ 1.  , 1. , 1. ]
            ] )

        elif cname == 'sigma':
            M = nmp.array( [
                [ 1.0 , 1.0 , 1.0 ], # white
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 0.2 , 1.0 , 0.0 ], # green
                [ 0.1 , 0.5 , 1.0 ], # light blue
                [ 0.0 , 0.0 , 0.4 ], # dark blue
                [ 0.6 , 0.0 , 0.8 ]  # violet
            ] )

        elif cname == 'sigma0':
            M = nmp.array( [
                [ 0.2 , 0.3 , 0.1 ], # dark redish brown
                [ 1.0 , 0.0 , 0.0 ], # red
                [ 1.0 , 1.0 , 0.0 ], # yellow
                [ 0.2 , 1.0 , 0.0 ], # green
                [ 0.1 , 0.5 , 1.0 ], # light blue
                [ 0.0 , 0.0 , 0.4 ], # dark blue
                [ 0.6 , 0.0 , 0.8 ], # violet
                [ 1.0 , 1.0 , 1.0 ]  # white
            ] )

        elif cname == 'mask':
            M = nmp.array( [
                [ 0.5 , 0.5 , 0.5 ], # gray
                [ 0.5 , 0.5 , 0.5 ]  # gray
            ] )

        else:
            print('ERROR: (''clprn_colmap.py) => unknown "climporn" colormap: '+cname)
            sys.exit(0)

        if lrev:
            # reverse colormap:
            my_cmap = __build_colormap__(M[::-1,:], log_ctrl=log_ctrl, exp_ctrl=exp_ctrl)
        else:
            my_cmap = __build_colormap__(M, log_ctrl=log_ctrl, exp_ctrl=exp_ctrl)

        return my_cmap

#=======================================================================

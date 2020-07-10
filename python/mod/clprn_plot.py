# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
########################################################################
#
#   B a r a K u d a
#
#   Ploting functions and utilities
#
## Authors:
#  --------
# 2010-2015: Laurent Brodeau (original primitive code)
#      2016: Saeed Falahat   (update to fancy grown-up coding!) :D
#
#######################################################################

import os
import sys
import numpy as nmp

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

reload(sys)
sys.setdefaultencoding('utf8')

# Some defaults:
#WDTH_DEF     = 10.
#HGHT_DEF     =  4.
WDTH_DEF     = 9.
HGHT_DEF     = 3.6
FIG_SIZE_DEF = ( WDTH_DEF , HGHT_DEF )
RAT_XY       = WDTH_DEF/10.
DPI_DEF      = 120
AXES_DEF     = [0.09, 0.082, 0.89, 0.84]

# Colors for line:    (https://www.daftlogic.com/projects-hex-colour-tester.htm)
b_blu = '#2C558A'
b_red = '#AD0000'
b_gre = '#3DE079'   ; #b_gre = '#3A783E' ; # b_gre = '#42BD82'
b_prp = '#8C008C'
b_org = '#ED7C4C'
clr_inf_box='0.9'


v_dflt_colors = [b_blu, b_red, b_gre, b_org, b_prp, 'pink', '0.5', 'b', 'g', 'brown', 'orange',
                 '0.25','0.75','k' ]
nmax_colors = len(v_dflt_colors)

# Some projections to use with BaseMap:
#
#         zone PROJ llcrnrlon llcrnrlat urcrnrlon urcrnrlat   lat1  lon0 mer/par continent-res
#                               (lcc = Lambert conformal conic)
projection_def = [
         ['nseas',   'lcc',  -55., 40., 55., 75.,           60., -20., 10., 'l' ],   # Nordic seas
         ['natarct', 'lcc', -60., 40., 80., 72.,            55., -32., 10., 'l' ],   # NATL + Arctic
         ['labir',   'lcc',  -62., 48., -10., 75.,          50., -30.,  5., 'l' ],
         ['labsp',   'lcc',  -60., 48., 50., 75.5,          50., -30., 10., 'l' ],
         ['npol',    'stere', -75., 45., 100., 60.,          80., -30., 10., 'l' ],
         ['npol2',   'stere', -55., 40., 145., 40.,          80.,  -5., 10., 'l' ],   # North Pole
         ['spstere', 'stere',  0.,  0.,  0., 0.,            -48.,  90., 10., 'l' ],   # South Pole Default matplotlib!
         ['matl' ,   'cyl',  -82.,-21.,  12., 79.,          30., -30., 15., 'l' ],   # Nordic seas
         ['atmed',   'lcc',  -18., 33.,  -2., 42.,          30., -10.,  5., 'h' ],
         ['kav7' ,   'kav',    0.,  0.,   0.,  0.,           0.,   0.,  0., 'l' ] ] # global map-monde







class plot :
    ''' This class encapsulates all the plot routines
    In order to use it you need to type as follows:

    plot(function name without the prefix __) (all the arguments of the function)
    for example for __vert_section we do as follows:

    plot("vert_section")(VX, VZ, XF, XMSK, rmin, rmax, dc, lkcont=True, cpal='jet',
        xmin=-80., xmax=85., dx=5, cfignm='fig', cbunit='', cxunit=' ',
        zmin = 0., zmax = 5000., l_zlog=False, cfig_type='png',
        czunit=' ', ctitle=' ', lforce_lim=False, i_cb_subsamp=1, l_z_increase=False )

    The reason to prefix all the function with double underscore __ is that all these function become private members of
    the plot class and they can not be accessed directly outside of the class. That is, you need to call these functios through the call wrapper
    as seen below in the function __call__.
    '''

    def __init__(self,splot) :

        self.splot = splot


    def __call__(self,*args, **kw) :

        if "_"+self.__class__.__name__+ "__" + self.splot in self.__class__.__dict__.keys() :

            self.__class__.__dict__["_"+self.__class__.__name__+ "__" + self.splot](self,*args, **kw)

        else :
            print "function " + "__" + self.splot + " does not exist"
            sys.exit()



    # Functions:


    def __vert_section(self, VX, VZ, XF, XMSK, rmin, rmax, dc, lkcont=False, cpal='jet', lzonal=True,
                       xmin=-80., xmax=85., dx=5, cfignm='fig', cbunit='', cxunit=' ',
                       zmin = 0., zmax = 5000., l_zlog=False, cfig_type='png',
                       czunit=' ', ctitle=' ', lforce_lim=False, i_cb_subsamp=1, l_z_increase=False ):

        import clprn_colmap as bcm

        font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=DPI_DEF)

        zVZ = __prepare_z_log_axis__(l_zlog, VZ)

        if lforce_lim: __force_min_and_max__(rmin, rmax, XF)

        i_x_sbsmp = 1
        dxgap = nmp.amax(VX) - nmp.amin(VX)
        if dxgap > 140.: dx = 5.; i_x_sbsmp = 2
        if dxgap < 50.:  dx = 2.
        if dxgap < 25.:  dx = 1.

        cbgcol='w'
        if not lkcont:
            cbgcol='k'
            XF = nmp.ma.masked_where(XMSK == 0, XF) ; # Masking where mask is zero!

        fig = plt.figure(num = 1, figsize=(WDTH_DEF , RAT_XY*5.), dpi=None, facecolor='w', edgecolor='k')
        ax  = plt.axes([0.07, 0.06, 0.98, 0.88], facecolor=cbgcol)
        vc  = __vcontour__(rmin, rmax, dc)

        # Colormap:
        colmap = bcm.chose_colmap(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)


        if lkcont:
            from clprn_tool import drown
            Xtmp = nmp.zeros(nmp.shape(XF))
            Xtmp[:,:] = XF[:,:]
            drown(Xtmp, XMSK, k_ew=2, nb_max_inc=5, nb_smooth=5)
            cf = plt.contourf(VX, zVZ, Xtmp, vc, cmap=colmap, norm=pal_norm, zorder=0.1)
            plt.contour(      VX, zVZ, Xtmp, vc, colors='k', linewidths=0.2, zorder=0.5)
            del Xtmp
        else:
            cf = plt.pcolormesh(VX, zVZ, XF, cmap=colmap, norm=pal_norm)

        # Colorbar:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, lkc=lkcont, cunit=cbunit, cfont=font_clb, fontsize=10)

        # Masking "rock":
        if lkcont:
            pal_lsm = bcm.chose_colmap('mask')
            norm_lsm = colors.Normalize(vmin = 0., vmax = 1., clip = False)
            prock = nmp.ma.masked_where(XMSK > 0.5, XMSK)
            cm = plt.pcolormesh(VX, zVZ, prock, cmap=pal_lsm, norm=norm_lsm, zorder=1)

        # X-axis:
        if lzonal:
            __nice_longitude_axis__(ax, plt, xmin, xmax, dx*i_x_sbsmp, axt='x')
        else:
            __nice_latitude_axis__(ax, plt, xmin, xmax, dx*i_x_sbsmp, axt='x')

        # Depth-axis:
        __nice_depth_axis__(ax, plt, zmin, zmax, l_log=l_zlog, l_z_inc=l_z_increase, cunit=czunit, cfont=font_xylb)

        plt.title(ctitle, **font_ttl)
        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=False) ; #vert_section
        print '  => '+cfignm+'.'+cfig_type+' created!'
        plt.close(1)

        return



    def __2d(self,VX, VY, XF, XMSK, rmin, rmax, dc, corca='ORCA1', lkcont=True, cpal='jet',
             cfignm='fig', cbunit='', ctitle=' ', lforce_lim=False, i_cb_subsamp=1, cb_orient='vertical',
             cfig_type='pdf', lat_min=-75., lat_max=75., lpix=False, vcont_spec = []):

        #
        # Plot nicely a field given on ORCA coordinates on 2D world map without using any projection
        #
        # if VX = [0] and VY = [0] => ignoring lon and lat...

        import clprn_tool   as bt
        import clprn_orca   as bo
        import clprn_colmap as bcm


        font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=DPI_DEF)

        i_lat_lon = 1
        if len(VX) == 1 or len(VY) == 1: i_lat_lon = 0 ; # no long. and lat. provided !
        if corca[:5] == 'eORCA':
            # don't know how to lat-lon 2d plot eORCA
            # so just plot without projections
            i_lat_lon = 0


        # First drowning the field:
        if not lpix:
            # Don't want to modify XF array, working with XFtmp:
            [ny, nx] = nmp.shape(XF)
            XFtmp = nmp.zeros((ny,nx))
            XFtmp[:,:] = XF[:,:]
            bt.drown(XFtmp, XMSK, k_ew=2, nb_max_inc=20, nb_smooth=10)
        else:
            XFtmp = XF

        ilon_ext = 32

        if lforce_lim: __force_min_and_max__(rmin, rmax, XFtmp)

        XMSK0 = bo.lon_reorg_orca(XMSK,  VX, ilon_ext=ilon_ext)
        XF0   = bo.lon_reorg_orca(XFtmp, VX, ilon_ext=ilon_ext)

        if i_lat_lon == 1:
            [ny, nx] = nmp.shape(XF0)
            dlong  = abs(VX[11] - VX[10])
            VX0 = nmp.arange(0.,nx)
            VX0 = VX0*dlong + dlong/2.

        if i_lat_lon == 1:
            vert_rat = (lat_max - lat_min)/(75. + 75.)
            fig_size = (WDTH_DEF , RAT_XY*4.76*vert_rat) ; #lolo 4.76 => 1080x520 when on 77S->75N
        else:
            fig_size = (WDTH_DEF , RAT_XY*float(nx)/float(ny)*5.)


        # FIGURE
        # ~~~~~~
        fig = plt.figure(num = 1, figsize=fig_size, dpi=None, facecolor='w', edgecolor='k')
        ax  = plt.axes([0.05, 0.06, 1., 0.86], facecolor = '0.5')

        vc = __vcontour__(rmin, rmax, dc)

        # Colmap:
        colmap = bcm.chose_colmap(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        if lpix:
            # Pixelized plot:
            XF0 = nmp.ma.masked_where(XMSK0 == 0, XF0)
            if i_lat_lon == 1:
                cf = plt.pcolormesh(VX0, VY, XF0, cmap = colmap, norm = pal_norm)
            else:
                cf = plt.imshow(             XF0, cmap = colmap, norm = pal_norm)
        else:
            # Contour fill plot:
            if i_lat_lon == 1:
                cf = plt.contourf(VX0, VY, XF0, vc, cmap = colmap, norm = pal_norm)
            else:
                cf = plt.contourf(         XF0, vc, cmap = colmap, norm = pal_norm)

            for c in cf.collections: c.set_zorder(0.15)

            if lkcont:
                if i_lat_lon == 1:
                    cfk = plt.contour(VX0, VY, XF0, vc, colors='k', linewidths = 0.2)
                else:
                    cfk = plt.contour(         XF0, vc, colors='k', linewidths = 0.2)

                for c in cfk.collections: c.set_zorder(0.25)

            # contour for specific values on the ploted field:
            if len(vcont_spec) >= 1:
                if i_lat_lon == 1:
                    cfs = plt.contour(VX0, VY, XF0, vcont_spec, colors='black', linewidths = 1.5)
                else:
                    cfs = plt.contour(         XF0, vcont_spec, colors='black', linewidths = 1.5)

                plt.clabel(cfs, inline=1, fmt='%4.1f', fontsize=10)
                for c in cfs.collections: c.set_zorder(0.35)

        if not lpix:
            # Putting land-sea mask on top of current plot, cleaner than initial masking...
            # because won't influence contours since they are done
            # field needs to be DROWNED prior to this though!!!
            idx_land = nmp.where(XMSK0[:,:] < 0.5)
            XF0 = nmp.ma.masked_where(XMSK0[:,:] > 0.5, XF0)
            XF0[idx_land] = 1000.
            if i_lat_lon == 1:
                cf0 = plt.pcolormesh(VX0, VY, XF0, cmap=bcm.chose_colmap("mask"))
            else:
                cf0 = plt.imshow(             XF0, cmap=bcm.chose_colmap("mask"))

        # Colorbar:
        ifsize = 14.*100./float(DPI_DEF)
        if i_lat_lon == 1: ifsize = int(ifsize*vert_rat); ifsize=max(ifsize,6)
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, lkc=lkcont,
                          cb_or=cb_orient, cunit=cbunit, cfont=font_clb, fontsize=ifsize)

        # X and Y nice ticks:
        if i_lat_lon == 1:
            [vvx, vvy, clon, clat] = __name_coor_ticks__(lon_ext=ilon_ext);
            plt.yticks(vvy,clat) ; plt.xticks(vvx,clon)
            plt.axis([ 0., 360.+ilon_ext-2., lat_min, lat_max])
        else:
            #ax.set_xlim(0., 360.+ilon_ext-2.)
            plt.axis([ 0., float(nx)+ilon_ext-2., 0, float(ny)])

        plt.title(ctitle, **font_ttl)

        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=False) ; #2d

        print '  => '+cfignm+'.'+cfig_type+' created!'
        plt.close(1)

        del XFtmp, XF0

        return



    def __2d_reg(self,VX, VY, XF, XMSK, rmin, rmax, dc, lkcont=False, cpal='jet',
                 cfignm='fig', cfig_type='pdf', cbunit=' ', ctitle='',
                 cb_orient='vertical', lat_min=-77., lat_max=77., i_cb_subsamp=1,
                 lpix=False, l_continent_pixel=True, colorbar_fs=14,
                 col_min='k', col_max='k', vcont_spec = []):

        import clprn_tool   as bt
        import clprn_colmap as bcm

        font_ttl, font_xylb, font_clb, font_inf = __font_unity__()


        # Don't want to modify XF array, working with XFtmp:
        [ny, nx] = nmp.shape(XF)
        XFtmp = nmp.zeros(ny*nx) ; XFtmp.shape = [ny, nx]
        XFtmp[:,:] = XF[:,:]

        # First drowning the field:
        bt.drown(XFtmp, XMSK, k_ew=0, nb_max_inc=20, nb_smooth=10)


        iskp = 28 ; iext = 32

        # Extending / longitude:
        VXe   = bt.extend_domain(VX,    iext, skp_west_deg=iskp) ; nxe = len(VXe)
        XFe   = bt.extend_domain(XFtmp, iext, skp_west_deg=iskp)
        XMSKe = bt.extend_domain(XMSK,  iext, skp_west_deg=iskp)



        # FIGURE

        rat_vert = 1. / ( ( 77. + 77. ) / ( lat_max - lat_min ) )


        if cb_orient == 'horizontal':
            # Horizontal colorbar!
            if ctitle == '':
                fig = plt.figure(num = 1, figsize=(12.4,7.*rat_vert), dpi=None, facecolor='w', edgecolor='k')
                ax = plt.axes([0.05, -0.01, 0.93, 1.], facecolor = 'white')
            else:
                fig = plt.figure(num = 1, figsize=(12.4,7.4*rat_vert), dpi=None, facecolor='w', edgecolor='k')
                ax = plt.axes([0.05, -0.01, 0.93, 0.96], facecolor = 'white')
        else:
            # Vertical colorbar!
            fig = plt.figure(num = 1, figsize=(12.4,6.*rat_vert), dpi=None, facecolor='w', edgecolor='k')
            ax = plt.axes([0.046, 0.06, 1.02, 0.88], facecolor = 'white')

        vc = __vcontour__(rmin, rmax, dc)

        # Colmap:
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)
        mpl.rcParams['contour.negative_linestyle'] = 'solid'
        plt.contour.negative_linestyle='solid'

        if lpix:
            cf = plt.pcolormesh(VXe, VY, XFe, cmap = cpal, norm = pal_norm)
        else:
            cf = plt.contourf(VXe, VY, XFe, vc, cmap = cpal, norm = pal_norm, extend="both")
            for c in cf.collections: c.set_zorder(0.15)
            cf.cmap.set_under(col_min)
            cf.cmap.set_over(col_max)

        # contour for specific values on the ploted field:
        if len(vcont_spec) >= 1:
            cfs = plt.contour(VXe, VY, XFe, vcont_spec, colors='w', linewidths = 1.)
            #plt.clabel(cfs, inline=1, fmt='%4.1f', fontsize=12)

        if lkcont:
            cfk = plt.contour(VXe, VY, XFe, vc, colors='k', linewidths = 0.2)
            for c in cfk.collections: c.set_zorder(0.25)

        # Putting land-sea mask on top of current plot, cleaner than initial masking...
        # because won't influence contours since they are done
        # field needs to be DROWNED prior to this though!!!

        if l_continent_pixel:
            idx_land = nmp.where(XMSKe[:,:] < 0.5)
            XFe = nmp.ma.masked_where(XMSKe[:,:] > 0.5, XFe)
            XFe[idx_land] = 1000.
            cf0 = plt.pcolor(VXe, VY, XFe, cmap = bcm.chose_colmap("mask"))
        else:
            # Masking with contour rather than pixel:
            cf0 = plt.contourf(VXe, VY, XMSKe, [ 0., 0.1 ], cmap = bcm.chose_colmap("mask"))
            for c in cf0.collections: c.set_zorder(5)
            plt.contour(VXe, VY, XMSKe, [ 0.25 ], colors='k', linewidths = 1.)

        # COLOR BAR:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, lkc=lkcont, cb_or=cb_orient, cunit=cbunit, cfont=font_clb, fontsize=colorbar_fs)

        # X and Y nice ticks:
        print "VXe[0], VXe[nxe-1] =>", VXe[0], VXe[nxe-1]
        rlon_min = round(VXe[0],0) ; rlon_max = round(VXe[nxe-1],0)
        print "rlon_min, rlon_max =>", rlon_min, rlon_max

        [vvx, vvy, clon, clat] = __name_coor_ticks__(lon_min=rlon_min, lon_max=rlon_max, dlon=30., lon_ext=iext-iskp)
        plt.yticks(vvy,clat) ; plt.xticks(vvx,clon)

        plt.axis([ rlon_min, rlon_max, lat_min, lat_max])

        if ctitle != ' ': plt.title(ctitle, **font_ttl)

        plt.savefig(cfignm+'.'+cfig_type, dpi=110, orientation='portrait', transparent=False)

        print '  => '+cfignm+'.'+cfig_type+' created!'
        plt.close(1)

        return


    def __2d_box(self,XF, XMSK, rmin, rmax, dc, lkcont=True,
                 cpal='jet', cfignm='fig', cbunit='', ctitle=' ', lforce_lim=False,
                 i_cb_subsamp=1, cfig_type='pdf', lcontours=True,
                 x_offset=0., y_offset=0., vcont_spec = [], lcont_mask=False):

        import clprn_colmap as bcm



        if lforce_lim: __force_min_and_max__(rmin, rmax, XF)

        [ ny , nx ] = XF.shape
        vert_rat = float(ny)/float(nx)
        print "Vert. ratio, nx, ny =", vert_rat, nx, ny

        # Masking field:
        if lcontours:
            idxm = nmp.where(XMSK[:,:] == 0); XF[idxm] = -9999.9  # c'est NaN qui merde!!!
        else:
            XF = nmp.ma.masked_where(XMSK == 0, XF)


        font_ttl, font_xylb, font_clb, font_inf = __font_unity__()


        # FIGURE
        # ~~~~~~
        fig = plt.figure(num = 1, figsize=(7.,6.*vert_rat), dpi=None, facecolor='w', edgecolor='k')

        ax = plt.axes([0.07, 0.05, 0.9, 0.9], facecolor = 'gray')

        vc = __vcontour__(rmin, rmax, dc); #print vc, '\n'

        # Colmap:
        colmap = bcm.chose_colmap(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        if lcontours:
            cf = plt.contourf(XF, vc, cmap = colmap, norm = pal_norm)
            for c in cf.collections: c.set_zorder(0.5)
        else:
            cf = plt.pcolor(XF, cmap = colmap, norm = pal_norm)




        # contour for specific values on the ploted field:
        if len(vcont_spec) >= 1:
            cfs = plt.contour(XF, vcont_spec, colors='white', linewidths = 1.)
            plt.clabel(cfs, inline=1, fmt='%4.1f', fontsize=10)



        if lkcont:
            cfk = plt.contour(XF, vc, colors='k', linewidths = 0.1)
            for c in cfk.collections: c.set_zorder(0.75)



        # contour for continents:
        if lcontours and lcont_mask:
            cfm = plt.contour(XMSK, [ 0.7 ], colors='k', linewidths = 1.)
            for c in cfm.collections: c.set_zorder(1.)

        # Colorbar:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, lkc=lkcont, cunit=cbunit, cfont=font_clb)

        if x_offset > 0 or y_offset > 0 :  __add_xy_offset__(plt, x_offset, y_offset)

        plt.axis([ 0., nx-1, 0, ny-1])

        plt.title(ctitle, **font_ttl)

        # Prevents from using scientific notations in axess ticks numbering:
        ax.get_xaxis().get_major_formatter().set_useOffset(False)

        plt.savefig(cfignm+'.'+cfig_type, dpi=100, orientation='portrait', transparent=True)

        print '  => '+cfignm+'.'+cfig_type+' created!'
        plt.close(1)

        return





    def __zonal(self,VYn, VZn, VY1=[0.],VZ1=[0.], VY2=[0.],VZ2=[0.], VY3=[0.],VZ3=[0.],
                cfignm='fig_zonal', zmin=-100., zmax=100., dz=25., i_z_jump=1,
                xmin=-90., xmax=90., dx=15., cfig_type='png', cxunit=r'Latitude ($^{\circ}$N)',
                czunit='', ctitle='', lab='', lab1='', lab2='', lab3='', box_legend=(0.6, 0.75),
                loc_legend='lower center', fig_size=FIG_SIZE_DEF):

        font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=DPI_DEF)

        ny = len(VYn)
        if len(VZn) != ny: print 'ERROR: plot_zonal.clprn_plot => VYn and VZn do not agree in size'; sys.exit(0)

        lp1=False ; lp2=False ; lp3=False
        if len(VZ1) > 1 and len(VZ1)==len(VY1): lp1=True
        if len(VZ2) > 1 and len(VZ2)==len(VY2): lp2=True
        if len(VZ3) > 1 and len(VZ3)==len(VY3): lp3=True

        if fig_size==FIG_SIZE_DEF: fig_size = (fig_size[0], 1.5*fig_size[1]) # extend height if == to default

        # Do we put the legend outside of the plot?
        l_legend_out = False ; y_leg = 0.
        if loc_legend == 'out':
            l_legend_out = True
            y_leg = 0.1 ; # Figure needs to be vertically extended in that case
            fig_size = (fig_size[0],(1.+y_leg)*fig_size[1])

        fig = plt.figure(num = 1, figsize=fig_size, facecolor='w', edgecolor='k')
        ax  = plt.axes([0.08, 0.075, 0.9, 0.85])

        plt.plot(VYn, VZn*0.0, 'k', linewidth=1)

        plt.plot(VYn, VZn, 'k', linewidth=3., label=lab)
        if lp1: plt.plot(VY1, VZ1, color=b_red, linewidth=2., label=lab1)
        if lp2: plt.plot(VY2, VZ2, color=b_blu, linewidth=2., label=lab2)
        if lp3: plt.plot(VY3, VZ3, color=b_gre, linewidth=2., label=lab3)

        # X-axis
        __nice_latitude_axis__(ax, plt, xmin, xmax, dx, axt='x')

        # Y-axis:
        __nice_y_axis__(ax, plt, zmin, zmax, dz, i_sbsmp=i_z_jump, cunit=czunit, cfont=font_xylb, dy_minor=0)

        # Legend:
        __fancy_legend__(ax, plt, loc_leg=loc_legend, ylg=y_leg, leg_out=l_legend_out, ncol=1)

        plt.title(ctitle, **font_ttl)

        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=False)
        plt.close(1)

        return




    def __nproj(self,czone, rmin, rmax, dc, xlon, xlat, XF,
                cfignm='fig', lkcont=False, cpal='jet', cbunit=' ',
                cfig_type='pdf', ctitle=' ', lforce_lim=False,
                cb_orient='vertical', i_cb_subsamp=1, dpi_fig=DPI_DEF, lpcont=True):

        # Plot projection with basemap...

        #===================================================================================
        # INPUT:
        #          xlon and xlat can be 1D or 2D !!!
        #
        #   lpcont=True  => do contourf
        #   lpcont=False => do pcolor
        #
        #===================================================================================



        from mpl_toolkits.basemap import Basemap
        from mpl_toolkits.basemap import shiftgrid
        import clprn_colmap as bcm

        font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=dpi_fig)

        # For projections :
        vp = __give_proj__(czone) ; # projection information


        # must work with XFtmp rather than XF, because sometimes XF is overwrited...
        [ny, nx] = nmp.shape(XF)
        XFtmp = nmp.zeros(ny*nx) ; XFtmp.shape = [ny, nx]
        XFtmp[:,:] = XF[:,:]


        if len(nmp.shape(xlat)) == 1 and len(nmp.shape(xlon)) == 1:
            if czone == 'kav7' and xlon[0] >= 0.:
                # Shifting data and longitude to be consistent with map projection
                XFtmp, xlon = shiftgrid(180.+xlon[0], XFtmp, xlon, start=False, cyclic=360.0)
            LON_2D, LAT_2D = nmp.meshgrid(xlon,xlat)
        else:
            LAT_2D = nmp.zeros(ny*nx) ; LAT_2D.shape = [ny, nx] ; LAT_2D[:,:] = xlat[:,:]
            LON_2D = nmp.zeros(ny*nx) ; LON_2D.shape = [ny, nx] ; LON_2D[:,:] = xlon[:,:]


        if lforce_lim: __force_min_and_max__(rmin, rmax, XFtmp)

        vc = __vcontour__(rmin, rmax, dc)

        # Colorbar position/size if horizontal
        vcbar = [0.1, 0.08, 0.86, 0.03]

        # Figure/canvas size:
        if cb_orient == 'horizontal':
            if czone == 'natarct':
                vfig_size = [ 5.8, 5.6 ]; vsporg = [0.08, 0.1, 0.9,  0.92]
                vcbar = [0.05, 0.08, 0.9, 0.03]
            if czone == 'npol2':
                vfig_size = [ 4.4, 5.6 ];  vsporg = [0.01, 0.15, 1., 0.8]
                vcbar = [0.05, 0.065, 0.92, 0.03]
            if czone == 'kav7':
                vfig_size = [ 8.1, 5.6 ];  vsporg = [0.001, 0.15, 1., 0.8]
                vcbar = [0.04, 0.08, 0.92, 0.03]

        else:
            # Vertical color bar on the rhs
            rw = 5.
            vfig_size                        = [ rw, rw ];        vsporg = [0.1, 0.1, 0.85, 0.85]
            if czone == 'nseas':   vfig_size = [ rw , 0.7*rw ];   vsporg = [0.08,  0.07, 0.85, 0.85]
            if czone == 'natarct': vfig_size = [ rw , rw ];       vsporg = [0.065, 0.04, 0.95, 0.92]
            if czone == 'spstere': vfig_size = [ rw , 0.76*rw ];  vsporg = [0.11, 0.05, 0.82, 0.89]
            if czone == 'npol2':   vfig_size = [ rw , 0.96*rw ] ; vsporg = [0.1,  0.05, 0.86, 0.9]


        fig = plt.figure(num = 1, figsize=(vfig_size), dpi=None, facecolor='w', edgecolor='k')
        ax  = plt.axes(vsporg, facecolor = 'w')


        ## Colmap:
        colmap = bcm.chose_colmap(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)
        mpl.rcParams['contour.negative_linestyle'] = 'solid'; plt.contour.negative_linestyle='solid'


        if vp[1] == 'lcc' or vp[1] == 'cyl' :
            carte = Basemap(llcrnrlon=vp[2],llcrnrlat=vp[3],urcrnrlon=vp[4],urcrnrlat=vp[5],\
                            resolution=vp[9],area_thresh=1000.,projection=vp[1],\
                            lat_1=vp[6],lon_0=vp[7])

        elif vp[1] == 'stere' :
            if vp[0] == 'spstere' or vp[0] == 'npstere':
                carte = Basemap(projection=vp[0], boundinglat=vp[6], lon_0=vp[7], resolution=vp[9])
            else:
                carte = Basemap(llcrnrlon=vp[2],llcrnrlat=vp[3],urcrnrlon=vp[4],urcrnrlat=vp[5],\
                              resolution=vp[9],area_thresh=1000.,projection='stere',\
                              lat_0=vp[6],lon_0=vp[7])
        elif vp[1] == 'kav' :
            print ' *** plot_nproj.clprn_plot => Projection '+vp[0]+' / '+str(vp[7])+' / '+vp[9]
            carte = Basemap(projection=vp[0],lon_0=vp[7],resolution=vp[9])

        else:
            print 'ERROR: clprn_plot.py => proj type '+vp[1]+' unknown!!!'; sys.exit(0)

        x0,y0 = carte(LON_2D,LAT_2D)

        if lpcont:
            cf = carte.contourf(x0, y0, XFtmp, vc, cmap = colmap, norm = pal_norm)
            # Black contours if needed :
            if lkcont:
                ckf = carte.contour(x0, y0, XFtmp, vc, colors='k', linewidths=0.5)
                if cpal != 'ice':
                    for c in cf.collections: c.set_zorder(0.5)   # Changing zorder so black cont. on top
                for c in ckf.collections: c.set_zorder(1.) # of filled cont. and under continents (zorder 1)

        else:
            cf = carte.pcolor(x0, y0, XFtmp, cmap = colmap, norm = pal_norm)




        carte.drawcoastlines() ; carte.fillcontinents(color='grey') ; carte.drawmapboundary()


        if vp[1] == 'lcc' or vp[1] == 'cyl':
            carte.drawmeridians(nmp.arange(-360,360,vp[8]), labels=[0,0,0,1])
            carte.drawparallels(nmp.arange( -90, 90,vp[8]), labels=[1,0,0,0])

        if vp[1] == 'stere':
            carte.drawmeridians(nmp.arange(-180,180,20), labels=[0,0,0,1])
            carte.drawparallels(nmp.arange( -90, 90,10), labels=[1,0,0,0])

        plt.title(ctitle, **font_ttl)

        # Colorbar:
        if cb_orient == 'horizontal':
            clbax = fig.add_axes(vcbar) # new axes for colorbar!
            __nice_colorbar__(cf, plt, vc, cax_other=clbax, i_sbsmp=i_cb_subsamp, lkc=(lkcont and lpcont), cb_or='horizontal', cunit=cbunit, cfont=font_clb, fontsize=10)
        else:
            __nice_colorbar__(cf, plt, vc,                  i_sbsmp=i_cb_subsamp, lkc=(lkcont and lpcont),                     cunit=cbunit, cfont=font_clb, fontsize=12)

        plt.savefig(cfignm+'.'+cfig_type, dpi=dpi_fig, orientation='portrait', transparent=False) ; #, transparent=True, acecolor='w', edgecolor='w',trans

        plt.close(1)

        print ' *** created figure '+cfignm+'.'+cfig_type+'\n'

        del LON_2D, LAT_2D, XFtmp

        return



    def __2d_box_2f(self,XF1, XF2, XMSK, rmin, rmax, dc, vcont_spec2, corca='ORCA1', lkcont=True,
                    cpal='jet', cfignm='fig', cbunit='', ctitle=' ', lforce_lim=False,
                    i_cb_subsamp=1, cfig_type='pdf', lcontours=True,
                    x_offset=0., y_offset=0., vcont_spec1 = []):

        # Take 2 fields as imput and shows contours of second field (vcont_spec2) on top of field 1

        import matplotlib.colors as colors   # colmap and co.
        import clprn_colmap as bcm


        if nmp.shape(XF1) != nmp.shape(XF2):
            print 'ERROR clprn_plot.plot_2d_box_2f: fields F1 and F2 dont have the same shape!'
            sys.exit(0)



        font_ttl, font_xylb, font_clb, font_inf = __font_unity__()


        if lforce_lim: __force_min_and_max__(rmin, rmax, XF1)

        [ ny , nx ] = XF1.shape
        vert_rat = float(ny)/float(nx)
        print "Vert. ratio, nx, ny =", vert_rat, nx, ny

        # Masking field:
        if lcontours:
            idxm = nmp.where(XMSK[:,:] == 0); XF1[idxm] = -9999.9  # c'est NaN qui merde!!!
        else:
            XF1 = nmp.ma.masked_where(XMSK == 0, XF1)



        # FIGURE
        # ~~~~~~
        fig = plt.figure(num = 1, figsize=(7.,6.*vert_rat), dpi=None, facecolor='w', edgecolor='k')

        ax = plt.axes([0.07, 0.05, 0.9, 0.9], facecolor = 'gray')

        vc = __vcontour__(rmin, rmax, dc); #print vc, '\n'

        # Colmap:
        colmap = bcm.chose_colmap(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        if lcontours:
            cf = plt.contourf(XF1, vc, cmap = colmap, norm = pal_norm)
            for c in cf.collections: c.set_zorder(0.5)
        else:
            cf = plt.pcolor(XF1, cmap = colmap, norm = pal_norm)

        # contour for specific values on the ploted field:
        if len(vcont_spec1) >= 1:
            cfs1 = plt.contour(XF1, vcont_spec1, colors='white', linewidths = 1.)
            plt.clabel(cfs1, inline=1, fmt='%4.1f', fontsize=10)

        # Contours of field F2:
        cfs2 = plt.contour(XF2, vcont_spec2, colors=b_red, linewidths = 1.3)
        #plt.clabel(cfs1, inline=1, fmt='%4.1f', fontsize=10)


        if lkcont:
            cfk = plt.contour(XF1, vc, colors='k', linewidths = 0.1)
            for c in cfk.collections: c.set_zorder(0.75)





        # contour for continents:
        if lcontours:
            cfm = plt.contour(XMSK, [ 0.7 ], colors='k', linewidths = 0.4)
            for c in cfm.collections: c.set_zorder(1.)

        # Colorbar:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, lkc=lkcont, cunit=cbunit, cfont=font_clb)

        if x_offset > 0 or y_offset > 0 :  __add_xy_offset__(plt, x_offset, y_offset)

        plt.axis([ 0., nx-1, 0, ny-1])

        plt.title(ctitle, **font_ttl)

        plt.savefig(cfignm+'.'+cfig_type, dpi=100, orientation='portrait', transparent=True)

        print '  => '+cfignm+'.'+cfig_type+' created!'
        plt.close(1)

        del Xtmp

        return






    def __trsp_sig_class(self,VT, vsigma_bounds, XF, rmin, rmax, dc,
                         lkcont=True, cpal='sigtr', dt=5., cfignm='fig',
                         cfig_type='pdf', ctitle='', vcont_spec1 = [],
                         i_cb_subsamp=2):

        # Plot transport by sigma class...
        if nmp.sum(XF) == 0.:
            print '\n  WARNING: plot_trsp_sig_class => doing nothing, arrays contains only 0!\n'
        else:
            
            import matplotlib.colors as colors   # colmap and co.
            import clprn_colmap as bcm
    
            font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=DPI_DEF)
    
            fig = plt.figure(num = 1, figsize=(WDTH_DEF , RAT_XY*6.), dpi=None, facecolor='w', edgecolor='k') ; #trsp_sig_class
            ax = plt.axes([0.075,  -0.025, 0.9, 0.98], facecolor = 'w')
    
            vc = __vcontour__(rmin, rmax, dc)
    
            nbins = len(vsigma_bounds) - 1
    
            # Colmap:
            colmap = bcm.chose_colmap(cpal)
            pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)
            mpl.rcParams['contour.negative_linestyle'] = 'solid'
            plt.contour.negative_linestyle='solid'
    
            #cf = plt.contourf(VT, vsigma_bounds[:-1], XF, vc, cmap = colmap, norm = pal_norm)
            cf = plt.pcolormesh(VT, vsigma_bounds[:-1], XF,  cmap = colmap, norm = pal_norm)
            if lkcont:
                cfc = plt.contour(VT, vsigma_bounds[:-1], XF, nmp.arange(-3.,3.,0.5), colors='k', linewidths=0.4)
    
            # contour for specific values on the ploted field:
            if len(vcont_spec1) >= 1:
                cfs1 = plt.contour(VT, vsigma_bounds[:-1], XF, vcont_spec1, colors='white', linewidths = 1.)
                plt.clabel(cfs1, inline=1, fmt='%4.1f', fontsize=11, manual=[(2080,2.)] )
    
            __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, cb_or='horizontal', cunit='Sv', cfont=font_clb, fontsize=10)
    
            # AXES:
            x1 = int(min(VT))  ; x2 = int(max(VT))+1
            plt.axis([x1, x2, vsigma_bounds[nbins], vsigma_bounds[0]])
    
            __nice_x_axis__(ax, plt, x1, x2, dt, cfont=font_xylb, dx_minor=__time_axis_minor_ticks__(dt))
    
            plt.yticks( nmp.flipud(vsigma_bounds) )
    
            label_big = { 'fontname':'Trebuchet MS', 'fontweight':'normal', 'fontsize':18 }
            plt.ylabel(r'$\sigma_0$', **label_big)
    
            plt.title(ctitle, **font_ttl)
            plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=False) ; #trsp_sig_class
            print '  => '+cfignm+'.'+cfig_type+' created!'
            plt.close(1)

        return




    def __vert_section_extra(self,VX, VZ, XF, XMSK, Vcurve, rmin, rmax, dc, lkcont=True, cpal='jet', xmin=-80., xmax=85.,
                             cfignm='fig', cbunit='', cxunit=' ', zmin = 0., zmax = 5000., l_zlog=False,
                             cfig_type='pdf', czunit=' ', ctitle=' ', lforce_lim=False, fig_size=(8.,8.) ):

        import matplotlib.colors as colors   # colmap and co.
        import clprn_colmap as bcm

        zVZ = __prepare_z_log_axis__(l_zlog, VZ)

        XF = nmp.ma.masked_where(XMSK == 0, XF)

        if lforce_lim: __force_min_and_max__(rmin, rmax, XF)
        # 
        font_ttl, font_xylb, font_clb, font_inf = __font_unity__()


        fig = plt.figure(num = 1, figsize=fig_size, dpi=None, facecolor='w', edgecolor='k')
        ax = plt.axes([0.1,  0.065,   0.92,       0.89], facecolor = 'gray')
        vc = __vcontour__(rmin, rmax, dc)

        # Colmap:
        colmap = bcm.chose_colmap(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        cf = plt.contourf(VX, zVZ, XF, vc, cmap = colmap, norm = pal_norm)
        if lkcont: plt.contour(VX, zVZ, XF, vc, colors='k', linewidths=0.2)

        # Colorbar:
        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, cunit=cbunit, cfont=font_clb, fontsize=10)

        # X-axis:
        __nice_x_axis__(ax, plt, xmin, xmax, dx, cunit=cxunit, cfont=font_xylb)

        plt.plot(VX,Vcurve, 'w', linewidth=2)

        for zz in zVZ[:]: plt.plot(VX,VX*0.+zz, 'k', linewidth=0.3)

        # Depth axis:
        __nice_depth_axis__(ax, plt, zmin, zmax, l_log=l_zlog, l_z_inc=False, cunit=czunit, cfont=font_xylb)

        plt.title(ctitle, **font_ttl)
        plt.savefig(cfignm+'.'+cfig_type, dpi=100, orientation='portrait', transparent=True)
        print '  => '+cfignm+'.'+cfig_type+' created!'
        plt.close(1)
        #
        return




    def __hovmoeller(self, VT, VY, XF, XMSK, rmin, rmax, dc, c_y_is='depth',
                     lkcont=False, cpal='jet', tmin=0., tmax=50., dt=5,
                     ymin=0., ymax=5000., dy=100., l_ylog=False,
                     cfignm='fig', cbunit='', ctunit=' ', cfig_type='png',
                     cyunit=' ', ctitle=' ', i_cb_subsamp=1,
                     l_y_increase=False ):
        #
        # c_y_is : 'depth', 'latitude'
        # lkcont : use contours rather than "pcolormesh"
        #

        import clprn_colmap as bcm

        font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=DPI_DEF)

        if c_y_is == 'depth':
            zVY = __prepare_z_log_axis__(l_ylog, VY)
            vax = [0.095, 0.06, 0.92, 0.88]
        else:
            zVY = VY
            vax = [0.05, 0.06, 0.98, 0.88]

        # Masking where mask is zero!
        XF = nmp.ma.masked_where(XMSK == 0, XF)

        fig = plt.figure(num = 1, figsize=(WDTH_DEF , RAT_XY*5.), dpi=None, facecolor='w', edgecolor='k')
        ax  = plt.axes(vax, facecolor='gray')
        vc  = __vcontour__(rmin, rmax, dc)

        # Colormap:
        colmap = bcm.chose_colmap(cpal)
        pal_norm = colors.Normalize(vmin = rmin, vmax = rmax, clip = False)

        if lkcont:
            cf = plt.contourf(VT, zVY, XF, vc, cmap = colmap, norm = pal_norm)
            #plt.contour( VT, zVY, XF, vc, colors='k', linewidths=0.2)
        else:
            cf = plt.pcolormesh(VT, zVY, XF, cmap = colmap, norm = pal_norm)

        __nice_colorbar__(cf, plt, vc, i_sbsmp=i_cb_subsamp, cunit=cbunit, cfont=font_clb, fontsize=10)

        # Time-axis:
        __nice_x_axis__(ax, plt, tmin, tmax, dt, cunit=ctunit, cfont=font_xylb, dx_minor=__time_axis_minor_ticks__(dt))

        # Y-axis:
        if c_y_is == 'depth':
            __nice_depth_axis__(ax, plt, ymin, ymax, l_log=l_ylog, l_z_inc=l_y_increase, cunit=cyunit, cfont=font_xylb)
        elif c_y_is == 'latitude':
            __nice_latitude_axis__(ax, plt, ymin, ymax, dy, axt='y')
        else:
            print 'ERROR: plot_hoevmoller.clprn_plot => axis "'+c_y_is+'" not supported!'; sys.exit(0)

        plt.title(ctitle, **font_ttl)
        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=False) ; #vert_section
        print '  => '+cfignm+'.'+cfig_type+' created!'
        plt.close(1)

        return



    def __oscillation_index(self, VT, VF, ymax=2.5, dy=0.5, yplusminus=0.,
                            tmin=0., tmax=0., dt=5,
                            cfignm='fig', cfig_type='png', cyunit='', ctitle=''):

        #--------------------------------------------------------------------------------------
        # Plot a ENSO / AMO / PDO -like graph from a time series VF that
        # has already been smoothed and detrended
        #--------------------------------------------------------------------------------------

        font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=DPI_DEF)

        Nt = len(VT)
        if len(VF) != Nt: print 'ERROR: oscillation_index.clprn_plot => VT and VF do not agree in size'; sys.exit(0)

        vf_plus = nmp.zeros(Nt) ; vf_mins = nmp.zeros(Nt)
        vf_plus[:] = VF[:]   ; vf_mins[:] = VF[:]
        vf_plus[nmp.where(VF[:] < 0. )] = 0.
        vf_mins[nmp.where(VF[:] > 0. )] = 0.
        vf_plus[0]  = 0. ; vf_mins[0]  = 0.
        vf_plus[-1] = 0. ; vf_mins[-1] = 0.

        t1 = tmin ; t2 = tmax
        if tmin == 0. and tmax == 0.:
            t1 = float(int(min(VT)))
            t2 = float(int(round(max(VT),0)))

        fig = plt.figure(num = 2, figsize=FIG_SIZE_DEF, facecolor='w', edgecolor='k')
        ax  = plt.axes(AXES_DEF)

        if yplusminus > 0.:
            plt.plot(VT, 0.*VT+yplusminus, 'r--', linewidth=1.5)
            plt.plot(VT, 0.*VT-yplusminus, 'b--', linewidth=1.5)

        plt.fill(VT, vf_plus, b_red, VT, vf_mins, b_blu, linewidth=0)
        plt.plot(VT, VF[:], 'k', linewidth=0.7)
        plt.plot(VT, 0.*VT, 'k', linewidth=0.7)

        __nice_x_axis__(ax, plt,    t1,   t2, dt,               cfont=font_xylb)
        __nice_y_axis__(ax, plt, -ymax, ymax, dy, cunit=cyunit, cfont=font_xylb)

        plt.title(ctitle, **font_ttl)
        plt.savefig(cfignm+'.'+cfig_type, dpi=DPI_DEF, orientation='portrait', transparent=True)
        plt.close(2)
        return



    def __1d_mon_ann(self,VTm, VTy, VDm, VDy, cfignm='fig', dt=5, cyunit='', ctitle='',
                     ymin=0, ymax=0, dy=0, i_y_jump=1, mnth_col='b', plt_m03=False, plt_m09=False,
                     cfig_type='png', l_tranparent_bg=True, fig_size=FIG_SIZE_DEF, y_cst_to_add=-9999.):

        # if you specify ymin and ymax you can also specify y increment (for y grid) as dy
        #
        # plt_m03 => plot march values on top in green
        # plt_m09 => plot september values on top in green

        font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=DPI_DEF)

        Nt1 = len(VTm) ; Nt2 = len(VTy)

        if len(VTm) != len(VDm): print 'ERROR: plot_1d_mon_ann.clprn_plot => VTm and VDm do not agree in size'; sys.exit(0)
        if len(VTy) != len(VDy): print 'ERROR: plot_1d_mon_ann.clprn_plot => VTy and VDy do not agree in size'; sys.exit(0)

        l_add_monthly = True
        if Nt1 == Nt2: l_add_monthly = False
        
        y_leg = 0.
        if plt_m03 or plt_m09:
            # We put the legend outside of the plot...
            y_leg = 0.075*2. ; # Figure needs to be vertically extended in that case
            fig_size = (fig_size[0],(1.+y_leg)*fig_size[1])

        fig = plt.figure(num = 1, figsize=fig_size, facecolor='w', edgecolor='k')

        ax = plt.axes(AXES_DEF) ; #1d_mon_ann

        if mnth_col == 'g': mnth_col = b_gre
        if mnth_col == 'b': mnth_col = b_blu

        if y_cst_to_add > -9000.:
            plt.plot(VTm, VTm*0.+y_cst_to_add, 'k', label=None, linewidth=1.8)
        if l_add_monthly:
            plt.plot(VTm, VDm, mnth_col, label=r'monthly', linewidth=1)
        plt.plot(VTy, VDy, b_red, label=r'annual', linewidth=2)

        ax.get_yaxis().get_major_formatter().set_useOffset(False); # Prevents from using scientific notations in axess ticks numbering

        if l_add_monthly:
            if plt_m03: plt.plot(VTm[2:Nt1:12], VDm[2:Nt1:12], 'orange', label=r'March',     linewidth=2)
            if plt_m09: plt.plot(VTm[8:Nt1:12], VDm[8:Nt1:12], 'orange', label=r'September', linewidth=2)
            if plt_m03 or plt_m09:
                box = ax.get_position()
                ax.set_position([box.x0, box.y0 + box.height*y_leg, box.width, box.height*(1.-y_leg)])
                plt.legend(bbox_to_anchor=(0.95, -0.075), ncol=2, shadow=True, fancybox=True)

        # Time bounds for t-axis:
        dcorr = 0.
        x1 = float(int(min(VTy)))
        #x2 = float(int(max(VTy)))
        x2 = round(max(VTy),0)
        #if min(VTy)-x1==0.5:
        #x1 = int(min(VTy)-dcorr)
        #x2 = int(max(VTy)+dcorr)

        mean_val = nmp.mean(VDy)
        df = max( abs(min(VDm)-mean_val), abs(max(VDm)-mean_val) )
        
        if ymin==0 and ymax==0:
            y1, y2, dy = __suitable_axis_dx__(min(VDm)-0.2*df, max(VDm)+0.2*df, nb_val=10.)
        elif dy == 0:
            y1, y2, dy = __suitable_axis_dx__(ymin,            ymax,            nb_val=10.)
        else:
            y1=ymin ; y2=ymax

        __nice_y_axis__(ax, plt, y1, y2, dy, i_sbsmp=i_y_jump, cunit=cyunit, cfont=font_xylb, dy_minor=0)

        __nice_x_axis__(ax, plt, x1, x2, dt, cfont=font_xylb, dx_minor=__time_axis_minor_ticks__(dt))

        plt.title(ctitle, **font_ttl)

        cf_fig = cfignm+'.'+cfig_type

        plt.savefig(cf_fig, dpi=DPI_DEF, orientation='portrait', transparent=l_tranparent_bg)
        print '  => '+cfignm+'.'+cfig_type+' created!'
        plt.close(1)








    def __1d_multi(self,vt, XD, vlabels, cfignm='fig', dt=5, i_t_jump=1, cyunit=None, ctitle='',
                   cfig_type='png', ymin=0, ymax=0, lzonal=False, xmin=0, xmax=0,
                   loc_legend='lower center', line_styles=[], fig_size=FIG_SIZE_DEF,
                   l_tranparent_bg=True, cxunit=None, lmask=True, cinfo='', y_cst_to_add=-9999.):

        # lzonal => zonally averaged curves...
        if lzonal:
            font_ttl, font_big_fixed, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=DPI_DEF, size='big')
        else:
            font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=DPI_DEF)

        # Number of lines to plot:
        [ nb_plt, nbt ] = XD.shape
        if len(vt) != nbt: print 'ERROR: plot_1d_multi.clprn_plot.py => vt and XD do not agree in shape! =>', len(vt), nbt,'\n'; sys.exit(0)
        if len(vlabels) != nb_plt: print 'ERROR: plot_1d_multi.clprn_plot.py => wrong number of labels...'; sys.exit(0)
        n0 = len(line_styles)
        if n0 > 0 and n0 != nb_plt: print 'ERROR: plot_1d_multi.clprn_plot.py => wrong number line styles!!!'; sys.exit(0)
        nb_col, nb_row = __nb_col_row_legend__(nb_plt) ; # nb of columns and rows for legend

        # Do we put the legend outside of the plot?
        l_legend_out = False ; y_leg = 0.
        if loc_legend == 'out':
            l_legend_out = True
            y_leg = 0.075*nb_row ; # Figure needs to be vertically extended in that case
            fig_size = (fig_size[0],(1.+y_leg)*fig_size[1])

        # Masking the time-series shorter than others (masked with -999.)
        if lmask: XD = nmp.ma.masked_where(XD < -900., XD)

        if lzonal:
            fig = plt.figure(num = 1, figsize=fig_size, facecolor='w', edgecolor='k')
            ax = plt.axes([0.08, 0.11, 0.88, 0.83])
        else:
            fig = plt.figure(num = 1, figsize=fig_size, facecolor='w', edgecolor='k')
            ax = plt.axes(AXES_DEF) ; #1d_multi

        if y_cst_to_add > -9000.:
            plt.plot(vt, vt*0.+y_cst_to_add, 'k', label=None, linewidth=1.8)

        if lzonal: plt.plot(vt, XD[0,:]*0., 'k', linewidth=1)

        
        if n0 <= 0 and nb_plt > nmax_colors:
            print 'ERROR: plot_1d_multi.clprn_plot => not enough colors defined in "v_dflt_colors", extend it!!!'
            sys.exit(0)

        for jp in range(nb_plt):
            if n0 > 0:
                plt.plot(vt, XD[jp,:], line_styles[jp],   label=vlabels[jp], linewidth=2)
            else:
                plt.plot(vt, XD[jp,:], v_dflt_colors[jp], label=vlabels[jp], linewidth=2)

        #ax.get_yaxis().get_major_formatter().set_useOffset(False) ; # Prevents from using scientific notations in axess ticks numbering

        if lzonal:
            dt = 15. ; # x-axis increment (latitude!)
            if xmin == 0 and xmax == 0:
                x1 = -90. ; x2 = 90.
            else:
                x1 = xmin ;  x2 = xmax
        else:
            if xmin == 0 and xmax == 0:
                x1 = int(vt[0])
                x2 = int(round(vt[len(vt)-1]+0.4))
            else:
                x1 = xmin ; x2 = xmax

        if ymin==0 and ymax==0:
            ymin = nmp.min(XD[:,:])
            ymax = nmp.max(XD[:,:])
        ymin, ymax, dy = __suitable_axis_dx__(ymin, ymax, nb_val=10.)

        if lzonal:
            __nice_x_axis__(ax, plt, x1, x2, 10., cunit=r'Latitude ($^{\circ}$N)', cfont=font_xylb, dx_minor=5.)
        else:
            __nice_x_axis__(ax, plt, x1, x2, dt, i_sbsmp=i_t_jump, cunit=cxunit, cfont=font_xylb, dx_minor=__time_axis_minor_ticks__(dt))

        __nice_y_axis__(ax, plt, ymin, ymax, dy, i_sbsmp=1, cunit=cyunit, cfont=font_xylb, dy_minor=0)

        plt.title(ctitle, **font_ttl)
        
        if cinfo != '':
            # Ading info:
            yp = 0.95
            if loc_legend != '0' and l_legend_out: yp = -0.1
            props = dict(boxstyle='round', facecolor='w') ;#, alpha=0.5)
            ax.text(0.05, yp, cinfo, transform=ax.transAxes,
                    verticalalignment='top', bbox=props, fontsize=10)
        
        __fancy_legend__(ax, plt, loc_leg=loc_legend, ylg=y_leg, leg_out=l_legend_out, ncol=nb_col)

        cf_fig = cfignm+'.'+cfig_type

        plt.savefig(cf_fig, dpi=DPI_DEF, orientation='portrait', transparent=l_tranparent_bg) ; #1d_multi

        plt.close(1)
        print '  => Multi figure "'+cf_fig+'" created!'




    def __1d(self,vt, VF, cfignm='fig', dt=5, i_t_jump=1, cyunit='', ctitle='',
                cfig_type='png', ymin=0, ymax=0, xmin=0, xmax=0,
                loc_legend='lower center', line_styles='-', fig_size=FIG_SIZE_DEF,
                l_tranparent_bg=False, cxunit='', lmask=True):

        font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=DPI_DEF)

        # Number of lines to plot:
        nbt = len(VF)

        if len(vt) != nbt: print 'ERROR: plot_1d.clprn_plot.py => vt and VF do not agree in shape!'; sys.exit(0)


        # Masking the time-series shorter than others (masked with -999.)
        if lmask: VF = nmp.ma.masked_where(VF < -900., VF)

        fig = plt.figure(num = 1, figsize=fig_size, facecolor='w', edgecolor='k')
        ax = plt.axes(AXES_DEF) ; #1d

        plt.plot(vt, VF[:], line_styles, linewidth=2)

        ax.get_yaxis().get_major_formatter().set_useOffset(False) ; # Prevents from using scientific notations in axess ticks numbering

        if xmin == 0 and xmax == 0:
            x1 = int(vt[0])
            x2 = int(round(vt[len(vt)-1]+0.4))
        else:
            x1 = xmin ; x2 = xmax

        if ymin==0 and ymax==0:
            mean_val = nmp.mean(VF[:])
            dA = max( abs(nmp.min(VF[:])-mean_val), abs(nmp.max(VF[:])-mean_val) )
            plt.axis( [x1, x2, nmp.min(VF[:])-0.2*dA, nmp.max(VF[:])+0.2*dA] )
        else:
            plt.axis([x1, x2, ymin,     ymax])

        print nmp.arange(x1, x2+dt, dt)

        __nice_x_axis__(ax, plt, x1, x2, dt, i_sbsmp=i_t_jump, cfont=font_xylb)

        if cyunit != '': plt.ylabel('('+cyunit+')', **font_xylb)
        if cxunit != '': plt.xlabel('('+cxunit+')', **font_xylb)

        plt.title(ctitle, **font_ttl)

        cf_fig = cfignm+'.'+cfig_type

        plt.savefig(cf_fig, dpi=DPI_DEF, orientation='portrait', transparent=l_tranparent_bg) ; #1d

        plt.close(1)
        print '   => Multi figure "'+cf_fig+'" created!'











    def __spectrum(self,vfrq, Vspec, cfignm='fig', cyunit='', log_x=True, log_y=False,
                      year_min=3., year_max = 50., rmax_amp = 10., rmin_amp = 0.,
                      cfig_type='png', vnoise=[ 0 ], vrci95=[ 0 ], lab95_xpos=0.5, lplot_1onF=False,
                      cnoise='White noise', lplot_freq_ax=True):

        l_do_ci95 = False ; l_do_ci95m = False

        nnoise = len(vnoise); nl95 = len(vrci95)

        if nnoise != 1 and nl95 != 1:
            if nl95 != len(Vspec) or nl95 != nnoise:
                print "ERROR: plot_spectrum.clprn_plot.py => length of 95 CI array and/or noise array doesnt match spectrum length!"
                sys.exit(0)
            l_do_ci95 = True
            l_do_ci95m = True

        font_ttl, font_xylb, font_clb, font_inf = __font_unity__()


        print "avant:", rmin_amp, rmax_amp
        if log_y:
            if rmin_amp <= 0.: rmin_amp = 0.01
            rmin_amp = 20.*nmp.log10(rmin_amp); rmax_amp = 20.*nmp.log10(rmax_amp)
        print "apres:", rmin_amp, rmax_amp


        # Spectral axis:
        x_min = 1./year_max ; x_max = 1./year_min ; # min and max in frequency!
        clbnd = str(int(round(year_min)))

        if log_x:
            cvee = [ '50','45','40','35','30','25','22','20','17','15','13','12','11','10','9','8','7','6','5','4','3' ]
            if year_max == 40.: cvee = cvee[2:]
            if year_max == 35.: cvee = cvee[3:]
            if year_min ==  5.: cvee = cvee[:-2]
        else:
            cvee = [ '50','30','20','15','12','10','9','8','7','6','5','4', '3' ]


        lvee = []
        civee = []
        for ce in cvee:
            lvee.append(float(ce))
            civee.append(str(round(1./float(ce),3)))
        vee  = nmp.asarray(lvee)

        print civee[:]


        rnoise = nmp.mean(vnoise[5:20])
        rrci95 = nmp.mean(vrci95[5:20])


        fig = plt.figure(num = 1, figsize=(8.,4.), facecolor='w', edgecolor='k')

        if lplot_freq_ax:
            ax = plt.axes([0.069, 0.13, 0.9, 0.8])
        else:
            ax = plt.axes([0.08, 0.13, 0.9, 0.83])

        if log_x:
            ifr1 = 1
            vfl = nmp.log10(vfrq[ifr1:])

            if not log_y:
                if l_do_ci95:
                    plt.plot(vfl, vnoise[ifr1:],               '--k'  , linewidth=1.8, label=cnoise)
                    plt.plot(vfl, vnoise[ifr1:]+vrci95[ifr1:], '0.4', linewidth=1.8, label='95% CI')
                    plt.plot(vfl, vnoise[ifr1:]-vrci95[ifr1:], '0.4', linewidth=1.8)
                plt.plot(vfl, Vspec[ifr1:],   '*-k', linewidth=2.)
                #if lplot_1onF: plt.plot(vfl, 1./vfl,   b_red, linewidth=2)
            else:
                if l_do_ci95:
                    plt.plot(vfl, 20.*nmp.log10(vnoise[ifr1:]),               '--k'  , linewidth=1.8, label=cnoise)
                    plt.plot(vfl, 20.*nmp.log10(vnoise[ifr1:]+vrci95[ifr1:]), '0.4', linewidth=1.8, label='95% CI')
                    plt.plot(vfl, 20.*nmp.log10(vnoise[ifr1:]-vrci95[ifr1:]), '0.4', linewidth=1.8)
                plt.plot(vfl, 20.*nmp.log10(Vspec[ifr1:]),   '*-k', linewidth=2.)

        else:
            if not log_y:
                if l_do_ci95:
                    plt.plot(vfrq, vnoise,        '--k'  , linewidth=1.8)
                    plt.plot(vfrq, vnoise+vrci95, '0.4', linewidth=1.8)
                    plt.plot(vfrq, vnoise-vrci95, '0.4', linewidth=1.8)
                plt.plot(vfrq, Vspec,   '*-k', linewidth=2)
                if lplot_1onF: plt.plot(vfrq[1:], 0.03*1./vfrq[1:],   b_red, linewidth=2)
            else:
                if l_do_ci95:
                    plt.plot(vfrq, 20.*nmp.log10(vnoise),        '--k'  , linewidth=1.8)
                    plt.plot(vfrq, 20.*nmp.log10(vnoise+vrci95), '0.4', linewidth=1.8)
                    plt.plot(vfrq, 20.*nmp.log10(vnoise-vrci95), '0.4', linewidth=1.8)
                plt.plot(vfrq, 20.*nmp.log10(Vspec),   '*-k', linewidth=2)

        plt.ylabel('Amplitude Spectrum ('+cyunit+')', color='k', **font_xylb)

        plt.xlabel('Period (years)', color='k', **font_xylb)

        if log_x:
            x1=nmp.log10(x_max) ; x2=nmp.log10(x_min)
            plt.axis([x1, x2, rmin_amp, rmax_amp])
            if lplot_freq_ax:
                plt.xticks(nmp.log10(1./vee[:]),cvee[:])
            else:
                print ''
                vee_n = nmp.arange(vee[0], vee[len(vee)-1]-1, -1.)
                print vee_n[:]
                cvee_n = []
                for rr in vee_n:
                    cr = str(int(rr))
                    if cr in cvee:
                        cvee_n.append(cr)
                    else:
                        cvee_n.append('')
                print 'cvee =>', cvee[:]
                print 'cvee_n =>', cvee_n[:]
                plt.xticks(nmp.log10(1./vee_n[:]),cvee_n[:])
        else:
            x1=x_max; x2=x_min
            plt.axis([x1, x2, rmin_amp, rmax_amp])
            plt.xticks(1./vee[:],cvee[:], color='k')

        ax.grid(color='0.4', linestyle='-', linewidth=0.3)

        plt.legend(loc='upper left', shadow=False, fancybox=True)

        if lplot_freq_ax:
            ax2 = ax.twiny()
            ax2.set_xlabel('Frequency (cpy)', color='k', **font_xylb)
            if log_x:
                plt.axis([x1, x2, rmin_amp, rmax_amp])
                for jp in range(1,18,2): civee[jp] = ''
                plt.xticks(nmp.log10(1./vee[:]),civee)
                #ax2.xaxis.set_ticks(nmp.log10(1./vee[:]))
            else:
                plt.axis([x1, x2, rmin_amp, rmax_amp])
                plt.xticks(1./vee[:],civee)

            for t in ax2.get_xticklabels(): t.set_fontsize(14)
            ax2.xaxis.labelpad = 12 ; # move label upwards a bit...

        plt.savefig(cfignm+'.'+cfig_type, dpi=100, facecolor='w', edgecolor='w', orientation='portrait', transparent=False)
        plt.close(1)

        def __del__(self) :        
            plot.__counter -=  1


            #
    def __pow_spectrum_ssh(self, vk1, vps1, clab1=None, clr1=b_gre, lw1=6, \
                           cfig_name='fig_spectrum_SSH.png', cinfo='', logo_on=True, \
                           L_min=7., L_max=5000., P_min_y=-6, P_max_y=6,    \
                           l_show_k11o3=False, l_show_k5=False, l_show_k4=False, l_show_k2=False, \
                           vk2=[], vps2=[], clab2=None, clr2=b_org, lw2=3, \
                           vk3=[], vps3=[], clab3=None, clr3=b_blu, lw3=4, \
                           vk4=[], vps4=[], clab4=None, clr4='0.5', lw4=3 ):
        #------------------------------------------------------------------
        ## L_min=7. ; L_max : min and max wave-length for x-axis (km)
        #------------------------------------------------------------------
        #
        # r2Pi = 2.*nmp.pi # k is in rad/[space unit]
        r2Pi = 1. # k is in cycle/[space unit]
        #
        font_ttl, font_xylb, font_clb, font_inf = __font_unity__(fig_dpi=80)
        #
        # x-axis (lambda):
        k_min = r2Pi/L_max ; k_max = r2Pi/L_min
        xdef_l = nmp.asarray([ 4000., 2500., 1500., 1000., 700., 500., 300., 200., 150., 100., 70., 50., 40., 25., 15., 10., 7., 5., 4., 3., 2. ])
        (idx1,) = nmp.where(xdef_l>L_max) ; (idx2,) = nmp.where(xdef_l<L_min)
        xtcks_l = nmp.delete(xdef_l,nmp.concatenate((idx1,idx2)))
        cxtcks_l = []
        for rr in xtcks_l: cxtcks_l.append(str(int(rr)))
        xtcks_k  = r2Pi/xtcks_l

        fig = plt.figure(num = 1, figsize=(9.,9.), facecolor='w', edgecolor='k')
        ax = plt.axes([0.1, 0.07, 0.875, 0.86])

        if len(vk3) > 1 and len(vps3) > 1:
            plt.plot(nmp.log10(vk3), nmp.log10(vps3), '-', color=clr3, linewidth=lw3, label=clab3, zorder=10)
        if len(vk2) > 1 and len(vps2) > 1:
            plt.plot(nmp.log10(vk2), nmp.log10(vps2), '-', color=clr2, linewidth=lw2, label=clab2, zorder=15)
        if len(vk4) > 1 and len(vps4) > 1:
            plt.plot(nmp.log10(vk4), nmp.log10(vps4), '-', color=clr4, linewidth=lw4, label=clab4, zorder=4)

        plt.plot(    nmp.log10(vk1), nmp.log10(vps1), '-', color=clr1, linewidth=lw1, label=clab1, zorder=5)

        nl = len(vk1)
        #i1=4 ; i2=nl-0.4*nl ; rfct = 1.
        i1=int(0.53*float(nl)) ; i2=nl ; rfct = 3.
        if l_show_k2:
            plt.plot(nmp.log10(vk1[i1:i2]), nmp.log10((vk1[i1:i2]**-2.)/(rfct*1.E7)), '--', color='k', linewidth=2, label=r'k$^\mathregular{-2}$', zorder=2)
        if l_show_k4:
            plt.plot(nmp.log10(vk1[i1:i2]), nmp.log10((vk1[i1:i2]**-4.)/(rfct*3.E8)), '-', color='0.6', linewidth=2, label=r'k$^\mathregular{-4}$', zorder=1)
        if l_show_k5:
            plt.plot(nmp.log10(vk1[i1:i2]), nmp.log10((vk1[i1:i2]**-5.)/(rfct*2.E10)), '--', color='0.6', linewidth=2, label=r'k$^\mathregular{-5}$', zorder=2)
        if l_show_k11o3:
            plt.plot(nmp.log10(vk1[i1:i2]), nmp.log10((vk1[i1:i2]**(-11./3.))/(rfct*1.E9)), '-.', color='0.6', linewidth=2, label=r'k$^\mathregular{-11/3}$', zorder=2)
            
        # Bottom X-axis:
        plt.xticks( nmp.log10(xtcks_k), cxtcks_l)
        ax.set_xlim(nmp.log10(k_min), nmp.log10(k_max))
        ax.grid(color='k', linestyle='-', linewidth=0.2)
        plt.xlabel('Wave-length [km]')
        #
        # Y-axis:
        ax.set_ylim(P_min_y,P_max_y)
        cytcks = []
        for ii in range(P_min_y,P_max_y+1): cytcks.append(r'$\mathregular{10^{'+str(ii)+'}}$')
        plt.yticks( nmp.arange(P_min_y,P_max_y+1,1) , nmp.asarray(cytcks))
        plt.ylabel(r'SSH PSD [$\mathregular{m^2}$/(cy/km)]', color='k')
        #
        if clab1 != None: plt.legend(loc='best', shadow=True, fancybox=True) ; #lulu
        #
        # Top X-axis:
        ax2 = ax.twiny()
        P_max_x = 1 ; P_min_x = -4
        cxtcks_k = []
        for ii in range(P_min_x,P_max_x+1): cxtcks_k.append(r'$\mathregular{10^{'+str(ii)+'}}$')
        plt.xticks( nmp.arange(P_min_x,P_max_x+1,1) , nmp.asarray(cxtcks_k))
        ax2.set_xlim(nmp.log10(k_min), nmp.log10(k_max))
        #ax2.grid(color='0.3', linestyle='--', linewidth=0.2)
        [t.set_color('0.3') for t in ax2.xaxis.get_ticklabels()]
        plt.xlabel('Wave-number [cy/km]', color='0.3')
        #
        if cinfo != '': ax2.annotate(cinfo, xy=(0.08, 0.24), xycoords='axes fraction',  bbox={'facecolor':clr_inf_box, 'alpha':1., 'pad':10}, zorder=100, **font_inf)
        #
        if logo_on:
            fon = { 'fontname':'Arial', 'fontweight':'normal', 'fontsize':10 }
            ax2.annotate('© Ocean Next, 2018', xy=(0.84, -0.06), xycoords='axes fraction', color='0.5', zorder=100, **fon)
        #
        plt.savefig(cfig_name, dpi=120, facecolor='w', edgecolor='w', orientation='portrait')
        plt.close(1)
        return 0
    


    
    
    
# LOCAL functions
# ===============

def __get_mat__(cf):
    f1 = open(cf, 'r')
    lines1=f1.readlines()
    f1.close()
    zm   = [] ; jy   = 0
    for l in lines1:
        if l[0] != '#':
            jy = jy + 1
            ls = l.split()
            zm.append([])
            for c in ls:
                zm[jy-1].append(float(c))
    zxm = array(zm)
    print 'Shape zxm = ',nmp.shape(zxm), '\n'
    return zxm



def __vcontour__(zmin, zmax, zdc):
    if (zmin,zmax) == (0.,0.) or (zmin,zmax) == (0,0):
        vcont = [0.]
    else:
        lngt = zmax - zmin
        ncont = lngt/zdc
        vcont = nmp.arange(zmin, zmax + zdc, zdc)
    return vcont


def __name_longitude_ticks__(lon_min=0., lon_max=360., dlon=30., lon_ext=0):
    #
    # Builds nice ticks for X (lon) axis!
    #
    # Arrange longitude axis !
    VX = nmp.arange(lon_min, lon_max+lon_ext+dlon, dlon); VX0 = nmp.arange(lon_min, lon_max+lon_ext+dlon, dlon);
    ivf = nmp.where(VX>180); VX0[ivf] = VX[ivf] - 360
    cn_lon = []
    for rlon in VX0:
        jlon = int(rlon)
        if jlon < 0:
            cn_lon.append(str(-jlon)+r'$^{\circ}$W')
        else:
            if jlon == 0:
                cn_lon.append(str(jlon)+r'$^{\circ}$')
            else:
                cn_lon.append(str(jlon)+r'$^{\circ}$E')
    return VX, cn_lon

def __name_latitude_ticks__(lat_min=-90., lat_max=90., dlat=15.):
    #
    # Builds nice ticks for Y (lat) axis!
    #
    # Arrange latitude axis !
    VY = nmp.arange(lat_min, lat_max+dlat, dlat)
    cn_lat = []
    for rlat in VY:
        jlat = int(rlat)
        if jlat < 0:
            cn_lat.append(str(-jlat)+r'$^{\circ}$S')
        else:
            if jlat == 0:
                cn_lat.append(str(jlat)+r'$^{\circ}$')
            else:
                cn_lat.append(str(jlat)+r'$^{\circ}$N')
    return VY, cn_lat


def __name_coor_ticks__(lon_min=0., lon_max=360., dlon=30., lat_min=-90., lat_max=90., dlat=15., lon_ext=0):
    # Builds nice ticks for X and Y (lon, lat) axes!
    VX, cn_lon = __name_longitude_ticks__(lon_min=lon_min, lon_max=lon_max, dlon=dlon, lon_ext=lon_ext)
    VY, cn_lat =  __name_latitude_ticks__(lat_min=lat_min, lat_max=lat_max, dlat=dlat)
    return VX, VY, cn_lon, cn_lat


def __give_proj__(cname):

    nb =nmp.shape(projection_def)[0]

    vproj = [ 'NC', 'NC', 0.,  0.,  0.,  0.,  0.,  0., 'NC' ]
    jb = 0
    while jb < nb :
        if projection_def[jb][0] == cname:
            break
        else :
            jb = jb + 1
    if jb == nb :
        print 'Zone "'+cname+'" does not exist!\n'
        print 'so far choice is :'
        for jb in range(nb): print projection_def[jb][0]
        sys.exit(0)
    vproj = projection_def[jb][:]
    return vproj


def __font_unity__(fig_dpi=100., size='normal'):

    rat = 100./float(fig_dpi)

    if size == 'big': rat = 1.25*rat

    params = { 'font.family':'Trebuchet MS',
               'font.size':       int(13.*rat),
               'legend.fontsize': int(13.*rat),
               'xtick.labelsize': int(13.*rat),
               'ytick.labelsize': int(13.*rat),
               'axes.labelsize':  int(13.*rat),
               'legend.facecolor': 'white',
               'figure.facecolor': 'white' }

    mpl.rcParams.update(params)

    title_fonts    = { 'fontname':'Trebuchet MS'  , 'fontweight':'normal', 'fontsize':int(15.*rat) }
    label_fonts    = { 'fontname':'Trebuchet MS'  , 'fontweight':'normal', 'fontsize':int(14.*rat) }
    colorbar_fonts = { 'fontname':'Trebuchet MS'  , 'fontweight':'normal', 'fontsize':int(13.*rat) }
    info_fonts     = { 'fontname':'Ubuntu Mono'   , 'fontweight':'normal', 'fontsize':int(13.*rat) }

    return title_fonts, label_fonts, colorbar_fonts, info_fonts




def __force_min_and_max__(rm, rp, Xin):
    idx_bad  = nmp.where(nmp.logical_not(nmp.isfinite(Xin)))
    Xin[idx_bad] = 0.
    idx1 = nmp.where(Xin <= rm); Xin[idx1] = rm + abs(rp-rm)*1.E-4
    idx2 = nmp.where(Xin >= rp); Xin[idx2] = rp - abs(rp-rm)*1.E-4
    Xin[idx_bad] = nmp.nan


def __subsample_colorbar__(i_sbsmp, vcc, clb_hndl, cb_or='vertical'):
    cb_labs = []

    # First checking if vcc countains integers or not...
    lcint = False
    vc = vcc.astype(nmp.int64)  ; # integer version of vcc
    if nmp.max(nmp.abs(vcc))>5 and nmp.sum(vcc-vc) == 0. : lcint=True

    cpt = 0
    nn = int(round(abs(vcc[-1]-vcc[0])/abs(vcc[0]-vcc[1]),0))
    if nn % 2 != 0: cpt = 1

    for rr in vcc:
        if cpt % i_sbsmp == 0:
            if lcint:
                cr = str(int(rr))
            else:
                cr = str(round(float(rr),6))
            cb_labs.append(cr)
        else:
            cb_labs.append(' ')
        cpt = cpt + 1
    if cb_or == 'horizontal':
        clb_hndl.ax.set_xticklabels(cb_labs)
    else:
        clb_hndl.ax.set_yticklabels(cb_labs)
    del cb_labs, vc



def __nice_colorbar__(fig_hndl, plt_hndl, vcc,
                      cax_other=None, i_sbsmp=1, lkc=False, cb_or='vertical', cunit=None, cfont=None, fontsize=0):

    if cb_or not in {'horizontal','vertical'}:
        print "ERROR: only 'vertical' or 'horizontal' can be specified for the colorbar orientation!"
        cb_or = 'vertical'

    if cb_or == 'horizontal':
        if cax_other is not None:
            clb = plt_hndl.colorbar(fig_hndl, cax=cax_other, ticks=vcc, drawedges=lkc, orientation='horizontal',
                                    pad=0.07, shrink=1., aspect=40, extend='both')
        else:
            clb = plt_hndl.colorbar(fig_hndl,                ticks=vcc, drawedges=lkc, orientation='horizontal',
                                    pad=0.07, shrink=1., aspect=40, extend='both')
    else:
        if cax_other is not None:
            clb = plt_hndl.colorbar(fig_hndl, cax=cax_other, ticks=vcc, drawedges=lkc,
                                    pad=0.03, extend='both')
        else:
            clb = plt_hndl.colorbar(fig_hndl,                ticks=vcc, drawedges=lkc,
                                    pad=0.03, extend='both')

    if i_sbsmp > 1: __subsample_colorbar__(i_sbsmp, vcc, clb, cb_or=cb_or)

    if not cunit is None:
        if cfont is None:
            clb.set_label(cunit)
        else:
            clb.set_label(cunit, **cfont)

    if fontsize > 0:
        if cb_or == 'horizontal':
            for t in clb.ax.get_xticklabels(): t.set_fontsize(fontsize) # Font size for colorbar ticks!
        else:
            for t in clb.ax.get_yticklabels(): t.set_fontsize(fontsize) # Font size for colorbar ticks!



def _add_xy_offset__(plt_hndl, ixo, iyo):
    if ( ixo != 0. ):
        locs, labels = plt_hndl.xticks() ; jl=0
        vlabs = []
        for ll in locs:
            clab = str(int(locs[jl])+int(ixo))
            vlabs.append(clab); jl=jl+1
        plt_hndl.xticks(locs,vlabs)
    if ( y_offset != 0. ):
        locs, labels = plt_hndl.yticks() ; jl=0; vlabs = []
        for ll in locs:
            clab = str(int(locs[jl])+int(iyo))
            vlabs.append(clab); jl=jl+1
        plt_hndl.yticks(locs,vlabs)
    del vlabs

def __subsample_axis__(plt_hndl, cax, i_sbsmp, icpt=1):
    ax_lab = []
    if   cax == 'x':
        locs, labels = plt_hndl.xticks()
    elif cax == 'y':
        locs, labels = plt_hndl.yticks()
    else:
        print ' Error: __subsample_axis__.clprn_plot => only "x" or "y" please'; sys.exit(0)
    cpt = icpt # with ipct = 1: tick priting will start at y1+dt on x axis rather than y1
    for rr in locs:
        if cpt % i_sbsmp == 0:
            if rr%1.0 == 0.:
                cr = str(int(rr))  # it's something like 22.0, we want 22 !!!
            else:
                cr = str(rr)
            ax_lab.append(cr)
        else:
            ax_lab.append(' ')
        cpt = cpt + 1
    if cax == 'x': plt_hndl.xticks(locs,ax_lab)
    if cax == 'y': plt_hndl.yticks(locs,ax_lab)
    del ax_lab



def __nice_x_axis__(ax_hndl, plt_hndl, x_0, x_H, dx, i_sbsmp=1, cunit=None, cfont=None, dx_minor=0):
    x_l = x_0
    if x_0%dx != 0.: x_l = float(int(x_0/dx))*dx
    if x_H%dx != 0.: x_H = float(int(x_H/dx)+1)*dx
    plt_hndl.xticks( nmp.arange(x_l, x_H+dx, dx) )
    locs, labels = plt_hndl.xticks()
    ax_hndl.get_xaxis().get_major_formatter().set_useOffset(False) ; # Prevents from using scientific notations in axess ticks numbering...
    if i_sbsmp > 1: __subsample_axis__( plt, 'x', i_sbsmp)
    if not cunit is None:
        if cfont is None:
            plt_hndl.xlabel(cunit)
        else:
            plt_hndl.xlabel(cunit, **cfont)
    # Add minor x-ticks and corresponding grid:
    if dx_minor > 0:
        locs, labels = plt_hndl.xticks()
        ax_hndl.set_xticks( nmp.arange(locs[0], locs[len(locs)-1] , dx_minor) , minor=True)
        ax_hndl.grid(which='both')
        ax_hndl.grid(which='minor', color='k', linestyle='-', linewidth=0.1)
    ax_hndl.grid(which='major', color='k', linestyle='-', linewidth=0.2)
    ax_hndl.set_xlim(x_l,x_H+dx/1000.)

def __nice_y_axis__(ax_hndl, plt_hndl, y_0, y_H, dy, i_sbsmp=1, cunit=None, cfont=None, dy_minor=0):
    y_l = y_0
    if y_0%dy != 0.: y_l = float(int(y_0/dy))*dy
    plt_hndl.yticks( nmp.arange(y_l, y_H+dy, dy) )
    locs, labels = plt_hndl.yticks()
    ax_hndl.get_yaxis().get_major_formatter().set_useOffset(False) ; # Prevents from using scientific notations in axess ticks numbering...
    if i_sbsmp > 1: __subsample_axis__( plt, 'y', i_sbsmp)
    if not cunit is None:
        if cfont is None:
            plt_hndl.ylabel(cunit)
        else:
            plt_hndl.ylabel(cunit, **cfont)
    # Add minor y-ticks and corresponding grid:
    if dy_minor > 0:
        locs, labels = plt_hndl.yticks()
        ax_hndl.set_yticks( nmp.arange(locs[0], locs[len(locs)-1] , dy_minor) , minor=True)
        ax_hndl.grid(which='both')
        ax_hndl.grid(which='minor', color='k', linestyle='-', linewidth=0.1)
    ax_hndl.grid(which='major', color='k', linestyle='-', linewidth=0.2)
    ax_hndl.set_ylim(y_0,y_H+dy/1000.)

def __nice_depth_axis__(ax_hndl, plt_hndl, z0, zK, l_log=False, l_z_inc=True, cunit=None, cfont=None):
    ax_hndl.get_yaxis().get_major_formatter().set_useOffset(False)
    if l_log:
        y_log_ofs = 10.
        vyview_list = [ 3. , 10. , 25., 50. , 100. , 250. , 500. , 1000. , 2500.,  5000. ]
        nd = len(vyview_list)
        vyview      = nmp.zeros(nd)
        for jn in range(nd): vyview[jn] = vyview_list[jn]
        vyview_log = nmp.log10(vyview + y_log_ofs)
        ylab = []
        for rr in vyview_list: ylab.append(str(int(rr)))
        z0 = nmp.log10(z0+y_log_ofs)
        zK = nmp.log10(zK+y_log_ofs)
        ax_hndl.set_yticks(vyview_log)
        ax_hndl.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
        ax_hndl.set_yticklabels(ylab)
    if l_z_inc:
        ax_hndl.set_ylim(z0,zK)
    else:
        ax_hndl.set_ylim(zK+(zK-z0)/50. , z0)
    ax_hndl.grid(color='k', linestyle='-', linewidth=0.5)
    if not cunit is None:
        if cfont is None: cfont  =  { 'fontname':'Ubuntu Mono', 'fontweight':'normal', 'fontsize':10, 'color':'k'}
        plt_hndl.ylabel(cunit, **cfont)

def __nice_latitude_axis__(ax_hndl, plt_hndl, lmin, lmax, dl, axt='x'):
    if axt == 'x':
        ax_hndl.get_xaxis().get_major_formatter().set_useOffset(False)
    elif axt == 'y':
        ax_hndl.get_yaxis().get_major_formatter().set_useOffset(False)
    else:
        print 'ERROR clprn_plot.__nice_latitude_axis__: only accepts "x" or "y" for axt!'
        sys.exit(0)
    [vvl, ctck] =  __name_latitude_ticks__(lat_min=lmin, lat_max=lmax, dlat=dl)
    if axt == 'x':
        plt_hndl.xticks(vvl,ctck)
        ax_hndl.set_xlim(lmin,lmax)
    else:
        plt_hndl.yticks(vvl,ctck)
        ax_hndl.set_ylim(lmin,lmax)


def __nice_longitude_axis__(ax_hndl, plt_hndl, lmin, lmax, dl, axt='x'):
    if axt == 'x':
        ax_hndl.get_xaxis().get_major_formatter().set_useOffset(False)
    elif axt == 'y':
        ax_hndl.get_yaxis().get_major_formatter().set_useOffset(False)
    else:
        print 'ERROR clprn_plot.__nice_longitude_axis__: only accepts "x" or "y" for axt!'
        sys.exit(0)
    [vvl, ctck] =  __name_longitude_ticks__(lon_min=lmin, lon_max=lmax, dlon=dl)
    if axt == 'x':
        plt_hndl.xticks(vvl,ctck)
        ax_hndl.set_xlim(lmin,lmax)
    else:
        plt_hndl.yticks(vvl,ctck)
        ax_hndl.set_ylim(lmin,lmax)



def __prepare_z_log_axis__(l_log, vz):
    import math
    nk  = len(vz)
    zvz = nmp.zeros(nk)
    if l_log:
        for jk in range(nk):
            zvz[jk] = math.log10(vz[jk])
    else:
        zvz= vz
    return zvz


def __nb_col_row_legend__(nn):
    if nn <= 3:
        nbc = 1 ; nbr = nn
    elif nn == 4:
        nbc = 2 ; nbr = 2
    elif nn > 4 and nn <= 6:
        nbc = 2 ; nbrfull = 2; nfull = nbc*nbrfull; nbr = nbrfull + nn/nfull
    elif nn > 6 and nn <= 9:
        nbc = 3 ; nbrfull = 2; nfull = nbc*nbrfull; nbr = nbrfull + nn/nfull
    elif nn > 9 and nn <= 16:
        nbc = 4 ; nbrfull = 3; nfull = nbc*nbrfull; nbr = nbrfull + nn/nfull
    else:
        nbc = 4 ; nbr = nn/nbc + 1
    return nbc, nbr

def __time_axis_minor_ticks__(dt):
    dt_mnr=0
    if ((dt>=2)  and (dt<10) and (dt%2 == 0)) or (dt==5) : dt_mnr=1
    if (dt>=10) and (dt<50) and (dt%5 == 0) : dt_mnr=5
    if (dt>=50) and (dt%50 == 0) : dt_mnr=10
    return dt_mnr


def __suitable_axis_dx__(hmin, hmax, nb_val=20, lsym0=False):
    if (hmin,hmax) == (0.,0.) or (hmin,hmax) == (0,0):
        dh = 0.
    else:
        dh = abs(hmax - hmin)/float(nb_val)
        lfound = False
        iexp = 20
        while not lfound :
            if dh%(10.**(iexp-1)) != dh: lfound = True
            iexp = iexp - 1
        if iexp < 1:
            dh = round(dh, -iexp)
        else:
            dh = round(dh,0)
            dh = round(dh/(10.**iexp),0)*10.**iexp
    
        dhi = dh*10.**(-iexp)
        if dhi == 3.:         dh = 2.5*10.**(iexp)
        if dhi in [4.,6.,7.]: dh =  5.*10.**(iexp)
        if dhi in [8.,9.]:    dh = 10.*10.**(iexp) ; iexp=iexp+1
    
        hmin = float(int(hmin*10.**(-iexp)))*10.**iexp
        hmax = float(int((hmax+dh)*10.**(-iexp)))*10.**iexp
    
        if lsym0:
            # Force symetry about 0 !
            hmax = max(abs(hmax),abs(hmin))
            if hmax%dh != 0.: hmax = float(int(hmax/dh))*dh
            hmin = -hmax

    return hmin, hmax, dh



def __fancy_legend__(ax_hndl, plt_hndl, loc_leg='0', ylg=0, leg_out=False, ncol=1):
    if loc_leg != '0':
        if leg_out:
            # Shrink Y axis's height by % on the bottom
            box = ax_hndl.get_position()
            ax_hndl.set_position([box.x0, box.y0 + box.height*ylg, box.width, box.height*(1.-ylg)])
            plt_hndl.legend(bbox_to_anchor=(0.95, -0.075), ncol=ncol, shadow=True, fancybox=True)
        else:
            plt_hndl.legend(loc=loc_leg, ncol=ncol, shadow=True, fancybox=True)

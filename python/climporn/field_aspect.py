from sys import exit

class field_aspect:

    ''' This is where you can define field specific plotting features such as colormap, bounds, etc...
        For now, only used by `nemo_movie_sclr.py`
    '''

    def __init__( self, cfield, cbox="" ):
        self.cfield = cfield
        self.cbox   = cbox

        CWHAT = self.cfield
        CBOX  = self.cbox

        # Defaults:
        color_top    = 'white'
        color_top_cb = 'white'

        cv_in  = CWHAT
        cv_out = CWHAT
        cv_msk = 'tmask'

        l_show_lsm = True
        l_show_ice = False ; #??? keep?

        l_log_field = False
        l_pow_field = False
        pow_field=0.5
        vc_fld_powlog = [ ]


        rmult = 1.  ; # multiplicator to apply to field read into netCDF file
        tmin = 0.
        tmax = 1.
        df   = 0.1
        cpal_fld = 'viridis'
        cb_jump = 1
        cb_extend = 'both'        
        cunit = ''

        l_apply_geov  = False
        l_apply_lap   = False
        l_apply_hgrad = False
        l_smooth  = False
        nb_smooth = 0


        if   CWHAT == 'MLD':
            cv_in = 'somxl010' ; cv_out = 'MLD'
            tmin=0. ;  tmax=1800.  ;  df = 50. ; cpal_fld='ncview_hotres' ;     cb_jump = 4
            cunit = r'MLD [m]'
            if CBOX == 'Nordic': tmin=0. ; tmax=1000. ; cb_jump = 2 ; cpal_fld='magma_r' ; color_top_cb='k'

        elif CWHAT == 'SST':
            cv_in = 'sosstsst' ; cv_out = CWHAT ; #in ['sosstsst','tos']:
            tmin=-2 ;  tmax=32.   ;  df = 1. ; cpal_fld='ncview_nrl' ;     cb_jump = 2
            cunit = r'SST ($^{\circ}$C)'
            if CBOX == 'EUROPA':   tmin=0. ;  tmax=25.
            if CBOX == 'EUROPAs':  tmin=6. ;  tmax=18.
            if CBOX == 'Med':      tmin=10.;  tmax=30. ; cb_jump = 1
            if CBOX == 'Med+BS':   tmin=7. ;  tmax=25.
            if CBOX == 'LabSea':   tmin=-2.;  tmax=15.
            if CBOX == 'Brittany': tmin=7. ;  tmax=13.
            if CBOX == 'GrlIcl':   tmax = 12.
            if CBOX in [ 'AzoresP','AzoresL','AzoresS']:  tmin = 15. ; tmax = 25. ; df=0.5
            if CBOX == 'Bretagne': tmin = 10. ; tmax = 22. ; df=1.
            if CNEMO == "CALEDO60": l_show_ice = False ; tmin=18. ;  tmax=30. ; cb_jump=1 ; df=1 ; cv_in = 'tos'

        elif CWHAT == 'T_1000':
            cv_in = 'votemper' ; cv_out = CWHAT ;
            tmin=0. ;  tmax=14.   ;  df = 1. ; cpal_fld='ncview_nrl' ; cb_jump = 1
            cunit = r'Potential temperature at 1000 m'

        elif CWHAT == 'T_60':
            cv_in = 'votemper' ; cv_out = CWHAT ;
            tmin=0. ;  tmax=14.   ;  df = 1. ; cpal_fld='ncview_nrl' ; cb_jump = 1
            cunit = r'Potential temperature at 60 m'
            if CBOX == 'BlackSea' : tmin=0. ; tmax=20. ;  df = 1.   ; cb_jump = 1 ; #cpal_fld='gist_stern_r'

        elif CWHAT == 'SSS':
            cv_in = 'sosaline' ; cv_out = CWHAT ; #in ['sosstsst','tos']:
            tmin=20. ;  tmax=40.   ;  df = 2. ; cpal_fld='ncview_ssec' ; cb_jump = 2
            cunit = r'Sea surface salinity'
            if CBOX == 'Med' :    tmin=33. ; tmax=39.5 ;  df = 0.25 ; cpal_fld='magma'
            if CBOX == 'LabSea' : tmin=28. ; tmax=35.5 ;  df = 1.   ; cb_jump = 1 ; cpal_fld='gist_stern_r'

        elif CWHAT == 'S_1000':
            cv_in = 'vosaline' ; cv_out = CWHAT ;
            tmin=33.5 ;  tmax=36.5   ;  df = 0.5 ; cpal_fld='ncview_helix2' ; cb_jump = 1
            cunit = r'Salinity at 1000 m'
            if CBOX == 'MeddiesW' : tmin=35. ;  tmax=36.6 ; df = 0.1

        elif CWHAT == 'GRAD_SST':
            cv_in = 'sosstsst' ; cv_out = CWHAT ;
            l_apply_hgrad = True
            l_smooth = True ; nb_smooth  = 5
            tmin=0. ;  tmax=0.001 ;  df = 0.0001 ; cpal_fld='ncview_hotres' ; cb_jump = 1
            cunit = r'$\left|\vec{\nabla}SST\right|$ (K/m)'

        elif CWHAT == 'SSU':
            cv_in = 'sozocrtx' ; cv_out = CWHAT ;
            tmin=-0.8 ;  tmax=-tmin  ;  df = 0.1 ; cpal_fld='bone' ; cb_jump = 1
            cunit = r'Zonal component of current speed [m/s]'
            cv_msk = 'umask'

        elif CWHAT == 'SSV':
            cv_in = 'somecrty' ; cv_out = CWHAT ;
            tmin=-1. ;  tmax=-tmin  ;  df = 0.2 ; cpal_fld='bone' ; cb_jump = 1
            cunit = r'Meridional component of current speed [m/s]'
            cv_msk = 'vmask'

        elif CWHAT == 'SSH':
            cv_in = 'sossheig' ; cv_out = CWHAT ;
            cpal_fld='RdBu_r' ; tmin=-3. ;  tmax=-tmin   ;  df = 0.5 ;
            cb_jump = 1
            cunit = r'SSH [m]'
            if CBOX == 'Med' or CBOX == 'Med+BS': tmin=-0.7; tmax=0.2   ; df = 0.1
            if CRUN[:4] == 'BLB0':                tmin=-1.2; tmax=-tmin ; df = 0.2
            if CBOX in [ 'Bretagne']:             tmin=-4.;  tmax=-tmin ; df = 0.5

        elif CWHAT == 'CURLOF':
            cv_in = 'socurloverf' ; cv_out = CWHAT ;
            #tmin=-0.8 ;  tmax=-tmin  ;  df = 0.1  ; cb_jump = 2 ;
            tmin=-0.7 ;  tmax=-tmin  ;  df = 0.1  ; cb_jump = 1
            cpal_fld='RdBu_r' ; color_top_cb='k' ; # cpal_fld='on2'
            cunit = r'$\zeta/f$'
            cv_msk = 'vmask'

        elif CWHAT == 'GEOSSV':
            # Geostrophic velocity speed out of SSH
            cv_in = 'sossheig' ; cv_out = CWHAT ;
            l_apply_geov = True
            cpal_fld='on3' ; tmin=0. ;  tmax=1.2   ;  df = 0.2 ; cb_extend = 'max'
            cb_jump = 1
            cunit = r'Surface geostrophic velocity speed [m/s]'
            l_save_nc = True

        elif CWHAT == 'LAP_SSH':
            cv_in = 'sossheig' ; cv_out = CWHAT ;
            l_apply_lap = True
            cpal_fld='on2' ; tmin=-1.2 ;  tmax=1.2   ;  df = 0.05 ;

        elif CWHAT == 'W_1000':
            cv_in = 'vovecrtz'  ; cv_out = CWHAT ;
            tmin=-0.01 ;  tmax=-tmin   ;  df = 0.005 ; cb_jump = 1
            cpal_fld='RdBu_r' ;    #cpal_fld='PiYG_r' ; #cpal_fld='BrBG_r'
            cunit = r'Vertical velocity at 1000 m [m/s]'
            if CBOX in [ 'AzoresP','AzoresL','AzoresS']: color_top = 'k'

        elif CWHAT == 'Amplitude':
            cv_in = 'r'     ; cv_out = cv_in ;
            cpal_fld='RdBu_r' ; tmin=-0.5 ;  tmax=-tmin   ;  df = 0.1 ; cb_jump = 1
            cunit = r'Amplitude [m]'

        elif CWHAT == 'Phase':
            cv_in = 'phi'     ; cv_out = cv_in ;
            cpal_fld='RdBu_r' ; tmin=-30. ;  tmax=-tmin   ;  df = 5. ; cb_jump = 1
            #
            cunit = r'Phase (deg.)'

        elif CWHAT == 'siconc':
            cv_in = 'siconc'  ; cv_out = cv_in ; color_top_cb='k'
            cpal_fld='ncview_ice' ; tmin=0. ;  tmax=1. ;  df = 0.1 ; cb_jump = 1 ; cb_extend = 'neither'
            cunit = 'Sea-ice concentration'

        elif CWHAT == 'sivolu':
            cv_in = 'sivolu'  ; cv_out = cv_in ; color_top_cb='k'
            cpal_fld='magma' ; tmin=0. ;  tmax=4. ;  df=1 ; cb_jump = 1; cb_extend = 'max'
            cunit = 'Sea-ice volume [m]'

        elif CWHAT == 'sithic':
            cv_in = 'sithic'  ; cv_out = cv_in ; color_top_cb='k'
            cpal_fld='magma' ; tmin=0. ;  tmax=4. ;  df=1 ; cb_jump = 1; cb_extend = 'max'
            cunit = 'Sea-ice thickness [m]'

        elif CWHAT in [ 'damage', 'damage-t', 'damage-f' ]:
            cv_in = CWHAT  ; cv_out = cv_in
            cpal_fld='magma' ; color_top_cb='k'
            tmin=0. ;  tmax=1. ; l_pow_field=True ; pow_field=7. ; cb_extend = 'neither'
            vc_fld_powlog = [ 0., 0.7, 0.8, 0.9, 0.95, 1. ]
            cunit = 'Damage@T'

        elif CWHAT in [ 'sidive', 'sidive-t', 'sidive-f' ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='k'
            rmult = 3600.*24.
            cpal_fld='RdBu_r' ; tmin=-0.2 ;  tmax=-tmin ;  df=0.05
            cunit = 'Divergence of sea-ice velocity [day$^{-1}$]'

        elif CWHAT in [ 'zfUu'  ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='k'
            #rmult = 3600.*24.
            cpal_fld='RdBu_r' ; tmin=-0.2 ;  tmax=-tmin ;  df=0.05
            cunit = r'$\partial ( h_t ~ \sigma_{11} ) / \partial x + \partial ( h_f ~ \sigma_{12} ) / \partial y$ [Pa]'

        elif CWHAT in [ 'sigI', 'sig_I', 'ice_sigI', 'ice_sig_I', 'ice_sigI-t', 'ice_sigI-f'  ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='k'
            cpal_fld='turbo_r' ; tmin=-50000. ;  tmax=10000. ;  df=10000.
            cunit = r'$\sigma_{I}$ [Pa]'

        elif CWHAT in [ 'sishea', 'sishea-t', 'sishea-f' ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='w'
            rmult = 3600.*24. ; color_top_cb='k'
            cpal_fld='inferno' ; tmin=0. ;  tmax=5.
            l_pow_field=True ; pow_field=0.25 ;  vc_fld_powlog= [tmin, 0.1, 0.5, 1., 3., tmax ]
            cunit = 'Shear of sea-ice velocity [day$^{-1}$]'

        else:
            print('ERROR: we do not know field '+str(CWHAT)+' !')
            exit(0)

            
        self.color_top     = color_top
        self.color_top_cb  = color_top_cb

        self.cv_in         = cv_in
        self.cv_out        = cv_out
        self.cv_msk        = cv_msk

        self.l_show_lsm    = l_show_lsm
        self.l_show_ice    = l_show_ice

        self.l_log_field   = l_log_field
        self.l_pow_field   = l_pow_field
        self.pow_field     = pow_field
        self.vc_fld_powlog = vc_fld_powlog

        self.rmult         = rmult
        self.tmin          = tmin
        self.tmax          = tmax
        self.df            = df
        self.cb_jump       = cb_jump
        self.cb_extend     = cb_extend
        self.cpal_fld      = cpal_fld
        self.cunit         = cunit

        self.l_apply_geov  = l_apply_geov
        self.l_apply_lap   = l_apply_lap
        self.l_apply_hgrad = l_apply_hgrad
        self.l_smooth      = l_smooth
        self.nb_smooth     = nb_smooth


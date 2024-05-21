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

        vc_fld_force = []
        
        l_log_field = False
        l_pow_field = False
        pow_field=0.5
        vc_fld_powlog = [ ]

        imask_no_ice_pc = 0 ; # integer, ice concntration in % below which we hide the field...
        nm_ice_conc = 'siconc'   ; # => in this case we need to find the ice concentration bearing this name !

        color_missing = '0.8'
        
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
        lSubstractMean = False
        l_smooth  = False
        nb_smooth = 0


        if   CWHAT == 'MLD':
            cv_in = 'somxl010' ; cv_out = 'MLD'
            tmin=0. ;  tmax=1800.  ;  df = 50. ; cpal_fld='magma_r' ;     cb_jump = 4
            cunit = r'MLD [m]'
            if CBOX in ['Nordic']:
                tmin=0. ; tmax=1000. ; cb_jump = 2 ; cpal_fld='magma_r' ; color_top_cb='k'
            if CBOX in ['ALL','ALLFR','MNATFR']:
                tmin=0. ; tmax=1000. ; cb_jump = 2 ; cpal_fld='magma_r' ; color_top_cb='w'

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
            if CBOX == 'BlackSea' : tmin=0. ; tmax=20. ;  df = 1.   ; #cpal_fld='gist_stern_r'

        elif CWHAT == 'SSS':
            cv_in = 'sosaline' ; cv_out = CWHAT ; #in ['sosstsst','tos']:
            tmin=20. ;  tmax=40.   ;  df = 2. ; cpal_fld='ncview_ssec' ; cb_jump = 2
            cunit = r'Sea surface salinity'
            if CBOX == 'Med' :    tmin=33. ; tmax=39.5 ;  df = 0.25 ; cpal_fld='magma'
            if CBOX == 'LabSea' : tmin=28. ; tmax=35.5 ;  df = 1.   ; cpal_fld='gist_stern_r'

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
            if CBOX in [ 'Bretagne']:             tmin=-4.;  tmax=-tmin ; df = 0.5
            if CBOX == 'Bahamas':
                lSubstractMean = True
                tmin=-0.4;  tmax=-tmin ; df = 0.1
                color_top    = 'k'
                color_top_cb = 'k'



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
            cv_in = 'siconc'  ; cv_out = cv_in
            #l_pow_field=True ; pow_field=2. ; vc_fld_powlog = [ 0., 0.5, 0.7, 0.8, 0.9, 0.95, 1. ]
            cpal_fld='ncview_ice' ; tmin=0. ;  tmax=1. ;  df = 0.1 ; cb_extend = 'neither' ; color_top_cb='k'; color_top='k'
            #cpal_fld='gray' ; tmin=0. ;  tmax=1. ;  df = 0.1 ; cb_extend = 'neither' ; color_top_cb='k'; color_top='k'
            l_pow_field=True ; pow_field=6. ; vc_fld_powlog = [ 0., 0.7, 0.8, 0.9, 0.95, 1. ]
            cunit = 'Sea-ice concentration'

        elif CWHAT in [ 'sithic', 'sivolu' ]:
            cv_in = CWHAT  ; cv_out = cv_in
            cpal_fld='cmocean_ice' ; color_top_cb='w' ; color_top = 'w'
            imask_no_ice_pc = 5 ; #color_missing = 'k' ; # we hide the field where A<`imask_no_ice_pc`%
            tmin=0. ;  tmax=5. ;  df=1 ; cb_jump = 1; cb_extend = 'max'
            cunit = 'Sea-ice thickness [m]'
            if CBOX in ["CentralArctic"]: imask_no_ice_pc=0.1 ; color_missing='k'
            if CBOX in ["Baffin"]:        tmin=0; tmax=2;  df=0.5; imask_no_ice_pc=0.1 ; color_missing='k'
            if CBOX in ["Spitzberg"]:     tmin=0; tmax=4;  df=1;   imask_no_ice_pc=0.1 ; color_top_cb='k'
            
        elif CWHAT in [ 'damage', 'damage-t', 'damage-f' ]:
            cv_in = CWHAT  ; cv_out = cv_in            
            #
            imask_no_ice_pc = 5 ; color_missing = '#2D4B87' ; # we hide the field where A<`imask_no_ice_pc`%
            #
            #cpal_fld='bone_r' ; color_top_cb='w'; tmin=0. ;  tmax=1. ; l_pow_field=True ; pow_field=5.; vc_fld_powlog=[ 0., 0.7, 0.8, 0.9, 0.95, 1. ]
            #cpal_fld='gray_r' ; color_top_cb='w'; tmin=0. ;  tmax=1. ; l_pow_field=True ; pow_field=7.; vc_fld_powlog=[ 0., 0.8, 0.9, 0.95, 1. ]; color_missing = None
            #cpal_fld='ncview_ice' ; color_top_cb='w'; tmin=0. ;  tmax=1. ; l_pow_field=True ; pow_field=5.; vc_fld_powlog=[ 0., 0.7, 0.8, 0.9, 0.95, 1. ]
            #
            cpal_fld='bone_r' ; color_top_cb='w'; color_top='w'; tmin=0. ;  tmax=1. ; l_pow_field=True ; pow_field=5.; vc_fld_powlog=[ 0., 0.7, 0.8, 0.9, 0.95, 1. ]
            #cpal_fld='gray_r' ; color_top_cb='k'; color_top='k'; tmin=0.78 ;  tmax=1. ; l_pow_field=False ; pow_field=10.; vc_fld_force=[  0.8, 0.85, 0.9, 0.95, 1. ]
            #
            # Paper:
            #imask_no_ice_pc = 50 ; color_missing = 'w' ; # we hide the field where A<`imask_no_ice_pc`%
            #cpal_fld='gray' ; color_top_cb='k'; color_top='k'; tmin=0.5 ;  tmax=1. ; l_pow_field=True ; pow_field=5.; vc_fld_powlog=[ 0.5, 0.85, 0.9, 0.95, 0.98, 1. ]
            #
            cb_extend = 'neither'
            #
            cunit = 'Sea-ice damage'
            #
            if CBOX=="CentralArctic": cpal_fld='cmocean_gray_r'; pow_field=3.; color_top='k'; color_top_cb='k'


            
        elif CWHAT in [ 'sivelo', 'sivelo-t', 'sivelo-f' ]:
            cv_in = CWHAT  ; cv_out = CWHAT
            #cpal_fld='magma' ; color_top_cb='k' 
            #cpal_fld='gnuplot2' ; color_top_cb='k'
            #cpal_fld='cmocean_tempo' ; color_top_cb='w'
            cpal_fld='cmocean_dense' ; color_top_cb='k'            
            tmin=0. ; tmax=0.6 ; df=0.1 ; cb_extend = 'max'
            imask_no_ice_pc = 5 ; # we hide the field where A<10%
            cunit = 'Sea-ice velocity [m/s]'

        elif CWHAT in [ 'sidive', 'sidive-t', 'sidive-f' ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='k'
            rmult = 3600.*24.
            cpal_fld='RdBu_r' ; tmin=-0.2 ;  tmax=-tmin ;  df=0.05
            cunit = 'Divergence of sea-ice velocity [day$^{-1}$]'

        elif CWHAT in [ 'sishear', 'sishear-t', 'sishear-f' ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='k'
            rmult = 3600.*24.
            cpal_fld='RdBu_r' ; tmin=-0.2 ;  tmax=-tmin ;  df=0.05
            imask_no_ice_pc = 5 ; #color_missing = 'k' ; # we hide the field where A<`imask_no_ice_pc`%
            cunit = 'Shear of sea-ice velocity [day$^{-1}$]'

        elif CWHAT in [ 'sivort', 'sivort-t', 'sivort-f' ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='k'
            rmult = 3600.*24.
            cpal_fld='RdBu_r' ; tmin=-0.3 ;  tmax=-tmin ;  df=0.1
            #cpal_fld='cmocean_balance_r' ; tmin=-0.3 ;  tmax=-tmin ;  df=0.1
            #cpal_fld='BrBG_r' ; tmin=-0.3 ;  tmax=-tmin ;  df=0.1
            imask_no_ice_pc = 5 ; #color_missing = 'k' ; # we hide the field where A<`imask_no_ice_pc`%
            cunit = 'Vorticity of sea-ice velocity [day$^{-1}$]'

        elif CWHAT in [ 'zfUu'  ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='k'
            #rmult = 3600.*24.
            cpal_fld='RdBu_r' ; tmin=-0.2 ;  tmax=-tmin ;  df=0.05
            cunit = r'$\partial ( h_t ~ \sigma_{11} ) / \partial x + \partial ( h_f ~ \sigma_{12} ) / \partial y$ [Pa]'

        elif CWHAT in [ 'sigI', 'sig_I', 'ice_sigI', 'ice_sig_I', 'ice_sigI-t', 'ice_sigI-f'  ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='k'
            cpal_fld='turbo_r' ; tmin=-50000. ;  tmax=10000. ;  df=10000.
            cunit = r'$\sigma_{I}$ [Pa]'

        elif CWHAT in [ 'sishea', 'sishea-t', 'sishea-f', 'simxshr-t', 'simxshr-f' ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='w'
            rmult = 3600.*24.
            cpal_fld='inferno' ; tmin=0. ;  tmax=1. ; cb_extend = 'max'
            imask_no_ice_pc = 5 ; #color_missing = 'k' ; # we hide the field where A<`imask_no_ice_pc`%
            l_pow_field=True ; pow_field=0.3 ;  vc_fld_powlog= [tmin, 0.01, 0.1, 0.5, tmax ]
            cunit = 'Maximum shear of sea-ice velocity [day$^{-1}$]'

        elif CWHAT in [ 'sidefo', 'sidefo-t', 'sidefo-f' ]:
            cv_in = CWHAT  ; cv_out = cv_in ; color_top_cb='w'
            rmult = 3600.*24.
            #cpal_fld='inferno' ; tmin=0. ;  tmax=1.5 ; cb_extend = 'max'
            cpal_fld='gray_r' ; tmin=0. ;  tmax=1.5 ; cb_extend = 'max' ; color_top_cb='k'; color_top='k'
            #cpal_fld='viridis' ; tmin=0. ;  tmax=1. ; cb_extend = 'max'
            #imask_no_ice_pc = 50 ; #color_missing = 'k' ; # we hide the field where A<`imask_no_ice_pc`%
            l_pow_field=True ; pow_field=0.3 ;  vc_fld_powlog= [tmin, 0.01, 0.1, 0.5, tmax ]
            cunit = 'Total deformation [day$^{-1}$]'

        elif CWHAT in [ 'dadvs12t', 'dadvs12f' ]:
            cv_in = CWHAT  ; cv_out = CWHAT
            cpal_fld='RdBu_r' ; tmin=-500. ; tmax=-tmin ; df=100. ; color_top_cb='k'
            cunit = 'Prather increment for stress [Pa]'

        elif CWHAT in [ 'dadvdmgt', 'dadvdmgf' ]:
            cv_in = CWHAT  ; cv_out = CWHAT
            cpal_fld='RdBu_r' ; tmin=-0.005 ; tmax=-tmin ; df=0.005 ; color_top_cb='k'
            cunit = 'Prather increment for stress [Pa]'

        elif CWHAT in [ 'Qns', 'nshfls', 'qns_oce', 'qns_oce_si' ]:
            cv_in = CWHAT ; cv_out = 'qns_oce'
            imask_no_ice_pc = 5 ; # we hide the field where A<10%
            cpal_fld = 'ncview_parula_r'; color_top_cb='k'
            #cpal_fld = 'viridis_r'; color_top_cb='k'
            #cpal_fld = 'magma_r'; color_top_cb='k'
            tmin=-100.; tmax=0.; df=25.
            #l_pow_field=True; pow_field=5.
            vc_fld_powlog = [ tmin, -75., -50., -25., tmax ]
            cunit = r'Non-solar heat flux to the ocean [$W/m^{2}$]'

        elif CWHAT in [ 'qns_atmo' ]:
            cv_in = CWHAT ; cv_out = 'qns_atmo'
            imask_no_ice_pc = 5 ; # we hide the field where A<10%
            #cpal_fld = 'ncview_parula'; color_top_cb='k'
            #cpal_fld = 'viridis'; color_top_cb='k'
            #cpal_fld = 'inferno'; color_top_cb='k'
            cpal_fld = 'cmocean_haline'; color_top='w' ; color_top_cb='w' ; color_missing ='k'
            tmin=0.; tmax=100.; df=25.
            #l_pow_field=True; pow_field=5.
            vc_fld_powlog = [ tmin, 25., 50., 75., tmax ]
            cunit = r'Non-solar heat flux to the atmosphere [$W/m^{2}$]'


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

        self.vc_fld_force = vc_fld_force
        
        self.l_log_field   = l_log_field
        self.l_pow_field   = l_pow_field
        self.pow_field     = pow_field
        self.vc_fld_powlog = vc_fld_powlog

        self.imask_no_ice_pc = imask_no_ice_pc
        self.nm_ice_conc     = nm_ice_conc

        self.color_missing   = color_missing
        
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
        self.lSubstractMean = lSubstractMean
        self.l_smooth      = l_smooth
        self.nb_smooth     = nb_smooth


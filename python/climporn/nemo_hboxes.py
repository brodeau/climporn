from sys import exit

class nemo_hbox:

    ''' Bla bla bla
    '''
        
    def __init__(self, conf, box):
        self.conf = conf
        self.box  = box

    def size(self):
        config = self.conf
        if config == 'eNATL60':
            Ni0 = 8354
            Nj0 = 4729
        elif config == 'NATL60':
            Ni0 = 5422
            Nj0 = 3454
        elif config == 'eNATL4':
            Ni0 = 559
            Nj0 = 313
        elif config == 'ORCA36':
            Ni0 = 12962
            Nj0 = 9173
        elif config == 'EORCA12.L75':
            Ni0 = 4322
            Nj0 = 3147
        elif config == 'TROPICO05_NST':
            Ni0 = 297
            Nj0 = 156
        elif config == 'CALEDO10':
            Ni0 = 138
            Nj0 = 153
        elif config == 'CALEDO60':
            Ni0 = 788
            Nj0 = 853
        elif config == 'NANUK4':
            Ni0 = 492
            Nj0 = 566
        else:
            print('ERROR: we do not know NEMO config "'+str(config)+'" !')
            exit(0)
            #
        return  (Ni0,Nj0)

                    
    def idx(self):

        config    = self.conf
        (Ni0,Nj0) = self.size()
        box       = self.box

        rfact_zoom=1.
        
        l_show_exp = False
        l_show_cb  = True
        l_fill_holes_k=False
        l_annotate_name=False
        l_show_clock = True
        l_show_sign = False
        l_add_logo=True
        l_add_logo_ige=False
        l_add_logo_prc=False
        
        c_imshow_interp = 'none'

        l_add_quiver = False
        n_subsamp_qvr = 1

        cf_logo_on  = 'ocean-next_trans_white_281x191.png'
        cf_logo_ige = 'IGE_blanc_notext.png'
        cf_logo_prc = 'PRACE_blanc.png'

        x_clock = 1600 ; y_clock = 200
        x_sign = 1600 ; y_sign = 200
        x_logo = 2200  ; y_logo  = 50
        x_exp = 40     ; y_exp = 980


        
        
        if   [ config, box ] == [ 'eNATL60', 'ALL']:
            i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1440./float(Nj0) ; vcb=[0.61, 0.1, 0.36, 0.018]  ; font_rat=8.*rfact_zoom
            x_clock = 1600 ; y_clock = 200 ; x_logo=2200 ; y_logo=1200
            l_show_exp = True ; x_exp = 1750 ; y_exp = 300
            l_fill_holes_k=True

        if   [ config, box ] == [ 'eNATL60', 'ALL4K']:
            l_add_logo=False
            i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=2400./float(Nj0) ; vcb=[0.61, 0.1, 0.36, 0.018]  ; font_rat=8.*rfact_zoom
            x_clock = 2666 ; y_clock = 333
            l_show_exp = False ; x_exp = 2916 ; y_exp = 500
            l_show_sign = True ; x_sign = 3450 ; y_sign = 80
            l_fill_holes_k=True

        elif [ config, box ] == [ 'ORCA36', 'ALL']:
            i1=0   ; j1=900    ; i2=Ni0 ; j2=Nj0-982  ; rfact_zoom=1920./float(Ni0) ; vcb=[0.02, 0.91, 0.23, 0.018]  ; font_rat=12.*rfact_zoom
            l_show_clock=False; x_clock = 1600 ; y_clock = 200
            l_add_logo=False ; x_logo=2200 ; y_logo=1200
            l_show_exp = False ; x_exp = 1750 ; y_exp = 300
            l_fill_holes_k=False


            
        elif [ config, box ] == [ 'TROPICO05_NST', 'ALL']:
            i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=5 ; font_rat=0.2*rfact_zoom
            l_show_cb  = False ; vcb=[0.61, 0.1, 0.36, 0.018]
            l_show_clock=True ; x_clock = 120 ; y_clock = 110
            l_add_logo=True ; x_logo=520  ; y_logo=100; cf_logo_on  = 'ocean-next_trans_white_80x54.png'
            l_show_exp =False
            l_fill_holes_k=True
            #
        elif [ config, box ] == [ 'CALEDO10', 'ALL']:
            i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1 ; font_rat=8.*rfact_zoom
            l_show_cb  = False ; vcb=[0.61, 0.1, 0.36, 0.018]
            l_show_clock=False
            l_add_logo=False
            l_show_exp =False
            l_fill_holes_k=True

        elif [ config, box ] == [ 'CALEDO60', 'ALL']:
            i1=4 ; j1=18  ;  i2=Ni0-4 ; j2=Nj0-15  ; rfact_zoom=1 ; font_rat=1.5*rfact_zoom
            l_show_cb  = True ; vcb=[0.15, 0.07, 0.7, 0.018]
            l_show_clock=True ; x_clock = 50 ; y_clock = 800
            l_add_logo=False
            #l_annotate_name=True
            l_show_exp=False ; x_exp = 50 ; y_exp = 800
            l_fill_holes_k=True

        elif [ config, box ] == [ 'CALEDO60', 'ALL']:
            i1=4 ; j1=18  ;  i2=Ni0-4 ; j2=Nj0-15  ; rfact_zoom=1 ; font_rat=1.5*rfact_zoom
            l_show_cb  = True ; vcb=[0.15, 0.07, 0.7, 0.018]
            l_show_clock=True ; x_clock = 50 ; y_clock = 800
            l_add_logo=False
            #l_annotate_name=True
            l_show_exp=False ; x_exp = 50 ; y_exp = 800
            l_fill_holes_k=True

        elif [ config, box ] == [ 'NANUK4', 'ALL']:
            i1=0 ; j1=0  ;  i2=Ni0 ; j2=Nj0  ; rfact_zoom=1 ; font_rat=1.*rfact_zoom
            l_show_cb  = True ; vcb=[0.15, 0.08, 0.7, 0.018]
            l_show_clock=True ; x_clock = 320 ; y_clock = 540
            l_add_logo=False
            l_annotate_name=True
            l_show_exp=False ; x_exp = 50 ; y_exp = 80
            l_fill_holes_k=False

            
            
        elif [ config, box ] == [ 'ORCA36', 'ALLFRX']:
            # FR = FullRes !
            i1=0   ; j1=1000    ; i2=Ni0 ; j2=Nj0-1500  ; rfact_zoom=1. ; vcb=[0.02, 0.97, 0.23, 0.018]  ; font_rat=12.*rfact_zoom
            l_show_clock=False; x_clock = 1600 ; y_clock = 200
            l_add_logo=False ; x_logo=2200 ; y_logo=1200
            l_show_exp = False ; x_exp = 1750 ; y_exp = 300
            l_fill_holes_k=False
            
        elif   [ config, box ] == [ 'eNATL60', 'ALLC']:
            # pour la com!!!
            i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1440./float(Nj0) ; vcb=[0.61, 0.1, 0.36, 0.018]  ; font_rat=8.*rfact_zoom
            x_clock = 1600 ; y_clock = 200 ; l_add_logo=False ; l_add_logo_ige=False ; l_add_logo_prc=False
            l_show_exp = False ; x_exp = 1750 ; y_exp = 300
            l_fill_holes_k=True
        
        elif [ config, box ] == [ 'eNATL60', 'EUROPA']:
            i2=6400 ; j2=4000 ; i1=i2-2*1920; j1=j2-2*1080; rfact_zoom=0.5   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
            x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16
            
        elif [ config, box ] == [ 'eNATL60', 'EUROPAs']:
            i2=6500 ; j2=3600 ; i1=i2-1920; j1=j2-1080; rfact_zoom=1.  ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
            x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16
            
        elif   [ config, box ] == [ 'eNATL60', 'ALLFR']:
            # FR = FullRes !
            i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.39, 0.018]  ; font_rat=8.*rfact_zoom
            x_clock = 4000 ; y_clock = 200 ; x_logo = 6000 ; y_logo  = 50; 
            l_fill_holes_k=True
            l_show_clock=False; l_add_logo=False ; l_add_logo_ige=False ; l_add_logo_prc=False; l_show_exp=False
    
        elif   [ config, box ] == [ 'eNATL60', 'SALL']:
            i1=0   ; j1=0    ; i2=Ni0 ; j2=Nj0  ; rfact_zoom=1080./float(Nj0) ; vcb=[0.59, 0.1, 0.39, 0.018]  ; font_rat=8.*rfact_zoom
            x_clock = 1600 ; y_clock = 200 ; x_logo = 2200 ; y_logo  = 50
    
        elif [ config, box ] == [ 'eNATL60', 'Med']:
            i2=8040; i1=i2-2560 ; j1=1525 ; j2=j1+1440 ; rfact_zoom=1. ; vcb=[0.025, 0.06, 0.4, 0.02] ; font_rat=2.*rfact_zoom
            l_show_exp = True ; x_exp = 100 ; y_exp = 240
            x_clock = 200 ; y_clock = 170 ; x_logo = 1650 ; y_logo  = 1200
    
        elif [ config, box ] == [ 'eNATL60', 'Meddies']:
            i2=5800; j1=1400; i1=i2-2560 ; j2=j1+1440 ; rfact_zoom=1. ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
    
        elif [ config, box ] == [ 'eNATL60', 'MeddiesW']:
            i2=5680; j1=1200; i1=i2-2560 ; j2=j1+1440 ; rfact_zoom=1. ; vcb=[0.1, 0.06, 0.5, 0.02] ; font_rat=2.*rfact_zoom
            x_clock = 100 ; y_clock = 1320
            l_show_exp = True ; x_exp = 100 ; y_exp = 1380
            l_add_logo_ige=False ; l_add_logo_prc=False

        elif [ config, box ] == [ 'eNATL60', 'Manche']:
            l_add_quiver = True ; n_subsamp_qvr=20
            i2=6300; j2=3700; i1=i2-1920 ; j1=j2-1080 ; rfact_zoom=1. ; vcb=[0.1, 0.06, 0.5, 0.02] ; font_rat=2.*rfact_zoom
            #i2=6300; j2=3700; i1=i2-200 ; j1=j2-100 ; rfact_zoom=7. ; vcb=[0.1, 0.06, 0.5, 0.02] ; font_rat=2.*rfact_zoom
            x_clock = 100 ; y_clock = 980
            l_show_exp = True ; x_exp = 100 ; y_exp = 1030
            l_add_logo_ige=False ; l_add_logo_prc=False
            cf_logo_on  = 'ocean-next_trans_white_210x142.png' ; x_logo = 1670 ; y_logo  = 30
    
        elif [ config, box ] == [ 'eNATL60', 'Nordic']:
            i1=2200; j2=4715; i2=i1+2560 ; j1=j2-1440 ; rfact_zoom=1. ; vcb=[0.25, 0.06, 0.5, 0.02] ; font_rat=2.*rfact_zoom
            x_clock = 1010 ; y_clock = 1300
            l_show_exp = True ; x_exp = 1010 ; y_exp = 1340
            l_add_logo_ige=False ; l_add_logo_prc=False
            x_logo = 50 ; y_logo  = 50
    
        elif [ config, box ] == [ 'eNATL60', 'Med+BS']:
            i2=Ni0 ; i1=5400; j1=1530;  j2=3310 ; rfact_zoom=1080./float(j2-j1)   ; vcb=[0.025, 0.06, 0.4, 0.02] ; font_rat=3.*rfact_zoom
            l_add_logo_ige=False ; l_add_logo_prc=False
            rcc = 1080./1440. ; x_clock = 100.*rcc ; y_clock = 170.*rcc ; x_logo = 2090.*rcc ; y_logo  = 10.*rcc
            cf_logo_on  = 'ocean-next_trans_white_210x142.png'
    
        elif [ config, box ] == [ 'eNATL60', 'LabSea']:
            # 1818,3600 -> 3630,4722
            #i2=3770 ; j2=Nj0 ; i1=i2-1920; j1=j2-1200;  rfact_zoom=1200./float(j2-j1)   ; vcb=[0.015, 0.07, 0.28, 0.02] ; font_rat=3.*rfact_zoom
            i2=3770 ; j2=Nj0 ; i1=i2-1920; j1=j2-1080;  rfact_zoom=1080./float(j2-j1)   ; vcb=[0.014, 0.09, 0.25, 0.02] ; font_rat=3.*rfact_zoom
            l_add_logo_ige=False ; l_add_logo_prc=False
            x_clock = 40 ; y_clock = 170 ; x_logo = 2090 ; y_logo  = 10
    
        elif [ config, box ] == [ 'eNATL60', 'BlackSea']:
            i2=Ni0 ; j2=3330 ; i1=i2-1200; j1=j2-1200;  rfact_zoom=1.   ; vcb=[0.35, 0.06, 0.61, 0.02] ; font_rat=2.*rfact_zoom
            x_clock = 500 ; y_clock = 130 ; x_logo = 950 ; y_logo = 160 ; l_add_logo_ige=False ; l_add_logo_prc=False
            l_show_exp = True ; x_exp = 50 ; y_exp = 1100
            cf_logo_on  = 'ocean-next_trans_white_210x142.png'
    
        elif [ config, box ] == [ 'eNATL60', 'eBlackSea']:
            i2=Ni0 ; j2=3330 ; i1=i2-1920; j1=j2-1080;  rfact_zoom=1.   ; vcb=[0.5, 0.875, 0.485, 0.02] ; font_rat=2.*rfact_zoom
            x_clock = 30 ; y_clock = 1040 ; x_logo = 1500 ; y_logo = 16

        elif [ config, box ] == [ 'eNATL60', 'Antilles']:
            i1=1400; i2=i1+1920; j1=20; j2=j1+1080 ; rfact_zoom=1.   ; vcb=[0.018, 0.07, 0.44, 0.02] ; font_rat=2.*rfact_zoom
            x_clock = 250 ; y_clock = 130 ; x_logo = 1700 ; y_logo = 930 ; l_add_logo_ige=False ; l_add_logo_prc=False
            l_show_exp = True ; x_exp = 70 ; y_exp = 1020
            cf_logo_on  = 'ocean-next_trans_white_210x142.png'


        elif [ config, box ] == [ 'eNATL60', 'Bretagne']:
            l_add_quiver = True ; n_subsamp_qvr=10
            i1=5250; j1=2850; i2=i1+640 ; j2=j1+360 ; rfact_zoom=3 ; vcb=[0.53, 0.25, 0.44, 0.02] ; font_rat=0.8*rfact_zoom
            x_clock = 100 ; y_clock = 980
            l_show_exp = True ; x_exp = 100 ; y_exp = 1030
            l_add_logo_ige=False ; l_add_logo_prc=False
            cf_logo_on  = 'ocean-next_trans_white_210x142.png' ; x_logo = 1670 ; y_logo  = 30

        
        elif [ config, box ] == [ 'eNATL60', 'Brest']:
            l_add_quiver = True ; n_subsamp_qvr=3
            i2=5700; j2=3040; i1=i2-320 ; j1=j2-180 ; rfact_zoom=6 ; vcb=[0.6, 0.35, 0.38, 0.02] ; font_rat=0.4*rfact_zoom
            x_clock = 1200 ; y_clock = 550
            l_show_exp = True ; x_exp = 1200 ; y_exp = 600
            l_add_logo_ige=False ; l_add_logo_prc=False
            cf_logo_on  = 'ocean-next_trans_white_210x142.png' ; x_logo=1680 ; y_logo=450
            c_imshow_interp = 'bicubic'
    
        elif [ config, box ] == [ 'eNATL60', 'Portrait']:
            i1=2760; j1=1000; i2=4870; j2=4000 ; rfact_zoom=1.     ; vcb=[0.59, 0.1, 0.38, 0.018]  ; font_rat=1.*rfact_zoom
            l_show_clock=False
    
        elif [ config, box ] == [ 'eNATL60', 'EATL']:
            i1=3100; j1=2160; i2=i1+1920; j2=j1+1200 ; rfact_zoom=1. ; vcb=[0.3, 0.06, 0.38, 0.018] ; font_rat = 2.
            x_clock = 1570 ; y_clock = 1150 ; x_logo = 1620 ; y_logo = 16
        elif [ config, box ] == [ 'eNATL60', 'EATLcom']:
            i1=3100; j1=2160; i2=i1+1920; j2=j1+1200 ; rfact_zoom=1. ; vcb=[0.3, 0.06, 0.38, 0.018] ; font_rat = 2.
            l_add_logo=False ; l_add_logo_ige=False ; l_add_logo_prc=False
            l_show_exp = False ; l_fill_holes_k=True
            l_show_cb=False ; l_show_clock=False
    
        elif [ config, box ] == [ 'eNATL60', 'EATL2']:
            i1=2740; j1=1600; i2=i1+2560; j2=j1+1440 ; rfact_zoom=1. ; vcb=[0.3, 0.06, 0.38, 0.018] ; font_rat = 2.5
            x_clock = 2140 ; y_clock = 1390 ; x_logo = 2220 ; y_logo = 16
            l_show_exp = True
        elif [ config, box ] == [ 'eNATL60', 'GrlIcl']:
            i1=3370; j1=3941; i2=5062; j2=Nj0 ; rfact_zoom=1. ; vcb=[0.3, 0.1, 0.38, 0.018] ; font_rat = 2.
            x_clock = 1350 ; y_clock = 750 ; x_logo = 1400 ; y_logo = 16
    
        elif [ config, box ] == [ 'eNATL60', 'AzoresP']:
            # Azores Portrait:
            # 785 x 1190 => comparison of two => 1600 x 1200 (5px for frame), image is: 780x1190
            ## use: CMD="montage -geometry 780x1190+10+5 -background black <img1> <img2> <img_montage>"
            # => montage is then 1600x1200
            i2=4410; j2=2240 ; i1=i2-780; j1=j2-1190; rfact_zoom=1. ; vcb=[0.05, 0.05, 0.9, 0.018] ; font_rat = 2.
            l_add_logo=False; l_add_logo_prc=False; l_add_logo_ige=False
            x_clock = 400 ; y_clock = 120
    
        elif [ config, box ] == [ 'eNATL60', 'AzoresL']:
            # Azores Landscape: (1920x1080)
            i2=5200; j2=2240 ; i1=i2-1920; j1=j2-1080; rfact_zoom=1. ; vcb=[0.57, 0.08, 0.4, 0.018] ; font_rat = 2.
            l_add_logo=False; l_add_logo_prc=False; l_add_logo_ige=False
            x_clock = 1400 ; y_clock = 120 ; l_show_exp = True ; x_exp = 40 ; y_exp = 1030
    
        elif [ config, box ] == [ 'eNATL60', 'AzoresS']:
            # Azores small square:
            l_show_cb = True ; l_show_clock = False
            i2=4300; j2=2140 ; i1=i2-360; j1=j2-360; rfact_zoom=1. ; vcb=[0.05, 0.12, 0.9, 0.018] ; font_rat = 1.5
            l_add_logo=False; l_add_logo_prc=False; l_add_logo_ige=False
            x_clock = 400 ; y_clock = 120
    
        elif [ config, box ] == [ 'eNATL60', 'Band']:
            i1=5100-1920; j1=2200; i2=5100; j2=j1+1080 ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.38, 0.018] ; font_rat = 2.          
            l_show_clock = False ; l_add_logo = False ; #x_clock = 1420 ; y_clock = 1030 ; x_logo = 1500 ; y_logo = 16
    
        elif [ config, box ] == [ 'eNATL60', 'Balear']:
            i1=5750; j1=1880; i2=6470; j2=2600 ; rfact_zoom=1. ; vcb=[0.59, 0.1, 0.38, 0.018] ; font_rat = 2.
            x_clock = 1420 ; y_clock = 1030 ; x_logo = 1500 ; y_logo = 16
    
        elif [ config, box ] == [ 'eNATL60', 'GulfMex']:
            i1=0; j1=600; i2=1920; j2=j1+1080 ; rfact_zoom=1. ; vcb=[0.05, 0.94, 0.38, 0.018] ; font_rat = 2.
            x_clock = 50 ; y_clock = 40
            cf_logo_on  = 'ocean-next_trans_white_210x142.png' ; x_logo = 1680 ; y_logo = 930 ; l_add_logo_ige=False ; l_add_logo_prc=False
        elif [ config, box ] == [ 'eNATL60', 'GulfMexNL']:
            i1=0; j1=600; i2=1920; j2=j1+1080 ; rfact_zoom=1. ; vcb=[0.05, 0.94, 0.38, 0.018] ; font_rat = 2.
            x_clock = 50 ; y_clock = 40
            l_add_logo=False ; l_add_logo_ige=False ; l_add_logo_prc=False

        #elif [ config, box ] == [ 'eNATL60', 'Florida']:
        #    i1=0; j1=600; i2=1920; j2=j1+1080 ; rfact_zoom=1. ; vcb=[0.05, 0.94, 0.38, 0.018] ; font_rat = 2.
        #    x_clock = 50 ; y_clock = 40
        #    l_add_logo=False ; l_add_logo_ige=False ; l_add_logo_prc=False


        elif [ config, box ] == [ 'EORCA12.L75', 'SATL']:
            i1=2520 ; j1=700 ; i2=i1+1280 ; j2=j1+720 ; rfact_zoom=1 ; vcb=[0.2, 0.1, 0.6, 0.024]  ; font_rat=2.*rfact_zoom
            l_show_clock=True; x_clock = 110 ; y_clock = 680
            #l_add_logo=True ; x_logo=1200 ; y_logo= 20; cf_logo_on  = 'ocean-next_trans_white_140x95.png'
            l_add_logo=True ; x_logo=1114 ; y_logo= 610; cf_logo_on  = 'ocean-next_trans_white_140x95.png'
            l_show_exp = False ; x_exp = 1750 ; y_exp = 300
            l_fill_holes_k=False


        elif [ config, box ] == [ 'ORCA36', 'ALL']:
            i1=0   ; j1=900    ; i2=Ni0 ; j2=Nj0-982  ; rfact_zoom=1920./float(Ni0) ; vcb=[0.02, 0.91, 0.23, 0.018]  ; font_rat=12.*rfact_zoom
            l_show_clock=False; x_clock = 1600 ; y_clock = 200
            l_add_logo=False ; x_logo=2200 ; y_logo=1200
            l_show_exp = False ; x_exp = 1750 ; y_exp = 300
            l_fill_holes_k=False

    
        else:
            print(' ERROR: unknow box "'+box+'" for config "'+config+'" !!!')
            exit(0)


        self.rfact_zoom      = rfact_zoom
        self.vcb             = vcb
        self.font_rat        = font_rat
        self.l_annotate_name = l_annotate_name
        self.l_show_clock    = l_show_clock
        self.clock           = (x_clock,y_clock)
        self.l_show_sign    = l_show_sign
        self.sign           = (x_sign,y_sign)
        self.l_add_logo      = l_add_logo
        self.logo            = (x_logo,y_logo)
        self.l_add_logo_ige  = l_add_logo_ige
        self.l_add_logo_prc  = l_add_logo_prc
        self.l_show_exp      = l_show_exp
        self.exp             = (x_exp,y_exp)
        self.l_fill_holes_k  = l_fill_holes_k
        self.l_show_cb       = l_show_cb

        self.l_add_quiver    = l_add_quiver
        self.n_subsamp_qvr   = n_subsamp_qvr

        self.c_imshow_interp = c_imshow_interp

        self.cf_logo_on  = cf_logo_on # 
        self.cf_logo_ige = cf_logo_ige
        self.cf_logo_prc = cf_logo_prc


        return (i1,j1, i2,j2)


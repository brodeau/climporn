#!/usr/bin/env python

#     CLIMPORN

import sys

import nemo_hboxes as nhb


CNEMO = 'eNATL60'
VBOX  = ['ALL', 'EUROPA', 'EUROPAs', 'ALLFR', 'SALL', 'Med', 'Meddies', 'Med+BS', 'LabSea', 'BlackSea', 'Brittany', 'Portrait', 'EATL', 'EATL2', 'GrlIcl', 'AzoresP', 'AzoresL', 'AzoresS', 'Band', 'Balear']


for CBOX in VBOX:
          
    print '\n\n'


    print '\n Box "'+CBOX+'" for NEMO config "'+CNEMO+'" :'



    nemo_box = nhb.nemo_hbox(CNEMO,CBOX)


    (Ni0,Nj0) = nemo_box.size()
    print " "+CNEMO+": Ni0,Nj0 => ", Ni0,Nj0

    #sys.exit(0)


    (i1,j1, i2,j2) = nemo_box.idx()
    print " i1,j1, i2,j2 => ", i1,j1, i2,j2
    
    if nemo_box.l_show_clock:
        (x_clock,y_clock) = nemo_box.clock
        print ' x_clock,y_clock =', x_clock,y_clock

    if nemo_box.l_show_exp:
        (x_exp,y_exp) = nemo_box.exp
        print ' x_exp,y_exp =', x_exp,y_exp

    if nemo_box.l_add_logo:
        (x_logo,y_logo) = nemo_box.logo
        print ' x_logo,y_logo =', x_logo,y_logo


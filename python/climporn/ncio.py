import numpy as nmp
import sys

from netCDF4 import Dataset, num2date

from os import path
from os import stat


def ToEpochTime( vt, units, calendar ):
    '''
    # INPUT:
    #   * vt: time vector provided as something like: "days since ..."
    #   * units, calendar:
    #
    # OUTPUT:
    #   * vet: time vector converted to UNIX epoch time,
    #          aka "seconds since 1970-01-01 00:00:00"
    #          => as FLOAT !! not INTEGER !!!
    '''
    cfrmt = '%Y-%m-%d %H:%M:%S'
    #
    lvect = not (nmp.shape(vt)==())
    #
    if lvect:
        nt = len(vt)
        t0 = vt[0]
    else:
        t0 = vt
    print(' *** [ToEpochTime()]: original t0 as "'+units+'" => ', t0)
    t0d = num2date( t0, units, calendar )
    print(' *** [ToEpochTime()]: intitial date in datetime format => ', t0d)

    # We need to round this to the nearest second, because our target format is Epoch time (seconds since 1970)
    # and we want an integer!
    rdec = t0d.microsecond*1.E-6
    #print('LOLO: rmicros ', t0d.microsecond ) ; print('LOLO: rdec =', rdec)
    # t0 as "float" UNIX time:
    t0E = float( timegm( dtm.strptime( t0d.strftime(cfrmt) , cfrmt ).timetuple() ) + rdec )
    #print(' *** [ToEpochTime()]: intitial date in "Epoch Time" => ', t0E)

    if lvect:
        # we are not going to convert the whole array but instead:
        if   units[0:10] == 'days since':
            vdt = (vt[1:] - vt[0])*86400.
        elif units[0:11] == 'hours since':
            vdt = (vt[1:] - vt[0])*3600.
        elif units[0:13] == 'seconds since':
            vdt =  vt[1:] - vt[0]
        else:
            MsgExit('[ToEpochTime()] => unknown time unit: '+units)
        vet = nmp.zeros(nt)
        vet[0]  = t0E
        vet[1:] = t0E + vdt[:]
        del vdt
    else:
        vet = t0E
    #
    return vet



def GetTimeInfo( ncfile ):
    '''
    # Inspect time dimension
    # Get number of time-records, first and last date
    # Return them + dates as UNIX epoch time, aka "seconds since 1970-01-01 00:00:00" (float)
    '''
    if ldebug: print(' *** [GetTimeInfo()] Getting calendar/time info in '+ncfile+' ...')
    id_f = Dataset(ncfile)
    list_dim = id_f.dimensions.keys()
    list_var = id_f.variables.keys()
    for cd in [ 'time', 'time_counter', 'TIME', 'record', 't', 'none' ]:
        if cd in list_dim: break
    if cd == 'none': MsgExit('found no time-record dimension in file '+ncfile)
    if ldebug: print('   => time/record dimension is "'+cd+'"')
    nt = id_f.dimensions[cd].size
    cv = cd ; # ASSUMING VAR == DIM !!! Is it bad???
    if not cv in list_var: MsgExit('name of time variable is different than name of time dimension')
    clndr = id_f.variables[cv]
    dt1 = num2date( clndr[0], clndr.units, clndr.calendar ) ; dt2 = num2date( clndr[nt-1], clndr.units, clndr.calendar )
    rt1 = ToEpochTime( clndr[0],    clndr.units, clndr.calendar )
    rt2 = ToEpochTime( clndr[nt-1], clndr.units, clndr.calendar )    
    id_f.close()
    #
    if ldebug: print('   => first and last time records: ',dt1,'--',dt2,' (UNIX epoch: ', rt1,'--',rt2,')\n')
    #
    return nt, (rt1,rt2)



def GetTimeEpochVector( ncfile, kt1=0, kt2=0, isubsamp=1, lquiet=False ):
    '''
    # Get the time vector in the netCDF file, (from index kt1 to kt2, if these 2 are != 0!)
    # returns it as:
    # isubsamp: subsampling !!!
    # lquiet: shut the f* up!
    #  => ivt: time as UNIX epoch time, aka "seconds since 1970-01-01 00:00:00" (integer)
    '''
    ltalk = ( ldebug and not lquiet )
    cv_t_test = [ 'time', 'time_counter', 'TIME', 'record', 't', 'none' ]
    id_f = Dataset(ncfile)
    list_var = id_f.variables.keys()
    for cv in cv_t_test:
        if cv in list_var:
            clndr = id_f.variables[cv]
            cunt  = clndr.units
            ccal  = clndr.calendar
            if ltalk: print(' *** [GetTimeEpochVector()] reading "'+cv+'" in '+ncfile+' and converting it to Epoch time...')
            if kt1>0 and kt2>0:
                if kt1>=kt2: MsgExit('mind the indices when calling GetTimeEpochVector()')
                #vdate = num2date( clndr[kt1:kt2+1:isubsamp], clndr.units, clndr.calendar )
                vdate = clndr[kt1:kt2+1:isubsamp]
                cc = 'read...'
            else:
                #vdate = num2date( clndr[::isubsamp],         clndr.units, clndr.calendar )
                vdate = clndr[::isubsamp]
                cc = 'in TOTAL!'
            break
    id_f.close()
    if cv == 'none': MsgExit('found no time-record variable in file '+ncfile+' (possible fix: "cv_t_test" in "GetTimeEpochVector()")')
    #
    # Create the Unix Epoch time version:
    rvte = ToEpochTime( vdate, cunt, ccal )
    #
    if ltalk: print('   => '+str(len(rvte))+' records '+cc+'\n')
    return rvte






def wrt_1d_series(vt, vd, cvar, cinfo,
                  cu_t='unknown', cu_d='unknown', cln_d='unknown', nsmooth=0,
                  vd2=[], vd3=[], vd4=[], vd5=[],
                  cvar2='', cvar3='', cvar4='', cvar5='',
                  cln_d2='', cln_d3='', cln_d4='', cln_d5=''):

    cf_o  = cvar+'_'+cinfo+'.nc'

    lsmooth = False

    if nsmooth > 0:
        import clprn_stat as bs
        lsmooth = True

        if nsmooth == 11:
            vd_sm = bs.running_mean_11(vd, l_fill_bounds=False)
        elif nsmooth == 5:
            vd_sm = bs.running_mean_5(vd, l_fill_bounds=False)
        else:
            print('ERROR: wrt_1d_series.clprn_ncio => smoothing with nsmooth='+str(nsmooth)+' not supported!'); sys.exit(0)


    f_o = Dataset(cf_o, 'w', format='NETCDF4')

    nt = len(vt)
    if len(vd) != nt:  print('ERROR: wrt_1d_series.clprn_ncio => data & time have different lengths!'); sys.exit(0)

    l_do_v2=False ; l_do_v3=False ; l_do_v4=False ; l_do_v5=False
    if len(vd2) == nt: l_do_v2=True
    if len(vd3) == nt: l_do_v3=True
    if len(vd4) == nt: l_do_v4=True
    if len(vd5) == nt: l_do_v5=True


    f_o.createDimension('time', None)
    id_t = f_o.createVariable('time','f4',('time',), zlib=True, complevel=8) ;  id_t.units = cu_t

    id_d = f_o.createVariable(cvar,'f4',('time',), zlib=True, complevel=8)
    id_d.units = cu_d ;  id_d.long_name = cln_d

    if l_do_v2: id_d2 = f_o.createVariable(cvar2,'f4',('time',), zlib=True, complevel=8); id_d2.units = cu_d; id_d2.long_name = cln_d2
    if l_do_v3: id_d3 = f_o.createVariable(cvar3,'f4',('time',), zlib=True, complevel=8); id_d3.units = cu_d; id_d3.long_name = cln_d3
    if l_do_v4: id_d4 = f_o.createVariable(cvar4,'f4',('time',), zlib=True, complevel=8); id_d4.units = cu_d; id_d4.long_name = cln_d4
    if l_do_v5: id_d5 = f_o.createVariable(cvar5,'f4',('time',), zlib=True, complevel=8); id_d5.units = cu_d; id_d5.long_name = cln_d5


    if lsmooth:
        id_sm = f_o.createVariable(cvar+'_'+str(nsmooth)+'yrm','f4',('time',), zlib=True, complevel=8)
        id_sm.units = cu_d ;  id_sm.long_name = str(nsmooth)+'-year running mean of '+cln_d


    for jt in range(nt):
        id_t[jt]   = vt[jt]
        id_d[jt]   = vd[jt]
        if lsmooth: id_sm[jt] = vd_sm[jt]
        if l_do_v2: id_d2[jt] = vd2[jt]
        if l_do_v3: id_d3[jt] = vd3[jt]
        if l_do_v4: id_d4[jt] = vd4[jt]
        if l_do_v5: id_d5[jt] = vd5[jt]

    f_o.Author = 'L. Brodeau (clprn_ncio.py of Climporn)'
    f_o.close()
    print(' * wrt_1d_series => '+cf_o+' written!\n')

    return 0



def read_1d_series(cf_i, cv_i, cv_t='time', l_return_time=True):

    if not path.exists(cf_i): print('ERROR: read_1d_series.clprn_ncio => '+cf_i+' not found!'); sys.exit(0)

    id_i = Dataset(cf_i)
    list_variables = id_i.variables.keys()
    
    if l_return_time: vt = id_i.variables[cv_t][:]
    
    if cv_i in list_variables[:]:
        vd = id_i.variables[cv_i][:]
    else:
        print('WARNING: read_1d_series.clprn_ncio => '+cv_i+' not found into '+cf_i+' !')
        if l_return_time:
            vd = nmp.zeros(len(vt))
        else:
            vd = [ 0. ]
    id_i.close()

    if l_return_time:
        print('  * read_1d_series => just read '+cv_t+' and '+cv_i+' into '+cf_i+'\n')
        return vt, vd
    else:
        print('  * read_1d_series => just read '+cv_i+' into '+cf_i+'\n')
        return vd






def read_1d_series_2(cf_i, cv_i, cv_t='time'):

    if not path.exists(cf_i): print('ERROR: read_1d_series_2.clprn_ncio => '+cf_i+' not found!'); sys.exit(0)

    id_i = Dataset(cf_i)
    vt   = id_i.variables[cv_t][:]
    nt   = len(vt)
    xout = nmp.zeros((2,nt))
    xout[0,:] = vt[:]

    xout[1,:] = id_i.variables[cv_i][:]

    id_i.close()

    print('  * read_1d_series_2 => just read '+cv_t+' and '+cv_i+' into '+cf_i+'\n')
    return xout





def wrt_appnd_1d_series(vt, vd, cf, cvar1,  cu_t='year', cu_d='unknown', cln_d='unknown',
                        vd2=[], cvar2='', cln_d2='', cun2='unknown',
                        vd3=[], cvar3='', cln_d3='', cun3='unknown',
                        vd4=[], cvar4='', cln_d4='', cun4='unknown',
                        vd5=[], cvar5='', cln_d5='', cun5='unknown',
                        vd6=[], cvar6='', cln_d6='', cun6='unknown',
                        vd7=[], cvar7='', cln_d7='', cun7='unknown',
                        vd8=[], cvar8='', cln_d8='', cun8='unknown',
                        vd9=[], cvar9='', cln_d9='', cun9='unknown',
                        vd10=[],cvar10='',cln_d10='',cun10='unknown' ):

    l_nc_is_new = not path.exists(cf)

    Nt = len(vt)
    if len(vd) != Nt:
        print('ERROR: wrt_appnd_1d_series.clprn_ncio => data & time have different lengths!')
        sys.exit(0)

    l_do_v2=False ; l_do_v3=False ; l_do_v4=False ; l_do_v5=False
    l_do_v6=False ; l_do_v7=False ; l_do_v8=False ; l_do_v9=False ; l_do_v10=False
    if len(vd2) == Nt: l_do_v2=True
    if len(vd3) == Nt: l_do_v3=True
    if len(vd4) == Nt: l_do_v4=True
    if len(vd5) == Nt: l_do_v5=True
    if len(vd6) == Nt: l_do_v6=True
    if len(vd7) == Nt: l_do_v7=True
    if len(vd8) == Nt: l_do_v8=True
    if len(vd9) == Nt: l_do_v9=True
    if len(vd10)== Nt: l_do_v10=True


    if l_nc_is_new:
        f_out = Dataset(cf, 'w', format='NETCDF4')
    else:
        f_out = Dataset(cf, 'a', format='NETCDF4')

    if l_nc_is_new:
        jrec2write = 0
        f_out.createDimension('time', None)
        id_t = f_out.createVariable('time','f4',('time',), zlib=True, complevel=8) ;      id_t.units = 'year'

        id_v01   = f_out.createVariable(cvar1 ,'f4',('time',), zlib=True, complevel=8)
        id_v01.units     = cu_d
        id_v01.long_name = cln_d
        if l_do_v2:
            id_v02   = f_out.createVariable(cvar2 ,'f4',('time',), zlib=True, complevel=8)
            id_v02.units = cu_d
            if not cun2=='unknown': id_v02.units = cun2
            id_v02.long_name = cln_d2
        if l_do_v3:
            id_v03   = f_out.createVariable(cvar3 ,'f4',('time',), zlib=True, complevel=8)
            id_v03.units     = cu_d
            if not cun3=='unknown': id_v03.units = cun3
            id_v03.long_name = cln_d3
        if l_do_v4:
            id_v04   = f_out.createVariable(cvar4 ,'f4',('time',), zlib=True, complevel=8)
            id_v04.units     = cu_d
            if not cun4=='unknown': id_v04.units = cun4
            id_v04.long_name = cln_d4
        if l_do_v5:
            id_v05   = f_out.createVariable(cvar5 ,'f4',('time',), zlib=True, complevel=8)
            id_v05.units = cu_d
            if not cun5=='unknown': id_v05.units = cun5
            id_v05.long_name = cln_d5
        if l_do_v6:
            id_v06   = f_out.createVariable(cvar6 ,'f4',('time',), zlib=True, complevel=8)
            id_v06.units     = cu_d
            if not cun6=='unknown': id_v06.units = cun6
            id_v06.long_name = cln_d6
        if l_do_v7:
            id_v07   = f_out.createVariable(cvar7 ,'f4',('time',), zlib=True, complevel=8)
            id_v07.units     = cu_d
            if not cun7=='unknown': id_v07.units = cun7
            id_v07.long_name = cln_d7
        if l_do_v8:
            id_v08   = f_out.createVariable(cvar8 ,'f4',('time',), zlib=True, complevel=8)
            id_v08.units = cu_d
            if not cun8=='unknown': id_v08.units = cun8
            id_v08.long_name = cln_d8
        if l_do_v9:
            id_v09   = f_out.createVariable(cvar9 ,'f4',('time',), zlib=True, complevel=8)
            id_v09.units     = cu_d
            if not cun9=='unknown': id_v09.units = cun9
            id_v09.long_name = cln_d9
        if l_do_v10:
            id_v10   = f_out.createVariable(cvar10 ,'f4',('time',), zlib=True, complevel=8)
            id_v10.units     = cu_d
            if not cun7=='unknown': id_v10.units = cun10
            id_v10.long_name = cln_d10




        jrw = 0
        for jt in range(Nt):
            jrw = jrec2write + jt
            id_t[jrw]   = vt[jt]
            id_v01[jrw] = vd[jt]
            if l_do_v2: id_v02[jrw] = vd2[jt]
            if l_do_v3: id_v03[jrw] = vd3[jt]
            if l_do_v4: id_v04[jrw] = vd4[jt]
            if l_do_v5: id_v05[jrw] = vd5[jt]
            if l_do_v6: id_v06[jrw] = vd6[jt]
            if l_do_v7: id_v07[jrw] = vd7[jt]
            if l_do_v8: id_v08[jrw] = vd8[jt]
            if l_do_v9: id_v09[jrw] = vd9[jt]
            if l_do_v10:id_v10[jrw] = vd10[jt]
        f_out.about = 'Diagnostics created with Climporn (https://github.com/brodeau/climporn)'
    else:
        vtr = f_out.variables['time']
        jrec2write = len(vtr)
        v01 = f_out.variables[cvar1]
        if l_do_v2: v02 = f_out.variables[cvar2]
        if l_do_v3: v03 = f_out.variables[cvar3]
        if l_do_v4: v04 = f_out.variables[cvar4]
        if l_do_v5: v05 = f_out.variables[cvar5]
        if l_do_v6: v06 = f_out.variables[cvar6]
        if l_do_v7: v07 = f_out.variables[cvar7]
        if l_do_v8: v08 = f_out.variables[cvar8]
        if l_do_v9: v09 = f_out.variables[cvar9]
        if l_do_v10:v10 = f_out.variables[cvar10]
        jrw = 0
        for jt in range(Nt):
            jrw = jrec2write + jt
            vtr[jrw]  = vt[jt]
            v01[jrw] = vd[jt]
            if l_do_v2: v02[jrw] = vd2[jt]
            if l_do_v3: v03[jrw] = vd3[jt]
            if l_do_v4: v04[jrw] = vd4[jt]
            if l_do_v5: v05[jrw] = vd5[jt]
            if l_do_v6: v06[jrw] = vd6[jt]
            if l_do_v7: v07[jrw] = vd7[jt]
            if l_do_v8: v08[jrw] = vd8[jt]
            if l_do_v9: v09[jrw] = vd9[jt]
            if l_do_v10:v10[jrw] = vd10[jt]
    f_out.close()

    print('   *** '+str(Nt)+' snapshots of fields saved into '+cf+' !\n')




# Write a Array(x,y,t) to netcdf:
def wrt_appnd_2dt_series(vx, vy, vt, xd, cf, cvar1, missing_val=-9999.,
                         cxdnm='x', cydnm='y', cxvnm='lon', cyvnm='lat',
                         cu_t='year', cu_d='unknown', cln_d='unknown',
                         xd2=[], cvar2='', cln_d2='', cun2='unknown',
                         xd3=[], cvar3='', cln_d3='', cun3='unknown',
                         xd4=[], cvar4='', cln_d4='', cun4='unknown',
                         xd5=[], cvar5='', cln_d5='', cun5='unknown',
                         xd6=[], cvar6='', cln_d6='', cun6='unknown',
                         xd7=[], cvar7='', cln_d7='', cun7='unknown',
                         xd8=[], cvar8='', cln_d8='', cun8='unknown',
                         xd9=[], cvar9='', cln_d9='', cun9='unknown',
                         xd10=[],cvar10='',cln_d10='',cun10='unknown' ):

    l_nc_is_new = not path.exists(cf)
    
    [ Nt, ny, nx ] = nmp.shape(xd)

    if Nt != len(vt):
        print('ERROR: wrt_appnd_2dt_series.clprn_ncio => data & time have different lengths!')
        sys.exit(0)
        
    l_do_v2=False ; l_do_v3=False ; l_do_v4=False ; l_do_v5=False
    l_do_v6=False ; l_do_v7=False ; l_do_v8=False ; l_do_v9=False ; l_do_v10=False
    if nmp.shape(xd2) == ( Nt, ny, nx ): l_do_v2=True
    if nmp.shape(xd3) == ( Nt, ny, nx ): l_do_v3=True
    if nmp.shape(xd4) == ( Nt, ny, nx ): l_do_v4=True
    if nmp.shape(xd5) == ( Nt, ny, nx ): l_do_v5=True
    if nmp.shape(xd6) == ( Nt, ny, nx ): l_do_v6=True
    if nmp.shape(xd7) == ( Nt, ny, nx ): l_do_v7=True
    if nmp.shape(xd8) == ( Nt, ny, nx ): l_do_v8=True
    if nmp.shape(xd9) == ( Nt, ny, nx ): l_do_v9=True
    if nmp.shape(xd10)== ( Nt, ny, nx ): l_do_v10=True


    if l_nc_is_new:
        f_out = Dataset(cf, 'w', format='NETCDF4')
    else:
        f_out = Dataset(cf, 'a', format='NETCDF4')

    if l_nc_is_new:
        jrec2write = 0

        f_out.createDimension('time', None)
        f_out.createDimension(cydnm , ny)
        f_out.createDimension(cxdnm , nx)

        id_t   = f_out.createVariable('time','f4',('time',), zlib=True, complevel=8) ; id_t.units = 'year'
        id_lat = f_out.createVariable(cyvnm ,'f4',(cydnm,), zlib=True, complevel=8)
        id_lon = f_out.createVariable(cxvnm ,'f4',(cxdnm,), zlib=True, complevel=8)

        id_lat[:] = vy[:]
        id_lon[:] = vx[:]
        
        id_x01   = f_out.createVariable(cvar1 ,'f4',('time',cydnm,cxdnm,), fill_value=missing_val, zlib=True, complevel=8)
        id_x01.units     = cu_d
        id_x01.long_name = cln_d
        if l_do_v2:
            id_x02   = f_out.createVariable(cvar2 ,'f4',('time',cydnm,cxdnm,), fill_value=missing_val, zlib=True, complevel=8)
            id_x02.units = cu_d
            if not cun2=='unknown': id_x02.units = cun2
            id_x02.long_name = cln_d2
        if l_do_v3:
            id_x03   = f_out.createVariable(cvar3 ,'f4',('time',cydnm,cxdnm,), fill_value=missing_val, zlib=True, complevel=8)
            id_x03.units     = cu_d
            if not cun3=='unknown': id_x03.units = cun3
            id_x03.long_name = cln_d3
        if l_do_v4:
            id_x04   = f_out.createVariable(cvar4 ,'f4',('time',cydnm,cxdnm,), fill_value=missing_val, zlib=True, complevel=8)
            id_x04.units     = cu_d
            if not cun4=='unknown': id_x04.units = cun4
            id_x04.long_name = cln_d4
        if l_do_v5:
            id_x05   = f_out.createVariable(cvar5 ,'f4',('time',cydnm,cxdnm,), fill_value=missing_val, zlib=True, complevel=8)
            id_x05.units = cu_d
            if not cun5=='unknown': id_x05.units = cun5
            id_x05.long_name = cln_d5
        if l_do_v6:
            id_x06   = f_out.createVariable(cvar6 ,'f4',('time',cydnm,cxdnm,), fill_value=missing_val, zlib=True, complevel=8)
            id_x06.units     = cu_d
            if not cun6=='unknown': id_x06.units = cun6
            id_x06.long_name = cln_d6
        if l_do_v7:
            id_x07   = f_out.createVariable(cvar7 ,'f4',('time',cydnm,cxdnm,), fill_value=missing_val, zlib=True, complevel=8)
            id_x07.units     = cu_d
            if not cun7=='unknown': id_x07.units = cun7
            id_x07.long_name = cln_d7
        if l_do_v8:
            id_x08   = f_out.createVariable(cvar8 ,'f4',('time',cydnm,cxdnm,), fill_value=missing_val, zlib=True, complevel=8)
            id_x08.units = cu_d
            if not cun8=='unknown': id_x08.units = cun8
            id_x08.long_name = cln_d8
        if l_do_v9:
            id_x09   = f_out.createVariable(cvar9 ,'f4',('time',cydnm,cxdnm,), fill_value=missing_val, zlib=True, complevel=8)
            id_x09.units     = cu_d
            if not cun9=='unknown': id_x09.units = cun9
            id_x09.long_name = cln_d9
        if l_do_v10:
            id_x10   = f_out.createVariable(cvar10 ,'f4',('time',cydnm,cxdnm,), fill_value=missing_val, zlib=True, complevel=8)
            id_x10.units     = cu_d
            if not cun7=='unknown': id_x10.units = cun10
            id_x10.long_name = cln_d10

        jrw = 0
        for jt in range(Nt):
            jrw = jrec2write + jt
            id_t[jrw]   = vt[jt]
            id_x01[jrw,:,:] = xd[jt,:,:]
            if l_do_v2: id_x02[jrw,:,:] = xd2[jt,:,:]
            if l_do_v3: id_x03[jrw,:,:] = xd3[jt,:,:]
            if l_do_v4: id_x04[jrw,:,:] = xd4[jt,:,:]
            if l_do_v5: id_x05[jrw,:,:] = xd5[jt,:,:]
            if l_do_v6: id_x06[jrw,:,:] = xd6[jt,:,:]
            if l_do_v7: id_x07[jrw,:,:] = xd7[jt,:,:]
            if l_do_v8: id_x08[jrw,:,:] = xd8[jt,:,:]
            if l_do_v9: id_x09[jrw,:,:] = xd9[jt,:,:]
            if l_do_v10:id_x10[jrw,:,:] = xd10[jt,:,:]
        f_out.about = 'Diagnostics created with Climporn (https://github.com/brodeau/climporn)'
        #
        #
    else:
        vtr = f_out.variables['time']
        jrec2write = len(vtr)
        x01 = f_out.variables[cvar1]
        if l_do_v2: x02 = f_out.variables[cvar2]
        if l_do_v3: x03 = f_out.variables[cvar3]
        if l_do_v4: x04 = f_out.variables[cvar4]
        if l_do_v5: x05 = f_out.variables[cvar5]
        if l_do_v6: x06 = f_out.variables[cvar6]
        if l_do_v7: x07 = f_out.variables[cvar7]
        if l_do_v8: x08 = f_out.variables[cvar8]
        if l_do_v9: x09 = f_out.variables[cvar9]
        if l_do_v10:x10 = f_out.variables[cvar10]
        jrw = 0
        for jt in range(Nt):
            jrw = jrec2write + jt
            vtr[jrw]  = vt[jt]
            x01[jrw,:,:] = xd[jt,:,:]
            if l_do_v2: x02[jrw,:,:] = xd2[jt,:,:]
            if l_do_v3: x03[jrw,:,:] = xd3[jt,:,:]
            if l_do_v4: x04[jrw,:,:] = xd4[jt,:,:]
            if l_do_v5: x05[jrw,:,:] = xd5[jt,:,:]
            if l_do_v6: x06[jrw,:,:] = xd6[jt,:,:]
            if l_do_v7: x07[jrw,:,:] = xd7[jt,:,:]
            if l_do_v8: x08[jrw,:,:] = xd8[jt,:,:]
            if l_do_v9: x09[jrw,:,:] = xd9[jt,:,:]
            if l_do_v10:x10[jrw,:,:] = xd10[jt,:,:]
    f_out.close()

    print('   *** '+str(Nt)+' snapshots of fields saved into '+cf+' !\n')


def write_2d_mask( cf_out, MSK, xlon=[], xlat=[], name='mask', clon='nav_lon', clat='nav_lat' ):

    (nj,ni) = nmp.shape(MSK)

    f_out = Dataset(cf_out, 'w', format='NETCDF4')

    # Dimensions:
    f_out.createDimension('x', ni)
    f_out.createDimension('y', nj)

    if (xlon != []) and (xlat != []):
        if (xlon.shape == (nj,ni)) and (xlon.shape == xlat.shape):
            id_lon  = f_out.createVariable(clon ,'f4',('y','x',), zlib=True, complevel=8)
            id_lat  = f_out.createVariable(clat ,'f4',('y','x',), zlib=True, complevel=8)
            id_lon[:,:] = xlon[:,:]
            id_lat[:,:] = xlat[:,:]
        
    id_msk  = f_out.createVariable(name ,'i1',('y','x',), zlib=True, complevel=8)
    id_msk[:,:] = MSK[:,:]

    f_out.about = 'Diagnostics created with Climporn (https://github.com/brodeau/climporn)'
    f_out.close()

    return


def dump_2d_field( cf_out, XFLD, xlon=[], xlat=[], name='field', unit='', long_name='', mask=[], clon='nav_lon', clat='nav_lat' ):
    #
    (nj,ni) = nmp.shape(XFLD)

    f_out = Dataset(cf_out, 'w', format='NETCDF4')

    # Dimensions:
    f_out.createDimension('y', nj)
    f_out.createDimension('x', ni)

    if (xlon != []) and (xlat != []):
        if (xlon.shape == (nj,ni)) and (xlon.shape == xlat.shape):
            id_lon  = f_out.createVariable(clon ,'f4',('y','x',), zlib=True, complevel=8)
            id_lat  = f_out.createVariable(clat ,'f4',('y','x',), zlib=True, complevel=8)
            id_lon[:,:] = xlon[:,:]
            id_lat[:,:] = xlat[:,:]
        
    id_fld  = f_out.createVariable(name ,'f4',('y','x',), zlib=True, complevel=8)
    if long_name != '': id_fld.long_name = long_name
    if unit      != '': id_fld.units     = unit

    if mask != []:
        xtmp = nmp.zeros((nj,ni))
        xtmp[:,:] = XFLD[:,:]
        idx_land = nmp.where( mask < 0.5)
        xtmp[idx_land] = nmp.nan
        id_fld[:,:] = xtmp[:,:]
        del xtmp
    else:
        id_fld[:,:] = XFLD[:,:]

    f_out.about = 'Diagnostics created with Climporn (https://github.com/brodeau/climporn)'
    f_out.close()

    return





def Dump2DMultiField( cf_out, XFLD, vnames, vlnames=[], vunits=[], vndim=[], xlon=[], xlat=[], vtime=[], \
                      clon='nav_lon', clat='nav_lat', dim_nm=['lat','lon'], cFieldPres='f8',
                      rfillval=None, cabout='' ):
    #
    if len(vtime)>0:
        l_add_time = True
        (nbfld, Nbt, nj, ni) = nmp.shape(XFLD)
        if Nbt != len(vtime): print('ERROR (Dump2DMultiField): array and time vector disagree!'); sys.exit(0)
    else:
        l_add_time = False
        (nbfld,      nj, ni) = nmp.shape(XFLD)

    l_add_long_names = ( len(vlnames) == nbfld )
    l_add_units      = ( len(vunits)  == nbfld )

    vnbdim = nmp.zeros(nbfld)
    if vndim == []:
        vnbdim[:] = 3 ; # default dim is 3 (time_counter,y,x)
    else:
        nn0 = len(vndim)
        if nbfld != nn0: print('ERROR (Dump2DMultiField): vndim and main array dont agree in shape!'); sys.exit(0)
        vnbdim[:] = vndim[:]
        
    nf = len(vnames)
    if nbfld != nf: print('ERROR (Dump2DMultiField): list of names and main array dont agree in shape!'); sys.exit(0)

    f_out = Dataset(cf_out, 'w', format='NETCDF4')

    l_coord_2d = False
    cnm_dim_y = dim_nm[0]
    cnm_dim_x = dim_nm[1]
    
    if (xlon != []) and (xlat != []):
        if (xlon.shape == (nj,ni)) and (xlon.shape == xlat.shape):
            l_coord_2d = True
            cnm_dim_x = 'x'
            cnm_dim_y = 'y'    
    
    # Dimensions:    
    f_out.createDimension(cnm_dim_x, ni)
    f_out.createDimension(cnm_dim_y, nj)
    if l_add_time: f_out.createDimension('time_counter', None)

    if (xlon != []) and (xlat != []):
        if l_coord_2d:
            id_lon  = f_out.createVariable(clon ,'f4',(cnm_dim_y,cnm_dim_x,), zlib=True, complevel=8)
            id_lat  = f_out.createVariable(clat ,'f4',(cnm_dim_y,cnm_dim_x,), zlib=True, complevel=8)
            id_lon[:,:] = xlon[:,:]
            id_lat[:,:] = xlat[:,:]
        else:
            id_lon  = f_out.createVariable(cnm_dim_x ,'f4',(cnm_dim_x,), zlib=True, complevel=8)
            id_lat  = f_out.createVariable(cnm_dim_y ,'f4',(cnm_dim_y,), zlib=True, complevel=8)
            id_lon[:] = xlon[:]
            id_lat[:] = xlat[:]

    if l_add_time:
        id_tim    = f_out.createVariable('time_counter' ,'f8',('time_counter',))
        id_tim[:] = vtime[:]
        
    #id_fld = nmp.zeros(nbfld, dtype=int)
    for jv in range(nbfld):
        if (not l_add_time) or (vnbdim[jv]==2):
            id_fld  = f_out.createVariable(vnames[jv] ,cFieldPres,(cnm_dim_y,cnm_dim_x,), fill_value=rfillval, zlib=True, complevel=8)
            if l_add_long_names: id_fld.long_name = vlnames[jv]
            if l_add_units:      id_fld.units     =  vunits[jv]
            #
            if l_add_time:
                id_fld[:,:] = XFLD[jv,0,:,:]
            else:
                id_fld[:,:] = XFLD[jv,:,:]
        else:
            id_fld  = f_out.createVariable(vnames[jv] ,cFieldPres,('time_counter',cnm_dim_y,cnm_dim_x,), fill_value=rfillval, zlib=True, complevel=8)
            id_fld[:,:,:] = XFLD[jv,:,:,:]
            
    if cabout == '':
        f_out.about = 'File generated with `ncio.Dump2DMultiField()` of `climporn` (https://github.com/brodeau/climporn)'
    else:
        f_out.about = cabout+' [ + `ncio.Dump2DMultiField()` of `climporn` ]'
    
    f_out.close()

    return






def dump_3d_multi_field( cf_out, XFLD, vnames, vndim=[], xlon=[], xlat=[], vdepth=[], vtime=[], \
                         clon='nav_lon', clat='nav_lat', dim_nm=['depth','lat','lon'], cFieldPres='f8' ):
    
    if len(vtime)>0:
        l_add_time = True
        (nbfld, Nbt, nk, nj, ni) = nmp.shape(XFLD)
        if Nbt != len(vtime): print('ERROR (dump_3d_multi_field): array and time vector disagree!'); sys.exit(0)
    else:
        l_add_time = False
        (nbfld,      nk, nj, ni) = nmp.shape(XFLD)

    vnbdim = nmp.zeros(nbfld)
    if vndim == []:
        vnbdim[:] = 4 ; # default dim is 4 (time_counter,z,y,x)
    else:
        nn0 = len(vndim)
        if nbfld != nn0: print('ERROR (dump_3d_multi_field): vndim and main array dont agree in shape!'); sys.exit(0)
        vnbdim[:] = vndim[:]
        
    nf = len(vnames)
    if nbfld != nf: print('ERROR (dump_3d_multi_field): list of names and main array dont agree in shape!'); sys.exit(0)

    f_out = Dataset(cf_out, 'w', format='NETCDF4')


    l_coord_2d = False
    cnm_dim_z = dim_nm[0]
    cnm_dim_y = dim_nm[2]
    cnm_dim_x = dim_nm[3]

    
    if (xlon != []) and (xlat != []):
        if (xlon.shape == (nj,ni)) and (xlon.shape == xlat.shape):
            l_coord_2d = True
            cnm_dim_x = 'x'
            cnm_dim_y = 'y'    

    cnm_dim_z = 'z'
                    
    # Dimensions:    
    f_out.createDimension(cnm_dim_x, ni)
    f_out.createDimension(cnm_dim_y, nj)
    f_out.createDimension(cnm_dim_z, nk)
    if l_add_time: f_out.createDimension('time_counter', None)

    if (xlon != []) and (xlat != []):
        if l_coord_2d:
            id_lon  = f_out.createVariable(clon ,'f4',(cnm_dim_y,cnm_dim_x,), zlib=True, complevel=8)
            id_lat  = f_out.createVariable(clat ,'f4',(cnm_dim_y,cnm_dim_x,), zlib=True, complevel=8)            
            id_lon[:,:] = xlon[:,:]
            id_lat[:,:] = xlat[:,:]
        else:
            id_lon  = f_out.createVariable(cnm_dim_x ,'f4',(cnm_dim_x,), zlib=True, complevel=8)
            id_lat  = f_out.createVariable(cnm_dim_y ,'f4',(cnm_dim_y,), zlib=True, complevel=8)
            id_lon[:] = xlon[:]
            id_lat[:] = xlat[:]

    if (vdepth != []):
        id_dep  = f_out.createVariable('deptht'  ,'f4',(cnm_dim_z,),           zlib=True, complevel=8)
        id_dep[:]   = vdepth[:]
        
    if l_add_time:
        id_tim    = f_out.createVariable('time_counter' ,'f8',('time_counter',))
        id_tim[:] = vtime[:]
        
    #id_fld = nmp.zeros(nbfld, dtype=int)
    for jv in range(nbfld):
        if (not l_add_time) or (vnbdim[jv]==3):
            id_fld  = f_out.createVariable(vnames[jv] ,cFieldPres,(cnm_dim_z,cnm_dim_y,cnm_dim_x,), zlib=True, complevel=8)
            if l_add_time:
                id_fld[:,:,:] = XFLD[jv,0,:,:,:]
            else:
                id_fld[:,:,:] = XFLD[jv,:,:,:]
        else:
            id_fld  = f_out.createVariable(vnames[jv] ,cFieldPres,('time_counter',cnm_dim_z,cnm_dim_y,cnm_dim_x,), zlib=True, complevel=8)
            id_fld[:,:,:,:] = XFLD[jv,:,:,:,:]
            
    f_out.about = 'Diagnostics created with Climporn (https://github.com/brodeau/climporn)'
    f_out.close()

    return



def SaveTimeSeries( ivt, xd, vvar, ncfile, time_units='unknown', vunits=[], vlnm=[], missing_val=-9999. ):
    '''
    #  * ivt: time vector of length Nt, unit: UNIX Epoch time            [integer]
    #         => aka "seconds since 1970-01-01 00:00:00"
    #  *  xd: 2D numpy array that contains Nf time series of length Nt   [real]
    #          => hence of shape (Nf,Nt)
    #  * vvar: vector of length Nf of the Nf variable names                         [string]
    #  * vunits, vlnm: vectors of length Nf of the Nf variable units and long names [string]
    #  * missing_val: value for missing values...                        [real]
    '''
    (Nf,Nt) = xd.shape
    if len(ivt) != Nt: MsgExit('SaveTimeSeries() => disagreement in the number of records between "ivt" and "xd"')
    if len(vvar)!= Nf: MsgExit('SaveTimeSeries() => disagreement in the number of fields between "vvar" and "xd"')
    l_f_units = (nmp.shape(vunits)==(Nf,)) ; l_f_lnm = (nmp.shape(vlnm)==(Nf,))
    #
    print('\n *** About to write file "'+ncfile+'"...')
    f_o = Dataset(ncfile, 'w', format='NETCDF4')
    f_o.createDimension('time', None)
    id_t = f_o.createVariable('time','f8',('time',))
    id_t.calendar = 'gregorian' ; id_t.units = time_units
    #
    id_d = []
    for jf in range(Nf):
        id_d.append( f_o.createVariable(vvar[jf],'f4',('time',), fill_value=missing_val, zlib=True, complevel=8) )
        if l_f_units: id_d[jf].units   = vunits[jf]
        if l_f_lnm:   id_d[jf].long_name = vlnm[jf]
    #
    print('   ==> writing "time"')
    id_t[:] = ivt.astype(nmp.float64)
    for jf in range(Nf):
        print('   ==> writing "'+vvar[jf]+'"')
        id_d[jf][:] = xd[jf,:]
    f_o.about = cabout_nc
    f_o.close()
    print(' *** "'+ncfile+'" successfully written!\n')
    return 0

#!/usr/bin/env python

import sys
import os
import numpy as nmp
from netCDF4 import Dataset

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors


fig_type='png'


cf_file_nemo = '/home/brodeau/tmp/test_python_nc/Meddies_votemper_eNATL60-BLBT02X_1h_20100610_20100624_gridT_20100619-20100619_jk158.nc'

# Opens and reads stuff into the netcdf file:

id_file = Dataset(cf_file_nemo)

list_var = id_file.variables.keys()

print '\n Variables in the files:'
for cvar in list_var:
    nb_dim = len(id_file.variables[cvar].dimensions)
    print ' *** '+cvar+' -> has '+str(nb_dim)+' dimensions'
print ''
    

if 'time_counter' in list_var:    
    vtime = id_file.variables['time_counter'][:]
    Nt = len(vtime)
    print '\n There is a "time_counter" in file '+cf_file_nemo+' !'
    print '   => '+str(Nt)+' snapshots!'
else:
    print '\n There is NO "time_counter" in file '+cf_file_nemo+' !'
    Nt = 0

# Read the depth vector:
vdepth = id_file.variables['deptht'][:]

print '\n Depth array:'
for jk in range(len(vdepth)):
    print ' Level #'+str(jk)+' => '+str(vdepth[jk])+' m'


# Read the 2D longitude and latitude:
Xlon = id_file.variables['nav_lon'][:,:]
Xlat = id_file.variables['nav_lat'][:,:]

# In case you would read the entire 4D (3D+time) temperature array, you would do: 
#Xtemp = id_file.variables['votemper'][:,:,:,:] ; # Xtemp is 4D array :(

# Read only hour # 12:
#Xtemp = id_file.variables['votemper'][11,:,:,:] ; # Xtemp is 3D array

# Read only hour # 12 at depth 1000m:
jk_1000 = 106
Xtemp = id_file.variables['votemper'][11,jk_1000,:,:] ; # Xtemp is 2D array
    
id_file.close()

(Nj,Ni) = nmp.shape(Xlon)

#############################


fig = plt.figure(num = 1, figsize=(8,8), dpi=None, facecolor='w', edgecolor='0.5')
ax  = plt.axes([0., 0., 1., 1.], axisbg = '0.5')

cf = plt.imshow(Xtemp[:,:], interpolation='none')

plt.axis([ 0, Ni, 0, Nj])

#plt.title('NEMO: '+cfield+', coupled '+CNEMO+', '+cday+' '+chour+':00', **cfont_title)


#ax2 = plt.axes([0.3, 0.08, 0.4, 0.025])


plt.savefig('figure.png', dpi=100, orientation='portrait', facecolor='k')
plt.close(1)


# when on unix, doing "ncdum -h Meddies_votemper_eNATL60-BLBT02X_1h_20100610_20100624_gridT_20100619-20100619_jk158"
#
#netcdf Meddies_votemper_eNATL60-BLBT02X_1h_20100610_20100624_gridT_20100619-20100619_jk158 {
#    dimensions:
#            deptht = 158 ;
#            axis_nbounds = 2 ;
#            y = 901 ;
#            x = 611 ;
#            time_counter = UNLIMITED ; // (24 currently)
#    variables:
#            float deptht(deptht) ;
#                    deptht:name = "deptht" ;
#                    deptht:long_name = "Vertical T levels" ;
#                    deptht:units = "m" ;
#                    deptht:axis = "Z" ;
#                    deptht:positive = "down" ;
#                    deptht:bounds = "deptht_bounds" ;
#            float deptht_bounds(deptht, axis_nbounds) ;
#            float nav_lat(y, x) ;
#                    nav_lat:standard_name = "latitude" ;
#                    nav_lat:long_name = "Latitude" ;
#                    nav_lat:units = "degrees_north" ;
#            float nav_lon(y, x) ;
#                    nav_lon:standard_name = "longitude" ;
#                    nav_lon:long_name = "Longitude" ;
#                    nav_lon:units = "degrees_east" ;
#            double time_centered(time_counter) ;
#                    time_centered:standard_name = "time" ;
#                    time_centered:long_name = "Time axis" ;
#                    time_centered:calendar = "gregorian" ;
#                    time_centered:units = "seconds since 1900-01-01 00:00:00" ;
#                    time_centered:time_origin = "1900-01-01 00:00:00" ;
#                    time_centered:bounds = "time_centered_bounds" ;
#            double time_centered_bounds(time_counter, axis_nbounds) ;
#            double time_counter(time_counter) ;
#                    time_counter:axis = "T" ;
#                    time_counter:standard_name = "time" ;
#                    time_counter:long_name = "Time axis" ;
#                    time_counter:calendar = "gregorian" ;
#                    time_counter:units = "seconds since 1900-01-01 00:00:00" ;
#                    time_counter:time_origin = "1900-01-01 00:00:00" ;
#                    time_counter:bounds = "time_counter_bounds" ;
#            double time_counter_bounds(time_counter, axis_nbounds) ;
#            float votemper(time_counter, deptht, y, x) ;
#                    votemper:standard_name = "sea_water_potential_temperature" ;
#                    votemper:long_name = "temperature" ;
#                    votemper:units = "degC" ;
#                    votemper:online_operation = "average" ;
#                    votemper:interval_operation = "40 s" ;
#                    votemper:interval_write = "1 h" ;
#                    votemper:cell_methods = "time: mean (interval: 40 s)" ;
#                    votemper:_FillValue = 1.e+20f ;
#                    votemper:missing_value = 1.e+20f ;
#                    votemper:coordinates = "time_centered deptht nav_lat nav_lon" ;
#
#    // global attributes:
#                    :name = "/scratch/tmp/5425687/eNATL60-BLBT02X_1h_20100610_20100704_gridT" ;
#                    :description = "ocean T grid variables" ;
#                    :title = "ocean T grid variables" ;
#                    :Conventions = "CF-1.6" ;
#                    :timeStamp = "2019-Apr-11 18:01:06 GMT" ;
#                    :uuid = "7ea4f92b-4cdc-42ad-ae1d-938a1866ed04" ;
#                    :ibegin = 0 ;
#                    :ni = 8354 ;
#                    :jbegin = 0 ;
#                    :nj = 9 ;
#                    :file_name = "eNATL60-BLBT02X_1h_20100610_20100704_gridT_20100619-20100619.nc" ;
#                    :TimeStamp = "12/04/2019 22:14:42 +0200" ;
#                    :history = "Wed Apr 17 16:36:16 2019: ncks -F -O -d x,5000,5610 -d y,1600,2500 -d deptht,1,158 -v votemper 01134001-01166400/eNATL60-BLBT02X_1h_20100610_20100624_gridT_20100619-20100619.nc -o /gpfs/projects/pr1egh00/pr1egh01/NOBACKUP/eNATL60/eNATL60-BLBT02X-S/ZOOMs/Meddies_votemper_eNATL60-BLBT02X_1h_20100610_20100624_gridT_20100619-20100619_jk158.nc" ;
#                    :NCO = "4.6.7" ;
#    }
#

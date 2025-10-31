import os
import argparse
import time
import numpy as np
import netCDF4 as nc
import csv
from datetime import datetime, timedelta
from datetime import datetime as dt
import sys
from scipy.interpolate import RegularGridInterpolator
import re
import ComputeMeshToMeshInterpWeights as mshint
#../../WW3/WW3/regtests/ww3_tp2.17/input/level.nc

#40767  0.00000  85.00000  3574.00830
#$EndNodes
#$Elements
#77699
#1  2  3  0  1  0  7625  7626  7527
#2  2  3  0  2  0  8121  8024  8023
#3  2  3  0  3  0  22360  22361  22210
def datenum(d):
    return 366 + d.toordinal() + (d - dt.fromordinal(d.toordinal())).total_seconds()/(24*60*60)

def loadWW3MeshCoords(fl):
    f=open(fl, 'r')
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of nodes
    nn = re.findall(r'\d+', header)
    nn=int(nn[0])
    print(nn)
    xi=np.zeros(nn)
    yi=np.zeros(nn)
    k=0
    for i in range(nn):
        A = f.readline()
        #print(A)
        values = A.split(" ")
        #print(values)
        xi[k]=values[2]
        yi[k]=values[4]
#        print(str(xi[k])+" "+str(yi[k]))
        k=k+1
    header = f.readline() 
    header = f.readline() 
    header = f.readline() # number of elements
    ne=int(header)
    print("ne=",str(ne))
    ei=np.zeros((ne,3), dtype=int)
    print(ei)
    k=0
    for i in range(ne):
        A = f.readline()
        #print(A)
        values = A.split(" ")
        ei[k,0]=int(values[12])
        ei[k,1]=int(values[14])
        ei[k,2]=int(values[16])
        k=k+1

    return xi, yi, ei


N=int(sys.argv[1])
mesh=sys.argv[2]
flout=sys.argv[3]
flin=sys.argv[4]
print("Combining "+str(N)+"interpolations to mesh: "+mesh+"from stofs file:",flin )
print("Output will be written to :"+flout)

#flout="etest.nc"
k=0    
fl=mesh+".files/zeta."+mesh+"."+str(k)+".txt"
data=np.loadtxt(fl)

for k in range(N):
    if k >0:
        fl=mesh+".files/zeta."+mesh+"."+str(k)+".txt"
        datal=np.loadtxt(fl)
        data=np.hstack((data,datal))
        print("reading data part: "+str(k)+" of "+str(N))
        max_index_flat = np.argmax( np.abs(datal) )
        rows, cols = datal.shape
        row = max_index_flat // cols
        col = max_index_flat % cols
        print("maxval="+ str(np.max(np.abs(datal))) + " at "+str(row)+" "+str(col))

# filter bad values - this should be done better

print("-------------------------------------------------")
print("data max vals:")
print(np.max(data))
print("-------------------------------------------------")
print("data min vals:")
print(np.min(data))
print("-------------------------------------------------")

data[np.where(np.abs(data[:])>10. )]=0.
print("-------------------------------------------------")
print("data max vals:")
print(np.max(data))
print("-------------------------------------------------")
print("data min vals:")
print(np.min(data))
print("-------------------------------------------------")

nt=data.shape[0]
nn=data.shape[1]
print("nn="+str(nn))
print("nt="+str(nt))
#datat=np.transpose(data)
flmsh="meshes/"+mesh+".msh"
xi,yi,ei=loadWW3MeshCoords(flmsh)
ne=ei.shape[0]
nn=len(xi)
print('nn ='+str(nn))
print('ne ='+str(ne))

#flin='20241201/stofs_2d_glo.t00z.fields.htp.nc'
data0 = nc.Dataset(flin,"r")
time=np.asarray(data0["time"][:])
nt=len(time)
timevar=data0.variables["time"]
base_date_flin=timevar.base_date
year=base_date_flin[0:4]
month=base_date_flin[5:7]
day=base_date_flin[8:10]
hr=base_date_flin[11:13]
print(year)
print(month)
print(day)
print(hr)
print(time)

dtime0 = datetime(int(year),int(month),int(day),int(hr), 0, 0)
dtime1 = datetime(1990,1,1,0,0,0)
"units         = 'seconds since 2024-04-04 12:00:00        !"
TimeDays=time[:]/86400.
dtime=dtime0.toordinal()+TimeDays - dtime1.toordinal()
print(base_date_flin)

with nc.Dataset(flout, 'w', format='NETCDF4') as nc_file:
    nc_file.title = 'Water level interpolated from STOFS'
    nc_file.description = 'Water level interpolated from STOFS'
    nc_file.source = 'Python netCDF4 library'

    #Global attributes
    nc_file.northernmost_latitude           = str(max(yi))
    nc_file.southernmost_latitude           = str(min(yi))
    nc_file.westernmost_longitude           = str(min(xi))
    nc_file.easternmost_longitude           = str(max(xi))

    # Create dimensions
#    nc_file.createDimension('time', None)  # Unlimited dimension
    nc_file.createDimension('time', nt)  # Unlimited dimension
    nc_file.createDimension('node', nn)
    nc_file.createDimension('element', ne)
    nc_file.createDimension('noel', 3)
    nc_file.createDimension('level', 1)

    # Create variables
    time_var = nc_file.createVariable('time', 'f8', ('time',))
    time_var.units = 'days since 1990-01-01 00:00:00'
    time_var.long_name     = 'julian day (UT)'
    time_var.standard_name = 'time'
    time_var.conventions   = 'relative julian days with decimal part (as parts of the day )'
    time_var.axis='T'
    time_var[:]=np.array(dtime[:])

    lon_var=nc_file.createVariable('longitude', 'f4', ('node',))
    lon_var.units         = 'degree_east'
    lon_var.long_name     = 'longitude'
    lon_var.standard_name = 'longitude'
    lon_var.valid_min     = -180
    lon_var.valid_max     = 360
    lon_var.axis          = 'X'

    lat_var=nc_file.createVariable('latitude', 'f4', ('node',))
    lat_var.units         = 'degree_north'
    lat_var.long_name     = 'latitude'
    lat_var.standard_name = 'latitude'
    lat_var.valid_min     = -0
    lat_var.valid_max     = 180
    lat_var.axis          = 'Y'

    wlv_var=nc_file.createVariable('wlv', 'f4', ('time','node'),fill_value=9.969209968386869e+36) 
    wlv_var.long_name     = 'sea surface height above sea level'
    wlv_var.standard_name = 'sea_surface_height_above_sea_level'
    wlv_var.globwave_name = 'sea_surface_height_above_sea_level'
    wlv_var.units         = 'm'
    #wlv_var._FillValue    = 9.969209968386869e+36
    wlv_var.scale_factor  = 1
    wlv_var.add_offset    = 0
    wlv_var.valid_min     = 0
    wlv_var.valid_max     = 10000

    MAPSTA_var=nc_file.createVariable('MAPSTA', 'i4', ('node',))
    MAPSTA_var.long_name     = 'status map'
    MAPSTA_var.standard_name = 'status map'
    MAPSTA_var.units         = '1'
    MAPSTA_var.valid_min     = -32
    MAPSTA_var.valid_max     = 32

    tri_var=  nc_file.createVariable('tri', 'i4', ('noel','element'))
    tri_var[:]=np.transpose(ei)    

    lon_var[:]=np.array(xi[:])
    lat_var[:]=np.array(yi[:])
    print("data.shape")
    print(data.shape)
#    wlv_var[:]=np.transpose(data) # should be [nn , nt]
    wlv_var[:]=np.array(data)
    print("water level, wlv")
    print(np.max(wlv_var))
    nc_file.close


#lat_var = nc_file.createVariable('latitude', 'f4', ('lat',))
#lon_var = nc_file.createVariable('longitude', 'f4', ('lon',))
#data_var = nc_file.createVariable('wlv', 'f4', ('node','time'))

    
    
    """
    
    
    
    filename = 'example_with_attributes.nc'

# Open a new NetCDF file in write mode
with netCDF4.Dataset(filename, 'w', format='NETCDF4') as nc_file:
    # Add global attributes
    nc_file.title = 'Example NetCDF File with Attributes'
    nc_file.description = 'This file demonstrates how to add global and variable attributes.'
    nc_file.source = 'Python netCDF4 library'

    # Create dimensions
    nc_file.createDimension('time', None)  # Unlimited dimension
    nc_file.createDimension('lat', 10)
    nc_file.createDimension('lon', 20)

    # Create variables
    time_var = nc_file.createVariable('time', 'f8', ('time',))
    lat_var = nc_file.createVariable('lat', 'f4', ('lat',))
    lon_var = nc_file.createVariable('lon', 'f4', ('lon',))
    data_var = nc_file.createVariable('temperature', 'f4', ('time', 'lat', 'lon'))

    # Add attributes to variables
    time_var.units = 'days since 2000-01-01 00:00:00'
    time_var.long_name = 'Time'

    lat_var.units = 'degrees_north'
    lat_var.long_name = 'Latitude'
    lat_var.valid_min = -90.0
    lat_var.valid_max = 90.0

    lon_var.units = 'degrees_east'
    lon_var.long_name = 'Longitude'
    lon_var.valid_min = -180.0
    lon_var.valid_max = 180.0

    data_var.units = 'Kelvin'
    data_var.long_name = 'Air Temperature'
    data_var.standard_name = 'air_temperature'

    # Write some dummy data (optional, but demonstrates a complete file)
    time_var[:] = np.arange(5)
    lat_var[:] = np.linspace(-45, 45, 10)
    lon_var[:] = np.linspace(-90, 90, 20)
    data_var[:] = np.random.rand(5, 10, 20) * 30 + 273.15 # Random temperatures in Kelvin

print(f"NetCDF file '{filename}' created successfully with attributes.")
    
    
    
    
    
Source:
           /mnt/sda/keston/ElevData/aws/level.nc
Format:
           netcdf4
Global Attributes:
           WAVEWATCH_III_version_number    = '6.03'
           WAVEWATCH_III_switches          = 'F90 NOGRB TRKNC NOPA LRB4 DIST MPI SCRIP MLIM PR3 UQ NC4 FLX0 LN1 ST4 STAB0 NL1 BT4 DB1 TR0 BS0 IS0 IC0 REF0 XX0 WCUR WNT2 WNX1 RWND CRT1 CRX1 O0 O1 O2 O2a O2b O2c O3 O4 O5 O6 O7'
           SIN4 namelist parameter BETAMAX = 1.33
           product_name                    = 'ww3.201512_wlv.nc'
           area                            = 'Inlet'
           latitude_resolution             = 'n/a'
           longitude_resolution            = 'n/a'
           southernmost_latitude           = '40.38446'
           northernmost_latitude           = '40.99023'
           westernmost_longitude           = '-72.92410'
           easternmost_longitude           = '-72.03251'
           minimum_altitude                = '-12000 m'
           maximum_altitude                = '9000 m'
           altitude_resolution             = 'n/a'
           start_date                      = '2015-12-14 00:00:00'
           stop_date                       = '2015-12-15 00:00:00'
Dimensions:
           level   = 1
           node    = 3070
           element = 5780
           time    = 25    (UNLIMITED)
           noel    = 3
Variables:
    longitude
           Size:       3070x1
           Dimensions: node
           Datatype:   single
           Attributes:
                       units         = 'degree_east'
                       long_name     = 'longitude'
                       standard_name = 'longitude'
                       valid_min     = -180
                       valid_max     = 360
                       axis          = 'X'
    latitude 
           Size:       3070x1
           Dimensions: node
           Datatype:   single
           Attributes:
                       units         = 'degree_north'
                       long_name     = 'latitude'
                       standard_name = 'latitude'
                       valid_min     = -90
                       valid_max     = 180
                       axis          = 'Y'
    time     
           Size:       25x1
           Dimensions: time
           Datatype:   double
           Attributes:
                       long_name     = 'julian day (UT)'
                       standard_name = 'time'
                       units         = 'days since 1990-01-01 00:00:00'
                       conventions   = 'relative julian days with decimal part (as parts of the day )'
                       axis          = 'T'
    tri      
           Size:       3x5780
           Dimensions: noel,element
           Datatype:   int32
    MAPSTA   
           Size:       3070x1
           Dimensions: node
           Datatype:   int16
           Attributes:
                       long_name     = 'status map'
                       standard_name = 'status map'
                       units         = '1'
                       valid_min     = -32
                       valid_max     = 32
    wlv      
           Size:       3070x25
           Dimensions: node,time
           Datatype:   single
           Attributes:
                       long_name     = 'sea surface height above sea level'
                       standard_name = 'sea_surface_height_above_sea_level'
                       globwave_name = 'sea_surface_height_above_sea_level'
                       units         = 'm'
                       _FillValue    = 9.969209968386869e+36
                       scale_factor  = 1
                       add_offset    = 0
                       valid_min     = 0
                       valid_max     = 10000
"""

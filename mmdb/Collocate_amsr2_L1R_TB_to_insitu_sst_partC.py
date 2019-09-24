#!/usr/bin/env python
# coding: utf-8

# # This is the in situ and SSS collocation code. 
# # this is the part A of the program that searches for L1R files that have any data where cruise is


import sys
import numpy as np
import xarray as xr
from glob import glob
from pyresample import image, geometry, load_area, save_quicklook, SwathDefinition, area_def2basemap
from pyresample.kd_tree import resample_nearest
import gzip
import shutil
from scipy import spatial
sys.path.append('./subroutines/')
from read_routines import read_usv

# # Define a function to read in insitu data
# - Read in the Saildrone USV file either from a local disc or using OpenDAP.
# - add room to write collocated data to in situ dataset
# input **********************************

# ## First let's figure out what orbital files actually have data in our area of interest.  To do this, use the pyresample software
#
# - read in the in situ data
# - calculate the in situ min/max dates to know what files to check
#
# Now we have our time of interest
#
# - loop through the satellite data
# - calculate the in situ min/max lat/lon on the same day to define a small box of interest
# - use pyresample to map the data onto a predefined 0.1 deg resolution spatial grid
# - subset the gridded map to the area of interest
# - see if there is any valid data in that area
# - if there is any valid data, save the filename into a list
#
#
input_iusv_start = int(input("Enter start cruise processing number 0-10: "))
input_iusv_end = int(input("Enter stop cruise processing number 0-10: "))
adir_usv = str(input("Enter directory for USV data: "))

#intialize grid
for iusv in range(input_iusv_start,input_iusv_end):
    ds_usv, usv_name = read_usv(adir_usv,iusv)
    print('usv', usv_name)
    fileout = adir_usv + 'mmdb_collocation_test/' + usv_name + 'AMSR2MMDB_usv2_testing.nc'
    ds_usv = xr.open_dataset(fileout)
# now collocation with orbital data is finished.  re-open file and create the mean values for each matchup so there aren't repeates

    fileout_norepeat = fileout[:-3]+'_norepeats_testing.nc'
 #   ds_usv = ds_usv.where(ds_usv.tb<10000,np.nan)
    ds_usv = ds_usv.where(ds_usv.amsr2_scan>10000,np.nan)  #keep valid scan nmbers, replace rest with nan
    ilen,index = ds_usv['insitu.time'].size,0
    ds_tem = ds_usv.copy(deep=True)
    dsst,dair,dwav,dpres,dsstu,dvwnd,duwnd,dsal,dchl,dgwnd,dlat,dlon,dut = [],[],[],[],[],[],[],[],[],[],[],[],np.empty((),dtype='datetime64')
    index=0
    while index <= ilen-2:
        index += 1
#        if np.isnan(ds_usv.tb[index]):
#            continue
        if np.isnan(ds_usv.amsr2_scan[index]):
            continue
        result = np.where((ds_usv.amsr2_scan == ds_tem.amsr2_cell[index].data) & (ds_usv.amsr2_scan == ds_tem.amsr2_cell[index].data))
        #duu=np.append(duu,ds_usv.smap_SSS[result[0][0]].data)
        #duu2=np.append(duu2,ds_usv.smap_iqc_flag[result[0][0]].data)
        dsst=np.append(dsst,ds_usv['insitu.sea_surface_temperature'][result].mean().data)
        dair=np.append(dair,ds_usv['insitu.air_temperature'][result].mean().data)
        dwav=np.append(dwav,ds_usv['insitu.sig_wave_height'][result].mean().data)
        dpres=np.append(dpres,ds_usv['insitu.baro_pres'][result].mean().data)
        dsstu=np.append(dsstu,ds_usv['insitu.sst_uncertainty'][result].mean().data)
        dvwnd=np.append(dvwnd,ds_usv['insitu.vwnd'][result].mean().data)
        duwnd=np.append(duwnd,ds_usv['insitu.uwnd'][result].mean().data)
        dsal=np.append(dsal,ds_usv['insitu.salinity'][result].mean().data)
        dchl=np.append(dchl,ds_usv['insitu.chlor'][result].mean().data)
        dgwnd=np.append(dgwnd,ds_usv['insitu.gust_wind'][result].mean().data)
        dlat=np.append(dlat,ds_usv['insitu.lat'][result].mean().data)
        dlon=np.append(dlon,ds_usv['insitu.lon'][result].mean().data)
        dut=np.append(dut,ds_usv['insitu.time'][result].mean().data)
        #ds_usv.tb[result]=np.nan
    dut2 = dut[1:]  #remove first data point which is a repeat from what array defined
    ds_new=xr.Dataset(data_vars={'insitu.sea_surface_temperature': ('time',dsst),
                                 'insitu.air_temperature': ('time',dair),
                                 'insitu.sig_wave_height':('time',dwav),
                                 'insitu.baro_pres':('time',dpres),
                                 'insitu.sst_uncertainty':('time',dsstu),
                                 'insitu.vwnd':('time',dvwnd),
                                 'insitu.uwnd':('time',duwnd),
                                 'insitu.salinity': ('time', dsal),
                                 'insitu.chlor': ('time', dchl),
                                 'insitu.gust_wind': ('time', dgwnd),
                                 'lon': ('time',dlon),
                                 'lat': ('time',dlat)},
                      coords={'time':dut2})
    ds_new.to_netcdf(fileout_norepeat)

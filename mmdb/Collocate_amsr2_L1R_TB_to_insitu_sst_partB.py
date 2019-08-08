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
adir_l1r = str(input("Enter directory for L1R data: "))

def read_usv(adir_usv,iusv):
    filename_usv_list = ['pmel_2015_sd126-ALL-1_min-v1.nc',
                         'pmel_2015_sd128-ALL-1_min-v1.nc',
                         'pmel_2016_sd126-ALL-1_min-v1.nc',
                         'pmel_2016_sd128-ALL-1_min-v1.nc',
                         'arctic_2019_sd1033-ALL-1_min-v1.nc',
                         'arctic_2019_sd1034-ALL-1_min-v1.nc',
                         'arctic_2019_sd1035-ALL-1_min-v1.nc',
                         'arctic_2019_sd1036-ALL-1_min-v1.nc',
                         'arctic_2019_sd1037-ALL-1_min-v1.nc',
                         'saildrone-gen_5-antarctica_circumnavigation_2019-sd1020-20190119T040000-20190803T043000-1440_minutes-v1.1564857794963.nc']
    name_usv_list = ['arctic2015_126',
                     'arctic2015_128',
                     'arctic2016_126',
                     'arctic2016_128',
                     'arctic2019_1033',
                     'arctic2019_1034',
                     'arctic2019_1035',
                     'arctic2019_1036',
                     'arctic2019_1037',
                    'antarctic2019']

    filename_usv = adir_usv + filename_usv_list[iusv]
    print('FILEIN:',filename_usv)
    ds_usv = xr.open_dataset(filename_usv)
    ds_usv.close()
#NEED TO FIND OUT IF wind_speed is to/from wind_direction ?
    if (iusv==0 or iusv==1):  #1033
        ds_usv = ds_usv.rename({'temp_air_mean':'TEMP_AIR_MEAN','rh_mean':'RH_MEAN','baro_pres_mean':'BARO_PRES_MEAN',
                                'sal_mean':'SAL_MEAN','temp_ctd_mean':'TEMP_CTD_MEAN','temp_o2_mean':'TEMP_O2_MEAN',
                                'chlor_mean':'CHLOR_MEAN','gust_wnd_mean':'GUST_WND_MEAN','temp_ctd_stddev':'TEMP_CTD_STDDEV'})
        tem_att=ds_usv.wind_speed_mean.attrs
        ds_usv['wind_speed_mean']=ds_usv.wind_speed_mean*.51444
        ds_usv.wind_speed_mean.attrs=tem_att
        ds_usv.wind_speed_mean.attrs['units']='m s-1'
        uwnd = ds_usv.wind_speed_mean*np.cos(np.deg2rad(ds_usv.wind_direction_mean))
        vwnd = ds_usv.wind_speed_mean*np.sin(np.deg2rad(ds_usv.wind_direction_mean))
        ds_usv['UWND_MEAN']=uwnd
        ds_usv.UWND_MEAN.attrs={'standard_name':'eastward_wind','long_name':'Eastward wind speed','units':'m s-1','installed_height':'5.2'}
        ds_usv['VWND_MEAN']=vwnd
        ds_usv.VWND_MEAN.attrs={'standard_name':'northward_wind','long_name':'Northward wind speed','units':'m s-1','installed_height':'5.2'}
        ilen = ds_usv.time.shape[0]
        ds_usv['WWND_MEAN'] = xr.DataArray(np.ones(ilen) * np.nan, coords={'time': ds_usv.time}, dims=('time'))
        ds_usv.WWND_MEAN.attrs={'standard_name':'upward_wind_velocity','long_name':'upward wind speed','units':'m s-1','installed_height':'5.2'}
    if (iusv==2 or iusv==3):  #1033
        ds_usv = ds_usv.rename({'temp_air_mean':'TEMP_AIR_MEAN','rh_mean':'RH_MEAN','baro_pres_mean':'BARO_PRES_MEAN',
                                'sal_mean':'SAL_MEAN','temp_ctd_mean':'TEMP_CTD_MEAN','temp_o2_mean':'TEMP_O2_MEAN',
                                'chlor_mean':'CHLOR_MEAN','gust_wnd_mean':'GUST_WND_MEAN','temp_ctd_stddev':'TEMP_CTD_STDDEV'})
        tem_att=ds_usv.wind_speed.attrs
        ds_usv['wind_speed']=ds_usv.wind_speed*.51444
        ds_usv.wind_speed.attrs=tem_att
        ds_usv.wind_speed.attrs['units']='m s-1'
        uwnd = ds_usv.wind_speed*np.cos(np.deg2rad(ds_usv.wind_direction))
        vwnd = ds_usv.wind_speed*np.sin(np.deg2rad(ds_usv.wind_direction))
        ds_usv['UWND_MEAN']=uwnd
        ds_usv.UWND_MEAN.attrs={'standard_name':'eastward_wind','long_name':'Eastward wind speed','units':'m s-1','installed_height':'5.2'}
        ds_usv['VWND_MEAN']=vwnd
        ds_usv.VWND_MEAN.attrs={'standard_name':'northward_wind','long_name':'Northward wind speed','units':'m s-1','installed_height':'5.2'}
        ilen = ds_usv.time.shape[0]
        ds_usv['WWND_MEAN']=xr.DataArray(np.ones(ilen)*np.nan,coords={'time':ds_usv.time},dims=('time'))
        ds_usv.WWND_MEAN.attrs={'standard_name':'upward_wind_velocity','long_name':'upward wind speed','units':'m s-1','installed_height':'5.2'}
    if iusv==4:  #1033
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_CTD_RBR_STDDEV':'TEMP_CTD_STDDEV','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN'})
    if iusv==5:  #1034
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_CTD_RBR_STDDEV':'TEMP_CTD_STDDEV','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN'})
    if iusv==6:  #1035
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_CTD_RBR_STDDEV':'TEMP_CTD_STDDEV','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN','WIND_MEASUREMENT_MEAN_HEIGHT':'WIND_MEAN_HEIGHT'})
    if iusv==7:  #1036
        ds_usv = ds_usv.isel(time=slice(100,-1))                                                                   #        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN'})
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_CTD_RBR_STDDEV':'TEMP_CTD_STDDEV','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN'})
    if iusv==8:  #1037
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_CTD_RBR_STDDEV':'TEMP_CTD_STDDEV','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN'})
    if iusv==9:  #1037
        ds_usv = ds_usv.isel(trajectory=0).swap_dims({'obs':'time'}).rename({'latitude':'lat','longitude':'lon','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN'})#TEMP_CTD_RBR_MEAN':'TEMP_
    if (iusv==9 or iusv<=3):
        ilen = ds_usv.time.shape[0]
        ds_usv['WIND_HEIGHT_MEAN']=xr.DataArray(np.ones(ilen)*np.nan,coords={'time':ds_usv.time},dims=('time'))
        ds_usv.WIND_HEIGHT_MEAN.attrs={'long_name':'Wind measurement height','units':'m','installed_height':'5.2'}
        ds_usv['WAVE_DOMINANT_PERIOD']=xr.DataArray(np.ones(ilen)*np.nan,coords={'time':ds_usv.time},dims=('time'))
        ds_usv.WAVE_DOMINANT_PERIOD.attrs={'standard_name':'sea_surface_wave_period_at_variance_spectral_density_maximum','long_name':'Dominant wave period','units':'s','installed_height':'0.34'}
        ds_usv['WAVE_SIGNIFICANT_HEIGHT']=xr.DataArray(np.ones(ilen)*np.nan,coords={'time':ds_usv.time},dims=('time'))
        ds_usv.WAVE_SIGNIFICANT_HEIGHT.attrs={'standard_name':'sea_surface_wave_significant_height','long_name':'Significant wave height','units':'m','installed_height':'0.34'}

    #add room to write collocated data information
    ilen = ds_usv.time.shape[0]
    ds_usv['insitu.dtime']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_name']=xr.DataArray(np.empty(ilen,dtype=str),coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_dist']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_scan']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_cell']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['insitu.id']=xr.DataArray(np.empty(ilen,dtype=str),coords={'time':ds_usv.time},dims=('time'))

    return ds_usv,name_usv_list[iusv]

#intialize grid
for num_usv in range(input_iusv_start,input_iusv_end):
    ds_usv, usv_name = read_usv(adir_usv,num_usv)
    filelist = adir_usv + usv_name + 'AMSR2MMDB_filesave2.nc'
    fileout = adir_usv + usv_name + 'AMSR2MMDB_usv2.nc'
    df = xr.open_dataset(filelist)
    for file2 in df.filenames.data:
        file = file2
        file.replace('\\', '/')
#replace drive
        ipos = file.find('amsr2')
        file = adir_l1r + file[ipos:]
        print(file[ipos + 1:])
        print('opening:',file)
        temp_file = 'c:/temp/tem_' + str(num_usv) + '.h5'
        if ds_usv.time.min().dt.year.data < 2018:  # early files gzipped
            with gzip.open(file, 'rb') as f_in:
                with open(temp_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            ds_l1r = xr.open_dataset(temp_file)
        else:
            ds_l1r = xr.open_dataset(file)
        ds_l1r.close()
        xlat = ds_l1r['Latitude of Observation Point for 89A'][:,::2]
        xlon = ds_l1r['Longitude of Observation Point for 89A'][:,::2]
        tb = ds_l1r['Brightness Temperature (res06,10.7GHz,H)']
        ph0 = ds_l1r['Brightness Temperature (res06,10.7GHz,H)'].phony_dim_0
        ph1 = ds_l1r['Brightness Temperature (res06,10.7GHz,H)'].phony_dim_1
        tem_time = np.datetime64('1993-01-01') + (ds_l1r['Scan Time'].data * 1000).astype('timedelta64[ms]')
        ds = xr.Dataset({'time': (['phony_dim_0'], tem_time),
                         'tb': (['phony_dim_0', 'phony_dim_1'], tb),
                         'lat': (['phony_dim_0', 'phony_dim_1'], xlat.data),
                         'lon': (['phony_dim_0', 'phony_dim_1'], xlon.data)},
                         coords={'phony_dim_0': (['phony_dim_0'], ph0),
                                 'phony_dim_1': (['phony_dim_1'], ph1)})

        ds2 = ds.stack(z=('phony_dim_0', 'phony_dim_1')).reset_index('z')
        # drop nan
        ds_drop = ds2.where(np.isfinite(ds2.lon), drop=True)
        lats = ds_drop.lat.data
        lons = ds_drop.lon.data
        inputdata = list(zip(lons.ravel(), lats.ravel()))
        tree = spatial.KDTree(inputdata)
        orbit_time = ds.time.max().data - np.timedelta64(1, 'D')
        orbit_time2 = ds.time.max().data + np.timedelta64(1, 'D')
        usv_subset = ds_usv.sel(time=slice(orbit_time, orbit_time2))
        ilen = ds_usv.time.size
        for iusv in range(ilen):
            if (ds_usv.time[iusv] < orbit_time) or (ds_usv.time[iusv] > orbit_time2):
                continue
            pts = np.array([ds_usv.lon[iusv], ds_usv.lat[iusv]])
            #        pts = np.array([ds_usv.lon[iusv]+360, ds_usv.lat[iusv]])
            tree.query(pts, k=1)
            i = tree.query(pts)[1]
            rdist = tree.query(pts)[0]
            # don't use matchups more than 25 km away
            if rdist > .25:
                continue
            # use .where to find the original indices of the matched data point
            # find by matching sss and lat, just randomly chosen variables, you could use any
            result = np.where((ds.tb == ds_drop.tb[i].data) & (ds.lat == ds_drop.lat[i].data))
            listOfCoordinates = list(zip(result[0], result[1]))
            if len(listOfCoordinates) == 0:
                continue
            ii, jj = listOfCoordinates[0][0], listOfCoordinates[0][1]
            deltaTa = ((ds_usv.time[iusv] - ds.time[ii]).data) / np.timedelta64(1, 'm')
            if np.abs(deltaTa) < np.abs(ds_usv['insitu.dtime'][iusv].data):
                ds_usv['insitu.dtime'][iusv] = deltaTa
                ds_usv.amsr2_name[iusv] = file2
                ds_usv.amsr2_dist[iusv] = rdist
                ds_usv.amsr2_scan[iusv] = ii
                ds_usv.amsr2_cell[iusv] = jj

    ds_usv = ds_usv.rename({'TEMP_CTD_MEAN':'insitu.sea_surface_temperature','TEMP_CTD_STDDEV':'insitu.sst_uncertainty',
                             'TEMP_AIR_MEAN':'insitu.air_temperature','VWND_MEAN':'insitu.vwnd','UWND_MEAN':'insitu.uwnd',
                             'WAVE_SIGNIFICANT_HEIGHT':'insitu.sig_wave_height','SAL_MEAN':'insitu.salinity','CHLOR_MEAN':'insitu.chlor',
                             'BARO_PRES_MEAN':'insitu.baro_pres','RH_MEAN':'insitu.rel_humidity','GUST_WND_MEAN':'insitu.gust_wind',
                             'lat':'insitu.lat','lon':'insitu.lon','time':'insitu.time'})

    ds_usv.to_netcdf(fileout)

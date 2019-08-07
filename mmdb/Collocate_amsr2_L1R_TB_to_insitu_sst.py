#!/usr/bin/env python
# coding: utf-8

# # This is the in situ and SSS collocation code. 
# 

# In[1]:

#import os
import sys
import numpy as np
#import matplotlib.pyplot as plt
#import datetime as dt
#import pandas as pd
import xarray as xr
#import scipy
from glob import glob
#import cartopy.crs as ccrs
#from pyresample.geometry import AreaDefinition
from pyresample import image, geometry, load_area, save_quicklook, SwathDefinition, area_def2basemap
from pyresample.kd_tree import resample_nearest
#from math import radians, cos, sin, asin, sqrt
#from scipy import spatial
#import os.path
#from os import path
import gzip
import shutil


# # Define a function to read in insitu data
# - Read in the Saildrone USV file either from a local disc or using OpenDAP.
# - add room to write collocated data to in situ dataset
# 

# In[5]:
input_iusv_start=int(str(sys.argv[1]))
input_iusv_end=int(str(sys.argv[2]))


def read_usv(iusv):
#    filename_usv_list = ['F:/data/cruise_data/saildrone/noaa_arctic/PMEL_2015/126/pmel_2015_sd126-ALL-1_min-v1.nc',
#                         'F:/data/cruise_data/saildrone/noaa_arctic/PMEL_2015/128/pmel_2015_sd128-ALL-1_min-v1.nc',
#                         'F:/data/cruise_data/saildrone/noaa_arctic/PMEL_2016/126/pmel_2016_sd126-ALL-1_min-v1.nc',
#                         'F:/data/cruise_data/saildrone/noaa_arctic/PMEL_2016/128/pmel_2016_sd128-ALL-1_min-v1.nc',
#                         'F:/data/cruise_data/saildrone/2019_arctic/daily_files/arctic_2019_sd1033-ALL-1_min-v1.nc',
#                         'F:/data/cruise_data/saildrone/2019_arctic/daily_files/arctic_2019_sd1034-ALL-1_min-v1.nc',
#                         'F:/data/cruise_data/saildrone/2019_arctic/daily_files/arctic_2019_sd1035-ALL-1_min-v1.nc',
#                         'F:/data/cruise_data/saildrone/2019_arctic/daily_files/arctic_2019_sd1036-ALL-1_min-v1.nc',
#                         'F:/data/cruise_data/saildrone/2019_arctic/daily_files/arctic_2019_sd1037-ALL-1_min-v1.nc',
#                        'F:/data/cruise_data/saildrone/antarctic/saildrone-gen_5-antarctica_circumnavigation_2019-sd1020-20190119T040000-20190803T043000-1440_minutes-v1.1564857794963.nc']
    filename_usv_list = ['C:/Users/gentemann/Google Drive/private/tem_saildrone/pmel_2015_sd126-ALL-1_min-v1.nc',
                         'C:/Users/gentemann/Google Drive/private/tem_saildrone/pmel_2015_sd128-ALL-1_min-v1.nc',
                         'C:/Users/gentemann/Google Drive/private/tem_saildrone/pmel_2016_sd126-ALL-1_min-v1.nc',
                         'C:/Users/gentemann/Google Drive/private/tem_saildrone/pmel_2016_sd128-ALL-1_min-v1.nc',
                         'C:/Users/gentemann/Google Drive/private/tem_saildrone/arctic_2019_sd1033-ALL-1_min-v1.nc',
                         'C:/Users/gentemann/Google Drive/private/tem_saildrone/arctic_2019_sd1034-ALL-1_min-v1.nc',
                         'C:/Users/gentemann/Google Drive/private/tem_saildrone/arctic_2019_sd1035-ALL-1_min-v1.nc',
                         'C:/Users/gentemann/Google Drive/private/tem_saildrone/arctic_2019_sd1036-ALL-1_min-v1.nc',
                         'C:/Users/gentemann/Google Drive/private/tem_saildrone/arctic_2019_sd1037-ALL-1_min-v1.nc',
                         'C:/Users/gentemann/Google Drive/private/tem_saildrone/saildrone-gen_5-antarctica_circumnavigation_2019-sd1020-20190119T040000-20190803T043000-1440_minutes-v1.1564857794963.nc']
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

    filename_usv = filename_usv_list[iusv]
#    if iusv==3:
#    elif iusv<3:
#    elif (iusv>3) & (iusv<8):
    ds_usv = xr.open_dataset(filename_usv)
    ds_usv.close()
#NEED TO FIND OUT IF wind_speed is to/from wind_direction ?
    if (iusv==0 or iusv==1):  #1033
        ds_usv = ds_usv.rename({'temp_air_mean':'TEMP_AIR_MEAN','rh_mean':'RH_MEAN','baro_pres_mean':'BARO_PRES_MEAN',
                                'sal_mean':'SAL_MEAN','temp_ctd_mean':'TEMP_CTD_MEAN','temp_o2_mean':'TEMP_O2_MEAN',
                                'chlor_mean':'CHLOR_MEAN','gust_wnd_mean':'GUST_WND_MEAN'})
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
        ds_usv['WWND_MEAN']=xr.DataArray(np.ones(ilen)*np.nan,coords={'time':ds_usv.time},dims=('time'))
        ds_usv.WWND_MEAN.attrs={'standard_name':'upward_wind_velocity','long_name':'upward wind speed','units':'m s-1','installed_height':'5.2'}
    if (iusv==2 or iusv==3):  #1033
        ds_usv = ds_usv.rename({'temp_air_mean':'TEMP_AIR_MEAN','rh_mean':'RH_MEAN','baro_pres_mean':'BARO_PRES_MEAN',
                                'sal_mean':'SAL_MEAN','temp_ctd_mean':'TEMP_CTD_MEAN','temp_o2_mean':'TEMP_O2_MEAN',
                                'chlor_mean':'CHLOR_MEAN','gust_wnd_mean':'GUST_WND_MEAN'})
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
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN'})
    if iusv==5:  #1034
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN'})
    if iusv==6:  #1035
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN'})
    if iusv==7:  #1036
        ds_usv = ds_usv.isel(time=slice(100,-1))                                                                   #        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN'})
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN'})
    if iusv==8:  #1037
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN'})
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
    ds_usv['deltaT']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_name']=xr.DataArray(np.empty(ilen,dtype=str),coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_dist']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_ydim']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_xdim']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))

    return ds_usv,name_usv_list[iusv]


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

# In[37]:


#intialize grid
for iusv in range(input_iusv_start,input_iusv_end):
    area_def = load_area('areas.cfg', 'pc_world')
    rlon=np.arange(-180,180,.1)
    rlat=np.arange(90,-90,-.1)

    ds_usv,name_usv = read_usv(iusv)

#    adir = 'C:/Users\gentemann/Google Drive/public/temp/'
    adir = 'd:/'
    if ds_usv.time.min().dt.year.data<2018:
        sat_directory=adir+'amsr2/L1r/v2/'
    else:
        sat_directory=adir+'amsr2_update/ftp.gportal.jaxa.jp/standard/GCOM-W/GCOM-W.AMSR2/L1R/2/'
    fileout = 'C:/Users\gentemann/Google Drive/private/tem_saildrone/'+name_usv+'AMSR2MMDB_filesave2.nc'
    file_end = '*.h5.gz'
    
#    if path.exists(fileout):
#        continue
    #init filelist
    file_save=[]

    #search usv data
    minday,maxday = ds_usv.time[0],ds_usv.time[-1]
    usv_day = minday
    print(minday.data,maxday.data)
    while usv_day<=maxday:
        
#        check_day = np.datetime64(str(usv_day.dt.year.data)+'-'+str(usv_day.dt.month.data).zfill(2)+'-'+str(usv_day.dt.day.data).zfill(2))
#        usv_day1 = usv_day + np.timedelta64(1,'D')
#        check_day1 = np.datetime64(str(usv_day1.dt.year.data)+'-'+str(usv_day1.dt.month.data).zfill(2)+'-'+str(usv_day1.dt.day.data).zfill(2))
#        ds_day = ds_usv.sel(time=slice(check_day,check_day1))

#while looping through USV data, look at data +-1 day 
        ds_day = ds_usv.sel(time=slice(usv_day-np.timedelta64(1,'D'),usv_day+np.timedelta64(1,'D')))
        ilen = ds_day.time.size
        if ilen<1:   #don't run on days without any data
            continue
        minlon,maxlon,minlat,maxlat = ds_day.lon.min().data,ds_day.lon.max().data,ds_day.lat.min().data,ds_day.lat.max().data
        #caluclate filelist
        syr=str(usv_day.dt.year.data)
        smon=str(usv_day.dt.month.data).zfill(2)
        sdy=str(usv_day.dt.day.data).zfill(2)
        #the more recent data is in daily directories, so easy to search
        #the older data, pre 2018 is in monthly directories so only search files for day
        if usv_day.dt.year.data<2018:
            adir_list=sat_directory+syr+'/'+smon+'/'+sdy+'/'+file_end
            filelist = glob(adir_list)
        else:
            adir_list=sat_directory+syr+'/'+smon+'/'+'/GW1AM2_'+syr+smon+sdy+file_end
            filelist = glob(adir_list)  
        #print(sat_directory+syr+'/'+smon+'/'+'/GW1AM2_'+syr+smon+sdy+file_end)
        print(usv_day.data,'numfiles:',len(filelist))
        print(adir_list)
        temp_file='c:/temp/tem_'+str(iusv)+'.h5'
        x,y,z = [],[],[]
        for file in filelist:
            file.replace('\\','/')
            with gzip.open(file, 'rb') as f_in:
                with open(temp_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            ds=xr.open_dataset(temp_file)
            ds.close()
            xlat=ds['Latitude of Observation Point for 89A'][:,::2]
            xlon=ds['Longitude of Observation Point for 89A'][:,::2]
            x = xlon.data
            y = xlat.data
            z = ds['Brightness Temperature (res06,6.9GHz,H)'].data*.01 
            lons,lats,data = x,y,z 
            swath_def = SwathDefinition(lons, lats)
            result1 = resample_nearest(swath_def, data, area_def, radius_of_influence=20000, fill_value=None)
            da = xr.DataArray(result1,name='tb6h',coords={'lat':rlat,'lon':rlon},dims=('lat','lon'))
            subset = da.sel(lat = slice(maxlat,minlat),lon=slice(minlon,maxlon))
            num_obs = np.isfinite(subset).sum()
            if num_obs>0:
                file_save = np.append(file_save,file)
        usv_day += np.timedelta64(1,'D')
        df = xr.DataArray(file_save,name='filenames')
        df.to_netcdf(fileout)



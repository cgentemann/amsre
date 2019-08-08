#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import numpy as np

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
    name_usv_list = ['pmel_2015_sd126',
                     'pmel_2015_sd128',
                     'pmel_2016_sd126',
                     'pmel_2016_sd128',
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


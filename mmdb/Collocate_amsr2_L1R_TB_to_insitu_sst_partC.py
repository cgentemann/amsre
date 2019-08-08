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
        ds_usv = ds_usv.rename({'TEMP_CTD_RBR_MEAN':'TEMP_CTD_MEAN','TEMP_O2_RBR_MEAN':'TEMP_O2_MEAN','SAL_RBR_MEAN':'SAL_MEAN','CHLOR_WETLABS_MEAN':'CHLOR_MEAN','WIND_MEASUREMENT_MEAN_HEIGHT':'WIND_MEAN_HEIGHT'})
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
    ds_usv['insitu.dtime']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_name']=xr.DataArray(np.empty(ilen,dtype=str),coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_dist']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_scan']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['amsr2_cell']=xr.DataArray(np.ones(ilen)*999999,coords={'time':ds_usv.time},dims=('time'))
    ds_usv['insitu.id']=xr.DataArray(np.empty(ilen,dtype=str),coords={'time':ds_usv.time},dims=('time'))

    return ds_usv,name_usv_list[iusv]

#intialize grid
for iusv in range(input_iusv_start,input_iusv_end):
    num_usv = 0
    ds_usv, usv_name = read_usv(adir_usv,num_usv)
    fileout = adir_usv + usv_name + 'AMSR2MMDB_usv2.nc'
    ds_usv = xr.open_dataset(fileout)
# now collocation with orbital data is finished.  re-open file and create the mean values for each matchup so there aren't repeates

    fileout_norepeat = fileout[:-3]+'_norepeats.nc'
    ds_usv = ds_usv.where(ds_usv.tb<10000,np.nan)
    ilen,index = ds_usv.dims['time'],0
    ds_tem = ds_usv.copy(deep=True)
    dsst,dair,dwav,dpres,dsstu,dvwnd,duwnd,dsal,dchl,dgwnd,dlat,dlon, dut = [],[],[],[],[],[],[],[],[],[],[],[],np.empty((),dtype='datetime64')
    index=0
    while index <= ilen-2:
        index += 1
        if np.isnan(ds_usv.tb[index]):
            continue
        if np.isnan(ds_usv.amsr2_scan[index]):
            continue
        result = np.where((ds_usv.amsr2_scan == ds_tem.amsr2_cell[index].data) & (ds_usv.amsr2_scan == ds_tem.amsr2_cell[index].data))
        #duu=np.append(duu,ds_usv.smap_SSS[result[0][0]].data)
        #duu2=np.append(duu2,ds_usv.smap_iqc_flag[result[0][0]].data)
        dsst=np.append(duv1,ds_usv['insitu.sea_surface_temperature'][result].mean().data)
        dair=np.append(duv1,ds_usv['insitu.air_temperature'][result].mean().data)
        dwav=np.append(duv1,ds_usv['insitu.sig_wave_height'][result].mean().data)
        dpres=np.append(duv1,ds_usv['insitu.baro_pres'][result].mean().data)
        dsstu=np.append(duv1,ds_usv['insitu.sst_uncertainty'][result].mean().data)
        dvwnd=np.append(duv1,ds_usv['insitu.vwnd'][result].mean().data)
        duwnd=np.append(duv1,ds_usv['insitu.uwnd'][result].mean().data)
        dsal=np.append(duv1,ds_usv['insitu.salinity'][result].mean().data)
        dchl=np.append(duv1,ds_usv['insitu.chlor'][result].mean().data)
        dgwnd=np.append(duv1,ds_usv['insitu.gust_wind'][result].mean().data)
        dlat=np.append(dlat,ds_usv['insitu.lat'][result].mean().data)
        dlon=np.append(dlon,ds_usv['insitu.lon'][result].mean().data)
        dut=np.append(dut,ds_usv['insitu.time'][result].mean().data)
        ds_usv.tb[result]=np.nan
    dut2 = dut[1:]  #remove first data point which is a repeat from what array defined
    ds_new=xr.Dataset(data_vars={'insitu.sea_surface_temperature': ('time',dsst),
                                 'insitu.air_temperature': ('time',ddair),
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

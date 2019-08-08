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
adir_l1r = str(input("Enter directory for L1R data: "))


#intialize grid
for iusv in range(input_iusv_start,input_iusv_end):
    area_def = load_area('areas.cfg', 'pc_world')
    rlon=np.arange(-180,180,.1)
    rlat=np.arange(90,-90,-.1)

    ds_usv,name_usv = read_usv(adir_usv,iusv)

#    adir = 'C:/Users\gentemann/Google Drive/public/temp/'
   # adir = 'd:/'
    if ds_usv.time.min().dt.year.data<2018:
        sat_directory = adir_l1r + 'amsr2/l1r/v2/'
        file_end = '*.h5.gz'
    else:
        sat_directory = adir_l1r + 'amsr2_update/ftp.gportal.jaxa.jp/standard/GCOM-W/GCOM-W.AMSR2/L1R/2/'
        file_end = '*.h5'
    fileout = adir_usv +name_usv+'AMSR2MMDB_filesave2.nc'

#    if path.exists(fileout):
#        continue
    #init filelist
    file_save=[]

    #search usv data
    minday,maxday = ds_usv.time[0],ds_usv.time[-1]
    usv_day = minday
    print(minday.data,maxday.data)
    while usv_day<=maxday:
        
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
            if ds_usv.time.min().dt.year.data < 2018:  #early files gzipped
                with gzip.open(file, 'rb') as f_in:
                    with open(temp_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                ds=xr.open_dataset(temp_file)
            else:
                ds=xr.open_dataset(file)
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



#!/usr/bin/env python
# coding: utf-8

# # This is the in situ and SSS collocation code. 
# # this is the part A of the program that searches for L1R files that have any data where cruise is


import sys
import numpy as np
import xarray as xr
from pyresample import image, geometry, load_area, save_quicklook, SwathDefinition, area_def2basemap
from pyresample.kd_tree import resample_nearest
from scipy import spatial

sys.path.append('./subroutines/')
from read_routines import read_usv, get_filelist_amsr2_l1r,get_orbital_data_amsr2

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
'''some definitions'''
area_def = load_area('areas.cfg', 'pc_world')
rlon = np.arange(-180, 180, .1)
rlat = np.arange(90, -90, -.1)

input_iusv_start = int(input("Enter start cruise processing number 0-10: "))
input_iusv_end = int(input("Enter stop cruise processing number 0-10: "))
adir_usv = str(input("Enter directory for USV data: "))
adir_l1r = str(input("Enter directory for L1R data: "))

#intialize grid
for iusv in range(input_iusv_start,input_iusv_end):

    file_save=[]  #initialize variable for each processing run

    ds_usv,name_usv = read_usv(adir_usv,iusv)
    fileout = adir_usv + 'mmdb_collocation_test' + name_usv +'AMSR2MMDB_combined.nc'

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
        #caluclate filelist ofr orbits on specific day
        filelist = get_filelist_amsr2_l1r(adir_l1r, usv_day)

        #the more recent data is in daily directories, so easy to search
        #the older data, pre 2018 is in monthly directories so only search files for day
        print(usv_day.data,'numfiles:',len(filelist))
        x,y,z = [],[],[]
        for file in filelist:
            xlat,xlon,sat_time,var_data=get_orbital_data_amsr2(iusv,file)
            x = xlon.data
            y = xlat.data
            z = var_data.data
            lons,lats,data = x,y,z
            swath_def = SwathDefinition(lons, lats)
            result1 = resample_nearest(swath_def, data, area_def, radius_of_influence=20000, fill_value=None)
            da = xr.DataArray(result1,name='sat_data',coords={'lat':rlat,'lon':rlon},dims=('lat','lon'))
            subset = da.sel(lat = slice(maxlat,minlat),lon=slice(minlon,maxlon))
            num_obs = np.isfinite(subset).sum()

            if num_obs<1:  #no collocations so go to next orbit
                continue

            # drop points outside of box
            usv_min_lon, usv_max_lon = ds_usv.lon.min().data - .5, ds_usv.lon.max().data + .5
            usv_min_lat, usv_max_lat = ds_usv.lat.min().data - .5, ds_usv.lat.max().data + .5
            cond = (xlon >= usv_min_lon) & (xlon <= usv_max_lon)
            sub_lon = xlon.where(cond)
            cond = (xlat >= usv_min_lat) & (xlat <= usv_max_lat)
            sub_lat = xlat.where(cond)

            ph0 = var_data.phony_dim_0
            ph1 = var_data.phony_dim_1
            tem_time = sat_time
            ds = xr.Dataset({'time': (['phony_dim_0'], tem_time),
                             'tb': (['phony_dim_0', 'phony_dim_1'], var_data.data),
                             'lat': (['phony_dim_0', 'phony_dim_1'], sub_lat.data),
                             'lon': (['phony_dim_0', 'phony_dim_1'], sub_lon.data)},
                            coords={'phony_dim_0': (['phony_dim_0'], ph0),
                                    'phony_dim_1': (['phony_dim_1'], ph1)})

            #stack xarray dataset then drop lon == nan
            ds2 = ds.stack(z=('phony_dim_0', 'phony_dim_1')).reset_index('z')
            #drop nan
            ds_dropa = ds2.where(np.isfinite(ds2.lon), drop=True)
            ds_drop = ds_dropa.where(np.isfinite(ds_dropa.lat), drop=True)

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
                    ds_usv.amsr2_name[iusv] = file
                    ds_usv.amsr2_dist[iusv] = rdist
                    ds_usv.amsr2_scan[iusv] = ii
                    ds_usv.amsr2_cell[iusv] = jj

        usv_day += np.timedelta64(1,'D')
#    df = xr.DataArray(file_save,name='filenames')
    ds_usv = ds_usv.rename({'TEMP_CTD_MEAN':'insitu.sea_surface_temperature','TEMP_CTD_STDDEV':'insitu.sst_uncertainty',
                             'TEMP_AIR_MEAN':'insitu.air_temperature','VWND_MEAN':'insitu.vwnd','UWND_MEAN':'insitu.uwnd',
                             'WAVE_SIGNIFICANT_HEIGHT':'insitu.sig_wave_height','SAL_MEAN':'insitu.salinity','CHLOR_MEAN':'insitu.chlor',
                             'BARO_PRES_MEAN':'insitu.baro_pres','RH_MEAN':'insitu.rel_humidity','GUST_WND_MEAN':'insitu.gust_wind',
                             'lat':'insitu.lat','lon':'insitu.lon','time':'insitu.time'})

    ds_usv.to_netcdf(fileout)



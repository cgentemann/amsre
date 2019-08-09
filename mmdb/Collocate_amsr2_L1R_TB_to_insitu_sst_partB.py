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
from read_routines import read_usv,get_orbital_data_amsr2

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
for num_usv in range(input_iusv_start,input_iusv_end):
    ds_usv, usv_name = read_usv(adir_usv,num_usv)
    filelist = adir_usv + usv_name + 'AMSR2MMDB_filesave2_testing.nc'
    fileout = adir_usv + usv_name + 'AMSR2MMDB_usv2_testing.nc'
    df = xr.open_dataset(filelist)
    for file2 in df.filenames.data:
        print(file2)
        xlat, xlon, sat_time, var_data = get_orbital_data_amsr2(num_usv, file2)
        ph0 = var_data.phony_dim_0
        ph1 = var_data.phony_dim_1
        tem_time = sat_time
        ds = xr.Dataset({'time': (['phony_dim_0'], tem_time),
                         'tb': (['phony_dim_0', 'phony_dim_1'], var_data.data),
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

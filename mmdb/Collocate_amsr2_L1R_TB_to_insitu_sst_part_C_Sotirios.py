#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 09:21:51 2020

@author: ssk
"""

import sys
import numpy as np
#import xarray as xr
#from glob import glob
#from pyresample import image, geometry, load_area, save_quicklook, SwathDefinition, area_def2basemap
#from pyresample.kd_tree import resample_nearest
import gzip
import shutil
#from scipy import spatial
#sys.path.append('./subroutines/')
from read_routines import read_usv
from netCDF4 import Dataset
import time

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
#input_iusv_start = int(input("Enter start cruise processing number 0-10: "))
#input_iusv_end = int(input("Enter stop cruise processing number 0-10: "))
#adir_usv = str(input("Enter directory for USV data: "))

input_iusv_start = 2
input_iusv_end = 4
adir_usv = '/data/nis/ocean/home/ssk/Projects/SAILDRONE_MATCHUP/saildrone_data/'

# Loop through Saildrone campaigns
# --------------------------------
for iusv in range(input_iusv_start,input_iusv_end):
    ds_usv, usv_name = read_usv(adir_usv,iusv)
    print('usv', usv_name)
    fileout = adir_usv + 'mmdb_collocation_test_L2P/' + usv_name + '_L2P_MMDB_combined.nc'
    ds_usv = Dataset(fileout,'r')
# now collocation with orbital data is finished.  re-open file and create the mean values for each matchup so there aren't repeates

#    fileout_norepeat = fileout[:-3]+'_norepeats_testing.nc'
    fileout_norepeat = '/media/ssk/Seagate Expansion Drive/saildrone_amsr/L2P_matchups/saildrone_' + usv_name + '_L2P_MMDB.nc'
    dsst,dair,dwav,dpres,dsstu,dvwnd,duwnd,dsal,dchl,dgwnd,dlat,dlon,dut,dfilename,dscan,dcell = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    
    # Get the filenames from the netCDF file
    files = ds_usv.variables['amsr2_name']
    file_list = []
    for i in range(0,len(files)):
        file_list.append(str(files[i]))
        
    file_list_all = np.array(file_list)
    # Get unique entries from filename list
    file_list_unique = np.unique(np.array(file_list))[1:-1]
    
    # Loop through filenames and get the corresponding variables
    # ----------------------------------------------------------
    for file in file_list_unique:
        idx_file = np.where(file == file_list_all)
        
        temp_cell = ds_usv.variables['amsr2_cell'][idx_file].data
        temp_scan = ds_usv.variables['amsr2_scan'][idx_file].data
        temp_dist = ds_usv.variables['amsr2_dist'][idx_file].data
        
        temp_sst = ds_usv.variables['insitu.sea_surface_temperature'][idx_file].data
        temp_sst[temp_sst < -1000] = np.nan
        temp_air = ds_usv.variables['insitu.air_temperature'][idx_file].data
        temp_air[temp_air < -1000] = np.nan
        temp_wav = ds_usv.variables['insitu.sig_wave_height'][idx_file].data
        temp_wav[temp_wav < -1000] = np.nan
        temp_pres = ds_usv.variables['insitu.baro_pres'][idx_file].data
        temp_pres[temp_pres < -1000] = np.nan
        temp_sstu = ds_usv.variables['insitu.sst_uncertainty'][idx_file].data
        temp_sstu[temp_sstu < -1000] = np.nan
        temp_vwnd = ds_usv.variables['insitu.vwnd'][idx_file].data
        temp_vwnd[temp_vwnd < -1000] = np.nan
        temp_uwnd = ds_usv.variables['insitu.uwnd'][idx_file].data
        temp_uwnd[temp_uwnd < -1000] = np.nan
        temp_sal = ds_usv.variables['insitu.salinity'][idx_file].data
        temp_sal[temp_sal < -1000] = np.nan
        temp_chl = ds_usv.variables['insitu.chlor'][idx_file].data
        temp_chl[temp_chl < -1000] = np.nan
        temp_gwnd = ds_usv.variables['insitu.gust_wind'][idx_file].data
        temp_gwnd[temp_gwnd < -1000] = np.nan
        temp_lat = ds_usv.variables['insitu.lat'][idx_file].data
        temp_lon = ds_usv.variables['insitu.lon'][idx_file].data
        temp_time = ds_usv.variables['insitu.time'][idx_file].data
        
        cell_start = int(min(temp_cell))
        cell_end = int(max(temp_cell))
        scan_start = int(min(temp_scan))
        scan_end = int(max(temp_scan))
        
#        x=1
        
        for i in range(cell_start,cell_end+1):
            for j in range(scan_start,scan_end+1):
                idx_obs = np.intersect1d(np.where(temp_cell == i),np.where(temp_scan == j))
                dsst.append(np.nanmean(temp_sst[idx_obs]))
                dair.append(np.nanmean(temp_air[idx_obs]))
                dwav.append(np.nanmean(temp_wav[idx_obs]))
                dpres.append(np.nanmean(temp_pres[idx_obs]))
                dsstu.append(np.nanmean(temp_sstu[idx_obs]))
                dvwnd.append(np.nanmean(temp_vwnd[idx_obs]))
                duwnd.append(np.nanmean(temp_uwnd[idx_obs]))
                dsal.append(np.nanmean(temp_sal[idx_obs]))
                dchl.append(np.nanmean(temp_chl[idx_obs]))
                dgwnd.append(np.nanmean(temp_gwnd[idx_obs]))
                dlat.append(np.nanmean(temp_lat[idx_obs]))
                dlon.append(np.nanmean(temp_lon[idx_obs]))
                dut.append(np.nanmean(temp_time[idx_obs]))
                dfilename.append(str(file))
                dcell.append(i)
                dscan.append(j)
    
    ds_usv.close()
    
    # With the method above, we will end up with NaNs in case there is not any data with the i,j while going throught the loops
    # So we need to clear our lists from NaNs
    # ---------------------------------------
    # First convert lists to arrays
    dsst = np.array(dsst)
    dair = np.array(dair)
    dwav = np.array(dwav)
    dpres = np.array(dpres)
    dsstu = np.array(dsstu)
    dvwnd = np.array(dvwnd)
    duwnd = np.array(duwnd)
    dsal = np.array(dsal)
    dchl = np.array(dchl)
    dgwnd = np.array(dgwnd) 
    dlat = np.array(dlat)
    dlon = np.array(dlon)
    dut = np.array(dut)
    dfilename = np.array(dfilename)
    dscan = np.array(dscan)
    dcell = np.array(dcell)
    
    # Then find the NaN indices
    idx_non_nan = np.squeeze(np.argwhere(~np.isnan(dut)))
    
    dsst = dsst[idx_non_nan]
    dair = dair[idx_non_nan]
    dwav = dwav[idx_non_nan]
    dpres = dpres[idx_non_nan]
    dsstu = dsstu[idx_non_nan]
    dvwnd = dvwnd[idx_non_nan]
    duwnd = duwnd[idx_non_nan]
    dsal = dsal[idx_non_nan]
    dchl = dchl[idx_non_nan]
    dgwnd = dgwnd[idx_non_nan] 
    dlat = dlat[idx_non_nan]
    dlon = dlon[idx_non_nan]
    dut = dut[idx_non_nan]
    dfilename = dfilename[idx_non_nan]
    dscan = dscan[idx_non_nan]
    dcell = dcell[idx_non_nan]
    
    # Then sort remaining observation by time
    idx_sort = np.argsort(dut)
    
    dsst = dsst[idx_sort]
    dair = dair[idx_sort]
    dwav = dwav[idx_sort]
    dpres = dpres[idx_sort]
    dsstu = dsstu[idx_sort]
    dvwnd = dvwnd[idx_sort]
    duwnd = duwnd[idx_sort]
    dsal = dsal[idx_sort]
    dchl = dchl[idx_sort]
    dgwnd = dgwnd[idx_sort] 
    dlat = dlat[idx_sort]
    dlon = dlon[idx_sort]
    dut = dut[idx_sort]
    dfilename = dfilename[idx_sort]
    dscan = dscan[idx_sort]
    dcell = dcell[idx_sort]    
    
    
    # Create AMSR2 matchups
    unc_adj = np.zeros_like(dut)
    l2p_flags = np.zeros_like(dut)
    unc_lsc = np.zeros_like(dut)
    l2p_lat = np.zeros_like(dut)
    l2p_lon = np.zeros_like(dut)
    ql = np.zeros_like(dut)
    sst = np.zeros_like(dut)
    sst_d = np.zeros_like(dut)
    sses_bias = np.zeros_like(dut)
    sses_std = np.zeros_like(dut)
#    sst_ddt = np.zeros_like(dut)
    sst_dtu = np.zeros_like(dut)
#    sst_dt = np.zeros_like(dut)
    unc_sc = np.zeros_like(dut)
    l2p_time = np.zeros_like(dut)
    unc_u = np.zeros_like(dut)
    ws = np.zeros_like(dut)
    ws_nwp = np.zeros_like(dut)
    
    
    fname = dfilename[0]
    start_flag = 1
    for i in range(0,len(dfilename)):
        
        if start_flag == 1:
            year = dfilename[i][0:4]
            month = dfilename[i][4:6]
            day = dfilename[i][6:8]
            file_tb = '/net/isilon/ifs/arch/home/sstdev/Projects/ESA_CCI2/AMSR/L2P/v2.0/AMSR2/'+year+'/'+month+'/'+day+'/'+dfilename[i]+'.0-fv01.0.nc'
            
#            if file_tb[-3:] == '.gz':  # early files gzipped
#                with gzip.open(file_tb, 'rb') as f_in:
#                    temp_file = '/data/ssk/temp/tem_' + str(iusv) + file[-6:-3]
#                    with open(temp_file, 'wb') as f_out:
#                        shutil.copyfileobj(f_in, f_out)
#                nc_amsr = Dataset(temp_file,'r')
#            else:
#                file_tb = '/net/isilon/ifs/arch/home/sstdev/Projects/ESA_CCI2/AMSR/L1R/AMSR2/2/'+year+'/'+month+'/'+dfilename[i]
            nc_amsr = Dataset(file_tb,'r')
            
            
            
            unc_adj[i] = nc_amsr.variables['adjustment_uncertainty'][0,dscan[i],dcell[i]]
            l2p_flags[i] = nc_amsr.variables['l2p_flags'][0,dscan[i],dcell[i]].data
            unc_lsc[i] = nc_amsr.variables['large_scale_correlated_uncertainty'][0,dscan[i],dcell[i]]
            l2p_lat[i] = nc_amsr.variables['lat'][dscan[i],dcell[i]].data
            l2p_lon[i] = nc_amsr.variables['lon'][dscan[i],dcell[i]].data
            ql[i] = nc_amsr.variables['quality_level'][0,dscan[i],dcell[i]].data
            sst[i] = nc_amsr.variables['sea_surface_temperature'][0,dscan[i],dcell[i]]
            sst_d[i] = nc_amsr.variables['sea_surface_temperature_depth'][0,dscan[i],dcell[i]]
            sses_bias[i] = nc_amsr.variables['sses_bias'][0,dscan[i],dcell[i]]
            sses_std[i] = nc_amsr.variables['sses_standard_deviation'][0,dscan[i],dcell[i]]
            sst_dtu[i] = nc_amsr.variables['sst_depth_total_uncertainty'][0,dscan[i],dcell[i]]
            unc_sc[i] = nc_amsr.variables['synoptically_correlated_uncertainty'][0,dscan[i],dcell[i]]
            l2p_time[i] = nc_amsr.variables['time'][:].data + nc_amsr.variables['sst_dtime'][0,dscan[i],dcell[i]].data
#            sst_dt[i] = nc_amsr.variables['time'][:].data + nc_amsr.variables['sst_depth_dtime'][0,dscan[i],dcell[i]].data
            unc_u[i] = nc_amsr.variables['uncorrelated_uncertainty'][0,dscan[i],dcell[i]]
            ws[i] = nc_amsr.variables['wind_speed'][0,dscan[i],dcell[i]]
            ws_nwp[i] = nc_amsr.variables['wind_speed_nwp'][0,dscan[i],dcell[i]]
            
            start_flag = 0
            fname = dfilename[i]
            
            x=1
            
        else:
            if dfilename[i] == fname:
                
                unc_adj[i] = nc_amsr.variables['adjustment_uncertainty'][0,dscan[i],dcell[i]]
                l2p_flags[i] = nc_amsr.variables['l2p_flags'][0,dscan[i],dcell[i]].data
                unc_lsc[i] = nc_amsr.variables['large_scale_correlated_uncertainty'][0,dscan[i],dcell[i]]
                l2p_lat[i] = nc_amsr.variables['lat'][dscan[i],dcell[i]].data
                l2p_lon[i] = nc_amsr.variables['lon'][dscan[i],dcell[i]].data
                ql[i] = nc_amsr.variables['quality_level'][0,dscan[i],dcell[i]].data
                sst[i] = nc_amsr.variables['sea_surface_temperature'][0,dscan[i],dcell[i]]
                sst_d[i] = nc_amsr.variables['sea_surface_temperature_depth'][0,dscan[i],dcell[i]]
                sses_bias[i] = nc_amsr.variables['sses_bias'][0,dscan[i],dcell[i]]
                sses_std[i] = nc_amsr.variables['sses_standard_deviation'][0,dscan[i],dcell[i]]
                sst_dtu[i] = nc_amsr.variables['sst_depth_total_uncertainty'][0,dscan[i],dcell[i]]
                unc_sc[i] = nc_amsr.variables['synoptically_correlated_uncertainty'][0,dscan[i],dcell[i]]
                l2p_time[i] = nc_amsr.variables['time'][:].data + nc_amsr.variables['sst_dtime'][0,dscan[i],dcell[i]].data
#                sst_dt[i] = nc_amsr.variables['time'][:].data + nc_amsr.variables['sst_depth_dtime'][0,dscan[i],dcell[i]].data
                unc_u[i] = nc_amsr.variables['uncorrelated_uncertainty'][0,dscan[i],dcell[i]]
                ws[i] = nc_amsr.variables['wind_speed'][0,dscan[i],dcell[i]]
                ws_nwp[i] = nc_amsr.variables['wind_speed_nwp'][0,dscan[i],dcell[i]]
                
                fname = dfilename[i]
            else:
                nc_amsr.close()
                year = dfilename[i][0:4]
                month = dfilename[i][4:6]
                day = dfilename[i][6:8]
                file_tb = '/net/isilon/ifs/arch/home/sstdev/Projects/ESA_CCI2/AMSR/L2P/v2.0/AMSR2/'+year+'/'+month+'/'+day+'/'+dfilename[i]+'.0-fv01.0.nc'
                nc_amsr = Dataset(file_tb,'r')
            
                unc_adj[i] = nc_amsr.variables['adjustment_uncertainty'][0,dscan[i],dcell[i]]
                l2p_flags[i] = nc_amsr.variables['l2p_flags'][0,dscan[i],dcell[i]].data
                unc_lsc[i] = nc_amsr.variables['large_scale_correlated_uncertainty'][0,dscan[i],dcell[i]]
                l2p_lat[i] = nc_amsr.variables['lat'][dscan[i],dcell[i]].data
                l2p_lon[i] = nc_amsr.variables['lon'][dscan[i],dcell[i]].data
                ql[i] = nc_amsr.variables['quality_level'][0,dscan[i],dcell[i]].data
                sst[i] = nc_amsr.variables['sea_surface_temperature'][0,dscan[i],dcell[i]]
                sst_d[i] = nc_amsr.variables['sea_surface_temperature_depth'][0,dscan[i],dcell[i]]
                sses_bias[i] = nc_amsr.variables['sses_bias'][0,dscan[i],dcell[i]]
                sses_std[i] = nc_amsr.variables['sses_standard_deviation'][0,dscan[i],dcell[i]]
                sst_dtu[i] = nc_amsr.variables['sst_depth_total_uncertainty'][0,dscan[i],dcell[i]]
                unc_sc[i] = nc_amsr.variables['synoptically_correlated_uncertainty'][0,dscan[i],dcell[i]]
                l2p_time[i] = nc_amsr.variables['time'][:].data + nc_amsr.variables['sst_dtime'][0,dscan[i],dcell[i]].data
#                sst_dt[i] = nc_amsr.variables['time'][:].data + nc_amsr.variables['sst_depth_dtime'][0,dscan[i],dcell[i]].data
                unc_u[i] = nc_amsr.variables['uncorrelated_uncertainty'][0,dscan[i],dcell[i]]
                ws[i] = nc_amsr.variables['wind_speed'][0,dscan[i],dcell[i]]
                ws_nwp[i] = nc_amsr.variables['wind_speed_nwp'][0,dscan[i],dcell[i]]            
                
                fname = dfilename[i]
    
    
    
                
    # Write to netCDF file
    # --------------------
    nc_out = Dataset(fileout_norepeat,'w',file_format='NETCDF4')
    nc_out.createDimension('matchup_count',len(dfilename))
    
    # Create global variables
    nc_out.description = "SAILDRONE USV in situ data with collocated L2P SST data"
    nc_out.project_leader = "Chelle Gentemann"
    nc_out.project_leader_email = "cgentemann@esr.org"
    nc_out.creator = "Sotirios Skarpalezos"
    nc_out.creator_email = "ssk@dmi.dk"
    nc_out.creator_institution = "Danish Meteorological Institute"
    nc_out.history = "Created " + time.ctime(time.time())
    
    # Create variables
    nc_sst = nc_out.createVariable('saildrone-sst.insitu.sea_surface_temperature', 'f8',('matchup_count'), fill_value = -1e34)
    nc_air = nc_out.createVariable('saildrone-sst.insitu.air_temperature', 'f8',('matchup_count'), fill_value = -1e34)
    nc_wav = nc_out.createVariable('saildrone-sst.insitu.sig_wave_height', 'f8',('matchup_count'), fill_value = -1e34)
    nc_pres = nc_out.createVariable('saildrone-sst.insitu.baro_press', 'f8',('matchup_count'), fill_value = -1e34)
    nc_sstu = nc_out.createVariable('saildrone-sst.insitu.sst_uncertainty', 'f8',('matchup_count'), fill_value = -1e34)
    nc_uwnd = nc_out.createVariable('saildrone-sst.insitu.5m_east_wind_component', 'f8',('matchup_count'), fill_value = -1e34)
    nc_vwnd = nc_out.createVariable('saildrone-sst.insitu.5m_north_wind_component', 'f8',('matchup_count'), fill_value = -1e34)
    nc_sal = nc_out.createVariable('saildrone-sst.insitu.sea_surface_salinity', 'f8',('matchup_count'), fill_value = -1e34)
    nc_chl = nc_out.createVariable('saildrone-sst.insitu.chlor', 'f8',('matchup_count'), fill_value = -1e34)
    nc_gwnd = nc_out.createVariable('saildrone-sst.insitu.gust_wind', 'f8',('matchup_count'), fill_value = -1e34)
    nc_lat = nc_out.createVariable('saildrone-sst.insitu.lat', 'f8',('matchup_count'), fill_value = -1e34)
    nc_lon = nc_out.createVariable('saildrone-sst.insitu.lon', 'f8',('matchup_count'), fill_value = -1e34)
    nc_ut = nc_out.createVariable('saildrone-sst.insitu.time', 'f8',('matchup_count'), fill_value = -1e34)
    nc_filename = nc_out.createVariable('saildrone-sst.insitu.file_name', 'S1',('matchup_count'), fill_value = -1e34)
    nc_scan = nc_out.createVariable('saildrone-sst.insitu.amsr2.matchup.line', 'f8',('matchup_count'), fill_value = -1e34)
    nc_cell = nc_out.createVariable('saildrone-sst.insitu.amsr2.matchup.elem', 'f8',('matchup_count'), fill_value = -1e34)
    
    nc_unc_adj = nc_out.createVariable('dmi.amsr2.adjustment_uncertainty', 'i2',('matchup_count'))
    nc_l2p_flags = nc_out.createVariable('dmi.amsr2.l2p_flags', 'i2',('matchup_count'))
    nc_unc_lsc = nc_out.createVariable('dmi.amsr2.large_scale_correlated_uncertainty', 'i2',('matchup_count'))
    nc_l2p_lat = nc_out.createVariable('dmi.amsr2.lat', 'f4',('matchup_count'))
    nc_l2p_lon = nc_out.createVariable('dmi.amsr2.lon', 'f4',('matchup_count'))
    nc_ql = nc_out.createVariable('dmi.amsr2.quality_level', 'i1',('matchup_count'))
    nc_sst_l2p = nc_out.createVariable('dmi.amsr2.sea_surface_temperature', 'i2',('matchup_count'))#, fill_value = nc_amsr.variables['sea_surface_temperature']._FillValue)
    nc_sst_d = nc_out.createVariable('dmi.amsr2.sea_surface_temperature_depth', 'i2',('matchup_count'))
    nc_sses_bias = nc_out.createVariable('dmi.amsr2.sses_bias', 'i1',('matchup_count'))
    nc_sses_std = nc_out.createVariable('dmi.amsr2.sses_standard_deviation', 'i2',('matchup_count'))
    nc_sst_dtu = nc_out.createVariable('dmi.amsr2.sst_depth_total_uncertainty', 'i2',('matchup_count'))
    nc_unc_sc = nc_out.createVariable('dmi.amsr2.synoptically_correlated_uncertainty', 'i2',('matchup_count'))
    nc_l2p_time = nc_out.createVariable('dmi.amsr2.time', 'i4',('matchup_count'))
    nc_unc_u = nc_out.createVariable('dmi.amsr2.uncorrelated_uncertainty', 'i2',('matchup_count'))
    nc_ws = nc_out.createVariable('dmi.amsr2.wind_speed', 'i4',('matchup_count'))
    nc_ws_nwp = nc_out.createVariable('dmi.amsr2.wind_speed_nwp', 'i4',('matchup_count'))
    
#    # Calculate offsets and scale factors
#    n = 16
#    n8 = 8
#    n32 = 32
#    
#    unc_adj_scale_factor = (max(unc_adj) - min(unc_adj)) / (2 ** n - 1)
#    unc_adj_add_offset = min(unc_adj) + 2 ** (n - 1) * unc_adj_scale_factor
#    
#    unc_lsc_scale_factor = (max(unc_lsc) - min(unc_lsc)) / (2 ** n - 1)
#    unc_lsc_add_offset = min(unc_lsc) + 2 ** (n - 1) * unc_lsc_scale_factor
#    
#    sst_scale_factor = (max(sst) - min(sst)) / (2 ** n - 1)
#    sst_add_offset = min(sst) + 2 ** (n - 1) * sst_scale_factor
#    
#    sst_d_scale_factor = (max(sst_d) - min(sst_d)) / (2 ** n - 1)
#    sst_d_add_offset = min(sst_d) + 2 ** (n - 1) * sst_d_scale_factor
#    
#    sses_bias_scale_factor = (max(sses_bias) - min(sses_bias)) / (2 ** n8 - 1)
#    sses_bias_add_offset = min(sses_bias) + 2 ** (n8 - 1) * sses_bias_scale_factor
#    
#    sses_std_scale_factor = (max(sses_std) - min(sses_std)) / (2 ** n - 1)
#    sses_std_add_offset = min(sses_std) + 2 ** (n - 1) * sses_std_scale_factor
#    
#    sst_dtu_scale_factor = (max(sst_dtu) - min(sst_dtu)) / (2 ** n - 1)
#    sst_dtu_add_offset = min(sst_dtu) + 2 ** (n - 1) * sst_dtu_scale_factor
#    
#    unc_sc_scale_factor = (max(unc_sc) - min(unc_sc)) / (2 ** n - 1)
#    unc_sc_add_offset = min(unc_sc) + 2 ** (n - 1) * unc_sc_scale_factor
#    
#    unc_u_scale_factor = (max(unc_u) - min(unc_u)) / (2 ** n - 1)
#    unc_u_add_offset = min(unc_u) + 2 ** (n - 1) * unc_u_scale_factor
#    
#    ws_scale_factor = (max(ws) - min(ws)) / (2 ** n32 - 1)
#    ws_add_offset = min(ws) + 2 ** (n32 - 1) * ws_scale_factor
#    
#    ws_nwp_scale_factor = (max(ws_nwp) - min(ws_nwp)) / (2 ** n32 - 1)
#    ws_nwp_add_offset = min(ws_nwp) + 2 ** (n32 - 1) * ws_nwp_scale_factor
    
    
    # Fill variables
    nc_sst.long_name = 'Seawater temperature'
    nc_sst.units = 'degrees_c'
    nc_sst[:] = dsst

    nc_air.long_name = 'Air temperature'
    nc_air.units = 'degrees_c'
    nc_air[:] = dair
    
    nc_wav.long_name = 'Significant wave height'
    nc_wav.units = 'm'
    nc_wav[:] = dwav
    
    nc_pres.long_name = 'Air pressure'
    nc_pres.units = 'hPa'
    nc_pres[:] = dpres
    
    nc_sstu.long_name = 'Seawater temperature SD'
    nc_sstu.units = 'degrees_c'
    nc_sstu[:] = dsstu
    
    nc_uwnd.long_name = 'Eastward wind speed'
    nc_uwnd.units = 'm s-1'
    nc_uwnd[:] = duwnd
    
    nc_vwnd.long_name = 'Northward wind speed'
    nc_vwnd.units = 'm s-1'
    nc_vwnd[:] = dvwnd
    
    nc_sal.long_name = 'Seawater salinity'
    nc_sal.units = 'PSS-78'
    nc_sal[:] = dsal
    
    nc_chl.long_name = 'Chlorophyll concentration'
    nc_chl.units = 'microgram L-1'
    nc_chl[:] = dchl
    
    nc_gwnd.long_name = 'Wind gust speed'
    nc_gwnd.units = 'knots'
    nc_gwnd[:] = dgwnd
    
    nc_lat.long_name = 'Latitude'
    nc_lat.units = 'degrees_north'
    nc_lat[:] = dlat
    
    nc_lon.long_name = 'Longitude'
    nc_lon.units = 'degrees_east'
    nc_lon[:] = dlon
    
    if (usv_name[11::] == '1033' or usv_name[11::] == '1035'):
        nc_ut.long_name = 'Time'
        nc_ut.units = 'minutes since 2019-05-14 00:00:00'
        nc_ut[:] = dut
    else:
        nc_ut.long_name = 'Time'
        nc_ut.units = 'seconds since 1970-01-01 00:00:00'
        nc_ut[:] = dut        
    
    nc_filename.long_name = 'AMSR2 L2 source file name'
    nc_filename[:] = dfilename
    
    nc_scan.long_name = 'AMSR2 L2 source file scan line'
    nc_scan[:] = dscan
    
    nc_cell.long_name = 'AMSR2 L2 source file cell'
    nc_cell[:] = dcell
    
    
    
    nc_unc_adj.long_name = nc_amsr.variables['adjustment_uncertainty'].long_name
    nc_unc_adj.units = nc_amsr.variables['adjustment_uncertainty'].units
    nc_unc_adj.add_offset = nc_amsr.variables['adjustment_uncertainty'].add_offset#unc_adj_add_offset
    nc_unc_adj.scale_factor = nc_amsr.variables['adjustment_uncertainty'].scale_factor#unc_adj_scale_factor
    nc_unc_adj.valid_min = nc_amsr.variables['adjustment_uncertainty'].valid_min
    nc_unc_adj.valid_max = nc_amsr.variables['adjustment_uncertainty'].valid_max
    nc_unc_adj.comment = nc_amsr.variables['adjustment_uncertainty'].comment
    nc_unc_adj.references = nc_amsr.variables['adjustment_uncertainty'].references
    nc_unc_adj.coordinates = nc_amsr.variables['adjustment_uncertainty'].coordinates
    nc_unc_adj.correlation_length_scale = nc_amsr.variables['adjustment_uncertainty'].correlation_length_scale
    nc_unc_adj.correlation_time_scale = nc_amsr.variables['adjustment_uncertainty'].correlation_time_scale
    nc_unc_adj[:] = unc_adj
    
    nc_l2p_flags.long_name = nc_amsr.variables['l2p_flags'].long_name
    nc_l2p_flags.valid_min = nc_amsr.variables['l2p_flags'].valid_min
    nc_l2p_flags.valid_max = nc_amsr.variables['l2p_flags'].valid_max
    nc_l2p_flags.flag_meanings = nc_amsr.variables['l2p_flags'].flag_meanings
    nc_l2p_flags.flag_masks = nc_amsr.variables['l2p_flags'].flag_masks
    nc_l2p_flags.comment = nc_amsr.variables['l2p_flags'].comment
    nc_l2p_flags.coordinates = nc_amsr.variables['l2p_flags'].coordinates
    nc_l2p_flags[:] = l2p_flags
    
    nc_unc_lsc.long_name = nc_amsr.variables['large_scale_correlated_uncertainty'].long_name
    nc_unc_lsc.units = nc_amsr.variables['large_scale_correlated_uncertainty'].units
    nc_unc_lsc.add_offset = nc_amsr.variables['large_scale_correlated_uncertainty'].add_offset#unc_lsc_add_offset
    nc_unc_lsc.scale_factor = nc_amsr.variables['large_scale_correlated_uncertainty'].scale_factor#unc_lsc_scale_factor
    nc_unc_lsc.valid_min = nc_amsr.variables['large_scale_correlated_uncertainty'].valid_min
    nc_unc_lsc.valid_max = nc_amsr.variables['large_scale_correlated_uncertainty'].valid_max
    nc_unc_lsc.comment = nc_amsr.variables['large_scale_correlated_uncertainty'].comment
    nc_unc_lsc.references = nc_amsr.variables['large_scale_correlated_uncertainty'].references
    nc_unc_lsc.coordinates = nc_amsr.variables['large_scale_correlated_uncertainty'].coordinates
    nc_unc_lsc[:] = unc_lsc
    
    nc_l2p_lat.least_significant_digit = nc_amsr.variables['lat'].least_significant_digit
    nc_l2p_lat.reference_datum = nc_amsr.variables['lat'].reference_datum
    nc_l2p_lat.long_name = nc_amsr.variables['lat'].long_name
    nc_l2p_lat.standard_name = nc_amsr.variables['lat'].standard_name
    nc_l2p_lat.units = nc_amsr.variables['lat'].units
    nc_l2p_lat.valid_min = nc_amsr.variables['lat'].valid_min
    nc_l2p_lat.valid_max = nc_amsr.variables['lat'].valid_max
    nc_l2p_lat[:] = l2p_lat
    
    nc_l2p_lon.least_significant_digit = nc_amsr.variables['lon'].least_significant_digit
    nc_l2p_lon.reference_datum = nc_amsr.variables['lon'].reference_datum
    nc_l2p_lon.long_name = nc_amsr.variables['lon'].long_name
    nc_l2p_lon.standard_name = nc_amsr.variables['lon'].standard_name
    nc_l2p_lon.units = nc_amsr.variables['lon'].units
    nc_l2p_lon.valid_min = nc_amsr.variables['lon'].valid_min
    nc_l2p_lon.valid_max = nc_amsr.variables['lon'].valid_max
    nc_l2p_lon[:] = l2p_lon
    
    nc_ql.long_name = nc_amsr.variables['quality_level'].long_name
    nc_ql.valid_min = nc_amsr.variables['quality_level'].valid_min
    nc_ql.valid_max = nc_amsr.variables['quality_level'].valid_max
    nc_ql.flag_meanings = nc_amsr.variables['quality_level'].flag_meanings
    nc_ql.flag_values = nc_amsr.variables['quality_level'].flag_values
    nc_ql.comment = nc_amsr.variables['quality_level'].comment
    nc_ql.coordinates = nc_amsr.variables['quality_level'].coordinates
    nc_ql[:] = ql
    
    nc_sst_l2p.long_name = nc_amsr.variables['sea_surface_temperature'].long_name
    nc_sst_l2p.standard_name = nc_amsr.variables['sea_surface_temperature'].standard_name
    nc_sst_l2p.units = nc_amsr.variables['sea_surface_temperature'].units
    nc_sst_l2p.add_offset = nc_amsr.variables['sea_surface_temperature'].add_offset#sst_add_offset
    nc_sst_l2p.scale_factor = nc_amsr.variables['sea_surface_temperature'].scale_factor#sst_scale_factor
    nc_sst_l2p.valid_min = nc_amsr.variables['sea_surface_temperature'].valid_min
    nc_sst_l2p.valid_max = nc_amsr.variables['sea_surface_temperature'].valid_max
    nc_sst_l2p.comment = nc_amsr.variables['sea_surface_temperature'].comment
    nc_sst_l2p.references = nc_amsr.variables['sea_surface_temperature'].references
    nc_sst_l2p.source = nc_amsr.variables['sea_surface_temperature'].source
    nc_sst_l2p.depth = nc_amsr.variables['sea_surface_temperature'].depth
    nc_sst_l2p.coordinates = nc_amsr.variables['sea_surface_temperature'].coordinates
    nc_sst_l2p[:] = sst


    nc_sst_d.long_name = nc_amsr.variables['sea_surface_temperature_depth'].long_name
    nc_sst_d.standard_name = nc_amsr.variables['sea_surface_temperature_depth'].standard_name
    nc_sst_d.units = nc_amsr.variables['sea_surface_temperature_depth'].units
    nc_sst_d.add_offset = nc_amsr.variables['sea_surface_temperature_depth'].add_offset#sst_d_add_offset
    nc_sst_d.scale_factor = nc_amsr.variables['sea_surface_temperature_depth'].scale_factor#sst_d_scale_factor
    nc_sst_d.valid_min = nc_amsr.variables['sea_surface_temperature_depth'].valid_min
    nc_sst_d.valid_max = nc_amsr.variables['sea_surface_temperature_depth'].valid_max
    nc_sst_d.comment = nc_amsr.variables['sea_surface_temperature_depth'].comment
    nc_sst_d.references = nc_amsr.variables['sea_surface_temperature_depth'].references
    nc_sst_d.source = nc_amsr.variables['sea_surface_temperature_depth'].source
    nc_sst_d.depth = nc_amsr.variables['sea_surface_temperature_depth'].depth
    nc_sst_d.coordinates = nc_amsr.variables['sea_surface_temperature_depth'].coordinates
    nc_sst_d[:] = sst_d
    
    nc_sses_bias.long_name = nc_amsr.variables['sses_bias'].long_name
    nc_sses_bias.units = nc_amsr.variables['sses_bias'].units
    nc_sses_bias.add_offset = nc_amsr.variables['sses_bias'].add_offset#sses_bias_add_offset
    nc_sses_bias.scale_factor = nc_amsr.variables['sses_bias'].scale_factor#sses_bias_scale_factor
    nc_sses_bias.valid_min = nc_amsr.variables['sses_bias'].valid_min
    nc_sses_bias.valid_max = nc_amsr.variables['sses_bias'].valid_max
    nc_sses_bias.comment = nc_amsr.variables['sses_bias'].comment
    nc_sses_bias.coordinates = nc_amsr.variables['sses_bias'].coordinates
    nc_sses_bias[:] = sses_bias
    
    nc_sses_std.long_name = nc_amsr.variables['sses_standard_deviation'].long_name
    nc_sses_std.units = nc_amsr.variables['sses_standard_deviation'].units
    nc_sses_std.add_offset = nc_amsr.variables['sses_standard_deviation'].add_offset#sses_std_add_offset
    nc_sses_std.scale_factor = nc_amsr.variables['sses_standard_deviation'].scale_factor#sses_std_scale_factor
    nc_sses_std.valid_min = nc_amsr.variables['sses_standard_deviation'].valid_min
    nc_sses_std.valid_max = nc_amsr.variables['sses_standard_deviation'].valid_max
    nc_sses_std.comment = nc_amsr.variables['sses_standard_deviation'].comment
    nc_sses_std.coordinates = nc_amsr.variables['sses_standard_deviation'].coordinates
    nc_sses_std[:] = sses_std
    
    nc_sst_dtu.long_name = nc_amsr.variables['sst_depth_total_uncertainty'].long_name
    nc_sst_dtu.units = nc_amsr.variables['sst_depth_total_uncertainty'].units
    nc_sst_dtu.add_offset = nc_amsr.variables['sst_depth_total_uncertainty'].add_offset#sst_dtu_add_offset
    nc_sst_dtu.scale_factor = nc_amsr.variables['sst_depth_total_uncertainty'].scale_factor#sst_dtu_scale_factor
    nc_sst_dtu.valid_min = nc_amsr.variables['sst_depth_total_uncertainty'].valid_min
    nc_sst_dtu.valid_max = nc_amsr.variables['sst_depth_total_uncertainty'].valid_max
    nc_sst_dtu.comment = nc_amsr.variables['sst_depth_total_uncertainty'].comment
    nc_sst_dtu.coordinates = nc_amsr.variables['sst_depth_total_uncertainty'].coordinates
    nc_sst_dtu[:] = sst_dtu
    
    nc_unc_sc.long_name = nc_amsr.variables['synoptically_correlated_uncertainty'].long_name
    nc_unc_sc.units = nc_amsr.variables['synoptically_correlated_uncertainty'].units
    nc_unc_sc.add_offset = nc_amsr.variables['synoptically_correlated_uncertainty'].add_offset#unc_sc_add_offset
    nc_unc_sc.scale_factor = nc_amsr.variables['synoptically_correlated_uncertainty'].scale_factor#unc_sc_scale_factor
    nc_unc_sc.valid_min = nc_amsr.variables['synoptically_correlated_uncertainty'].valid_min
    nc_unc_sc.valid_max = nc_amsr.variables['synoptically_correlated_uncertainty'].valid_max
    nc_unc_sc.comment = nc_amsr.variables['synoptically_correlated_uncertainty'].comment
    nc_unc_sc.correlation_length_scale = nc_amsr.variables['synoptically_correlated_uncertainty'].correlation_length_scale
    nc_unc_sc.correlation_time_scale = nc_amsr.variables['synoptically_correlated_uncertainty'].correlation_time_scale
    nc_unc_sc.references = nc_amsr.variables['synoptically_correlated_uncertainty'].references
    nc_unc_sc.coordinates = nc_amsr.variables['synoptically_correlated_uncertainty'].coordinates
    nc_unc_sc[:] = unc_sc
    
    nc_l2p_time.calendar = nc_amsr.variables['time'].calendar
    nc_l2p_time.standard_name = nc_amsr.variables['time'].standard_name
    nc_l2p_time.units = nc_amsr.variables['time'].units
    nc_l2p_time[:] = l2p_time
    
    nc_unc_u.long_name = nc_amsr.variables['uncorrelated_uncertainty'].long_name
    nc_unc_u.units = nc_amsr.variables['uncorrelated_uncertainty'].units
    nc_unc_u.add_offset = nc_amsr.variables['uncorrelated_uncertainty'].add_offset#unc_u_add_offset
    nc_unc_u.scale_factor = nc_amsr.variables['uncorrelated_uncertainty'].scale_factor#unc_u_scale_factor
    nc_unc_u.valid_min = nc_amsr.variables['uncorrelated_uncertainty'].valid_min
    nc_unc_u.valid_max = nc_amsr.variables['uncorrelated_uncertainty'].valid_max
    nc_unc_u.comment = nc_amsr.variables['uncorrelated_uncertainty'].comment
    nc_unc_u.references = nc_amsr.variables['uncorrelated_uncertainty'].references
    nc_unc_u.coordinates = nc_amsr.variables['uncorrelated_uncertainty'].coordinates
    nc_unc_u[:] = unc_u
    
    nc_ws.long_name = nc_amsr.variables['wind_speed'].long_name
    nc_ws.standard_name = nc_amsr.variables['wind_speed'].standard_name
    nc_ws.units = nc_amsr.variables['wind_speed'].units
    nc_ws.add_offset = nc_amsr.variables['wind_speed'].add_offset#ws_add_offset
    nc_ws.scale_factor = nc_amsr.variables['wind_speed'].scale_factor#ws_scale_factor
    nc_ws.valid_min = nc_amsr.variables['wind_speed'].valid_min
    nc_ws.valid_max = nc_amsr.variables['wind_speed'].valid_max
    nc_ws.comment = nc_amsr.variables['wind_speed'].comment
    nc_ws.references = nc_amsr.variables['wind_speed'].references
    nc_ws.time_offset = nc_amsr.variables['wind_speed'].time_offset
    nc_ws.height = nc_amsr.variables['wind_speed'].height
    nc_ws.coordinates = nc_amsr.variables['wind_speed'].coordinates
    nc_ws[:] = ws
    
    nc_ws_nwp.long_name = nc_amsr.variables['wind_speed_nwp'].long_name
    nc_ws_nwp.standard_name = nc_amsr.variables['wind_speed_nwp'].standard_name
    nc_ws_nwp.units = nc_amsr.variables['wind_speed_nwp'].units
    nc_ws_nwp.add_offset = nc_amsr.variables['wind_speed_nwp'].add_offset#ws_nwp_add_offset
    nc_ws_nwp.scale_factor = nc_amsr.variables['wind_speed_nwp'].scale_factor#ws_nwp_scale_factor
    nc_ws_nwp.valid_min = nc_amsr.variables['wind_speed_nwp'].valid_min
    nc_ws_nwp.valid_max = nc_amsr.variables['wind_speed_nwp'].valid_max
    nc_ws_nwp.comment = nc_amsr.variables['wind_speed_nwp'].comment
    nc_ws_nwp.references = nc_amsr.variables['wind_speed_nwp'].references
    nc_ws_nwp.source = nc_amsr.variables['wind_speed_nwp'].source
    nc_ws_nwp.time_offset = nc_amsr.variables['wind_speed_nwp'].time_offset
    nc_ws_nwp.height = nc_amsr.variables['wind_speed_nwp'].height
    nc_ws_nwp.coordinates = nc_amsr.variables['wind_speed_nwp'].coordinates
    nc_ws_nwp[:] = ws_nwp
    
    nc_out.close()
    

# coding: utf-8

# In[3]:



from netCDF4 import Dataset  
import datetime as dt
from datetime import datetime
import numpy as np
import math
import os
import sys
import pandas as pd
import xarray as xr
####################you will need to change some paths here!#####################
dir_mmdb='F:/data/mmd/mmd06c_re01_pp/drifter-sst_amsre-aq/' #'f:/data/mmd/new_mmd/';
dir_mmdb_out='f:/data/mmd/mmd06c_re01_pp/wind/' #'f:/data/mmd/mmd06c_wnd/';
dir_ccmp='F:/data/sat_data/ccmp/v02.0/Y' #'F:/data/sat_data/ccmp/v02.0/Y'
#################################################################################


input_year=int(str(sys.argv[1]))
input_month=int(str(sys.argv[2]))
#input_date=datetime(input_year,input_month,1)
#input_year=2003
#input_month=1
print('processing year and month below, month 13 means process all for that year')
print(input_year)
print(input_month)
#print input_date
#print input_date.month


# In[61]:


istart_flag = 0 
#filename_mmdb=dir_mmdb + 'mmd06c_sst_drifter-sst_amsre-aq_2002-152_2002-158.nc'
#filename_mmdb_out=dir_mmdb_out + 'mmd06c_sst_drifter-sst_amsre-aq_2002-152_2002-158.nc'
#import time

for root, dirs, files in os.walk(dir_mmdb, topdown=False):
    for name in files:
        if name.endswith(".nc"):
            filename_mmdb=os.path.join(root, name)
            mmdb_year=int(name[32:36])
            mmdb_jdy=int(name[37:40])
            mmdb_date=dt.datetime(mmdb_year, 1, 1) + dt.timedelta(mmdb_jdy - 1)
            mmdb_month=mmdb_date.month
            if mmdb_year!=input_year:
                continue
            if input_month<=12:
                if mmdb_month!=input_month:
                    continue
            filename_mmdb_out = dir_mmdb_out + 'ccmp_' + name
            fnc = Dataset(filename_mmdb,'r') 
            amsr_time=fnc.variables['amsre.time'][:,10,10]
            lats_amsr=fnc.variables['amsre.latitude'][:,10,10] 
            lons_amsr=fnc.variables['amsre.longitude'][:,10,10] 
            fnc.close()
            lons_amsr_360 = lons_amsr % 360  #put lons_amsr 0 to 360 for matching with ccmp which is 0-360
            ilen=len(amsr_time)
            print(filename_mmdb,' with ',ilen,' obs')
            date_amsre=pd.to_datetime(amsr_time.data, unit='s',origin='1993-01-01')
#            date_amsre=[0]*ilen
#            for i in range(0,ilen):
#                date_amsre[i]=dt.datetime(1993,1,1,0,0,0)+dt.timedelta(seconds=amsr_time[i])
#            rtime=amsr_time.tolist()
#            iarg=np.argsort(rtime)
            idysv=0
            istart=0
            col_wndu=[0]*ilen
            col_wndv=[0]*ilen
            date_start=min(date_amsre)
            if date_start.hour<=6:
                date_start=date_start-dt.timedelta(days=1)
            date_end=max(date_amsre)
            if date_end.hour>=18:
                date_end=date_end+dt.timedelta(days=1)
            date_start = dt.datetime(date_start.year,date_start.month,date_start.day)  #take off hours
            date_end = dt.datetime(date_end.year,date_end.month,date_end.day)  #take off hours
            date_incr=(date_end-date_start).days
            istart_flag=0
#initialize data, get CCMP winds 1 day plus/minus around day you are collocating in time
            for incr in range(0,date_incr+1):
                amsr_date2=date_start+dt.timedelta(days=incr)  
                syr, smon, sdym =str(amsr_date2.year).zfill(4), str(amsr_date2.month).zfill(2), str(amsr_date2.day).zfill(2)
                fname_tem='/CCMP_Wind_Analysis_' + syr + smon + sdym + '_V02.0_L3.0_RSS.nc'
                filename_ccmp = dir_ccmp + syr + '/M' + smon + fname_tem      
                print(filename_ccmp)
                ds = xr.open_dataset(filename_ccmp)
                if istart_flag == 0:
                    ds_ccmp = ds
                    istart_flag = 1
                    continue
                ds_ccmp = xr.concat([ds_ccmp, ds],'time')
#collocate data 
#            ccmp_data=ds_ccmp.interp(time=date_amsre[0:100],longitude=lons_amsr_360[0:100],latitude=lats_amsr[0:100],assume_sorted=False)                

#            start = time.time()
#            for i in range(0,100): #ilen):
#                if int(i/500)*500==i:
#                    print(ilen,i)
#                tem=ds_ccmp.interp(time=date_amsre[i],longitude=lons_amsr_360[i],latitude=lats_amsr[i]) #,assume_sorted=False)                
#                col_wndu[i]=tem.uwnd.data
#                col_wndv[i]=tem.vwnd.data
#            end = time.time()
#            print(end-start)
#            ccmp_data.to_netcdf(filename_mmdb_out)

#test interp function
#            start = time.time()
            t = xr.DataArray(date_amsre, dims='z')
            x = xr.DataArray(lons_amsr_360, dims='z')
            y = xr.DataArray(lats_amsr, dims='z')
            ccmp_data=ds_ccmp.interp(time=t,longitude=x,latitude=y)
            ccmp_data.to_netcdf(filename_mmdb_out)
#            end = time.time()
#            print(end-start)

#            break
#            fnc = Dataset(filename_mmdb_out,'w', format='NETCDF4') 
#            fnc.createDimension('t', ilen)
#            u_netcdf = fnc.createVariable('wndu', 'f4', ('t'))
#            v_netcdf = fnc.createVariable('wndv', 'f4', ('t'))
#            u_netcdf[:] = col_wndu
#            v_netcdf[:] = col_wndv
#            fnc.close()
            print('output data',filename_mmdb_out)



# In[50]:





# In[53]:





# In[58]:





# In[59]:






# coding: utf-8

# In[1]:


from netCDF4 import Dataset  
####################you will need to change some paths here!#####################
dir_mmdb='/group_workspaces/cems2/esacci_sst/mms_new/mmd/mmd06c_re01_pp/drifter-sst_amsre-aq/' #'f:/data/mmd/new_mmd/';
dir_mmdb_out='/group_workspaces/cems2/esacci_sst/mms_new/mmd/mmd06c_re01_pp/drifter-sst_amsre-aq/' #'f:/data/mmd/mmd06c_wnd/';
dir_ccmp='/group_workspaces/cems2/esacci_sst/scratch/jhoeyer001/CCMP_winds/ccmp/v02.0/Y2010/' #'F:/data/sat_data/ccmp/v02.0/Y'
#################################################################################
import datetime as dt
from datetime import datetime
import numpy as np
import math
import os


# In[2]:


istart_flag = 0 
for root, dirs, files in os.walk(dir_mmdb, topdown=False):
    for name in files:
        if name.endswith(".nc"):
            filename_mmdb=os.path.join(root, name)
            filename_mmdb_out = dir_mmdb_out + 'ccmp_' + name
            fnc = Dataset(filename_mmdb,'r') 
            amsr_time=fnc.variables['amsre.time'][:,10,10]
            lats_amsr=fnc.variables['amsre.latitude'][:,10,10] 
            lons_amsr=fnc.variables['amsre.longitude'][:,10,10] 
            fnc.close()
            ilen=len(amsr_time)
            print(filename_mmdb,' with ',ilen,' obs')
            date_amsre=[0]*ilen
            for i in range(0,ilen):
                date_amsre[i]=dt.datetime(1993,1,1,0,0,0)+dt.timedelta(seconds=amsr_time[i])
            rtime=amsr_time.tolist()
            iarg=np.argsort(rtime)
            idysv=0
            istart=0
            col_wndu=[0]*ilen
            col_wndv=[0]*ilen
            for i in range(0,ilen):
                amsr_date=date_amsre[iarg[i]]      
#initialize data, get CCMP winds 1 day plus/minus around day you are collocating in time
                if istart==0: 
                    for incr in range(-1,2):
                        amsr_date2=amsr_date+dt.timedelta(days=float(incr))  
                        syr=str(amsr_date2.year).zfill(4)
                        smon=str(amsr_date2.month).zfill(2)
                        sdym=str(amsr_date2.day).zfill(2)
                        sjdy=str(amsr_date2.timetuple().tm_yday).zfill(3)
                        fname_tem='/CCMP_Wind_Analysis_' + syr + smon + sdym + '_V02.0_L3.0_RSS.nc'
                        ccmp_filename = dir_ccmp + syr + '/M' + smon + fname_tem      
                        print(ccmp_filename)
                        nc_fid = Dataset(ccmp_filename, 'r')
                        wndu = nc_fid.variables['uwnd'][:]
                        wndv = nc_fid.variables['vwnd'][:]
                        mlat_ccmp = nc_fid.variables['latitude'][:]
                        mlon_ccmp = nc_fid.variables['longitude'][:]
                        t=nc_fid.variables['time'][:] 
                        time_ccmp=[0]*4
                        for itt in range(0,4):
                            time_ccmp[itt]=dt.datetime(1987,1,1,0,0,0)+dt.timedelta(hours=t[itt])
                        nc_fid.close()
                        if incr==-1:
                            wndu2=wndu
                            wndv2=wndv
                            time_ccmp2=time_ccmp
                        else:
                            wndu2 =  np.append(wndu2,wndu, axis=0)
                            wndv2 =  np.append(wndv2,wndv, axis=0)
                            time_ccmp2 = np.append(time_ccmp2,time_ccmp, axis = 0)
                        idysv=amsr_date.day
                        istart=1
                if amsr_date.day!=idysv:
                    amsr_date2=amsr_date+dt.timedelta(days=float(1))  #create new time array that can be queried for year etc
                    syr=str(amsr_date2.year).zfill(4)
                    smon=str(amsr_date2.month).zfill(2)
                    sdym=str(amsr_date2.day).zfill(2)
                    sjdy=str(amsr_date2.timetuple().tm_yday).zfill(3)
                    fname_tem='/CCMP_Wind_Analysis_' + syr + smon + sdym + '_V02.0_L3.0_RSS.nc'
                    ccmp_filename = dir_ccmp + syr + '/M' + smon + fname_tem      
                    print(ccmp_filename)
                    nc_fid = Dataset(ccmp_filename, 'r')
                    wndu = nc_fid.variables['uwnd'][:]
                    wndv = nc_fid.variables['vwnd'][:]
                    t=nc_fid.variables['time'][:] #units: hours since 1987-01-01 00:00:00
                    time_ccmp=[0]*4
                    for itt in range(0,4):
                        time_ccmp[itt]=dt.datetime(1987,1,1,0,0,0)+dt.timedelta(hours=t[itt])
                    nc_fid.close()
                    idysv=amsr_date.day
                    wndu2[0:8,:,:]=wndu2[4:12,:,:]
                    wndv2[0:8,:,:]=wndv2[4:12,:,:]
                    time_ccmp2[0:8]=time_ccmp2[4:12]
                    wndu2[8:12,:,:]=wndu[:]
                    wndv2[8:12,:,:]=wndv[:]
                    time_ccmp2[8:12]=time_ccmp[:]
                alat=lats_amsr[iarg[i]]
                alon=lons_amsr[iarg[i]]
                if alon<0:
                    alon=alon+360
                latli = np.argmin( np.abs( mlat_ccmp - alat ) )
                lonli = np.argmin( np.abs( mlon_ccmp - alon ) )
                timei = np.argmin( np.abs( time_ccmp2 - amsr_date ) )
                dttime=abs(time_ccmp2[timei] - amsr_date)
                f2=dttime.seconds/(6.*60*60)
                f1=abs(f2-1.)
                if time_ccmp2[timei]<amsr_date:
                    timej=timei+1
                if time_ccmp2[timei]>=amsr_date:
                    timej=timei-1
                tem=f1*wndu2[timei,latli,lonli]+f2*wndu2[timej,latli,lonli]
                col_wndu[iarg[i]]=tem
                col_wndv[iarg[i]]=f1*wndv2[timei,latli,lonli]+f2*wndv2[timej,latli,lonli]
            fnc = Dataset(filename_mmdb_out,'w', format='NETCDF4') 
            fnc.createDimension('t', ilen)
            u_netcdf = fnc.createVariable('wndu', 'f4', ('t'))
            v_netcdf = fnc.createVariable('wndv', 'f4', ('t'))
            u_netcdf[:] = col_wndu
            v_netcdf[:] = col_wndv
            fnc.close()
            print('output data',filename_mmdb_out)


# In[3]:


#get index order to process data in date ordered time
#print(min(date_amsre),max(date_amsre))
#print(date_amsre[iarg[0]],date_amsre[iarg[-1]])


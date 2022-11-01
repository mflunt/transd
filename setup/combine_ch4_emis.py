#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 15:32:09 2022

Combine wetland, EDGAR and GFED emissions

Also read in Chris' LST dataset

Will need to regrid all onto common 0.25x0.3125 grid.

Versions:
    1. EDGAR + wetcharts + gfed
    2. As 1 but replace wetcharts in 10x10 domain with Chris' file
    
Save at 0.25x0.3125 but also regrid to same scale as basis functions

@author: mlunt
"""

import numpy as np
import xarray 
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colorbar as cbar
import pandas as pd
import glob
from acrg_grid import regrid as acrg_regrid
import areagrid

def open_ds(fname):
    """
    Open netcdf file as xarray dataset. 
    Loads to memor then closes the file.
    """
    with xarray.open_dataset(fname) as ds:
        ds.load()
    return ds

def filenames(file_dir, file_str, start=None, end=None, freq='D'):
    """
    Output a list of available file names,for given directory and date range.
    """
    files = []
    # Convert into time format
    if (start != None) and (end != None):
        #days = pandas.DatetimeIndex(start = start, end = end, freq = "D").to_pydatetime()
        dates = pd.date_range(start = start, end = end, freq = freq)
        
        if freq =="D":
            yearmonthday = [str(d.year) + str(d.month).zfill(2) + str(d.day).zfill(2) for d in dates]
        elif freq =="MS":
            yearmonthday = [str(d.year) + str(d.month).zfill(2) for d in dates]
        elif freq =="YS":
            yearmonthday = [str(d.year) for d in dates]
    
        for ymd in yearmonthday:
            f=glob.glob(file_dir + "/" + file_str + "*" + ymd + "*")
            if len(f) > 0:
                files += f
        files.sort()
        
    else:
        f=glob.glob(file_dir + "/" + file_str + "*.nc")
        if len(f) > 0:
            files += f     # Not entirely sure if this will work - might be worth checking! - Yes it works
        files.sort()

    if len(files) == 0:
        raise ValueError("Can't find file: " + file_dir + "/" + file_str + "*")
                        
    return files

def read_reduced_netcdfs(files, lonmin,lonmax, latmin,latmax, dim = "time"):
    '''
    Use xray to open sequential netCDF files. 
    Makes sure that file is closed after open_dataset call.
    '''
    datasets = [open_ds(p) for p in sorted(files)]
    
    datasets_red=[]
    for ds in datasets:
        ds_dum = ds.sel(lon=slice(lonmin,lonmax), lat=slice(latmin,latmax))
        
        datasets_red.append(ds_dum)
    combined = xarray.concat(datasets_red, dim)
    return combined  

def read_reduced_single_nc(file, lonmin,lonmax, latmin,latmax):
    '''
    Use xray to open sequential netCDF files. 
    Makes sure that file is closed after open_dataset call.
    '''
    ds = open_ds(file)
    ds_red = ds.sel(lon=slice(lonmin,lonmax), lat=slice(latmin,latmax))
        
    return ds_red  

#%%
# Define data directories and file strings
    
wet_dir = "/home/mlunt/ceph/model_input_data/JPL_WetCHARTs/"
wet_str = "WetCHARTs_mean_v132_05x05_"   

edgar_dir = "/home/mlunt/ceph/model_input_data/EDGAR/"
edgar_str = "v6.0_CH4_"     #v6.0_CH4_2010_TOTALS.0.1x0.1.nc

gfed_dir =  "/home/mlunt/datastore/GFED/"
gfed_str = "GFED4s_CO2_CH4_025x025_monthly_"    #2010.nc

chris_dir = "/home/mlunt/datastore/EAF/emissions/"
#chris_fname = chris_dir + "ch4_lst_wetland_flux_sudd_025x025.nc"
#chris_fname = chris_dir + "ch4_lst_wetland_flux_sudd_v2_025x025.nc"
#chris_fname = chris_dir + "ch4_lst_wetland_flux_sudd_v2_Cs_025x025.nc"
#chris_fname = chris_dir + "ch4_lst_wetland_flux_EAfrica_v2_Cs_025x025.nc"

chris_fname = chris_dir + "ch4_lst_wetland_flux_merged_v2_Cs_030cm_025x025.nc"
#chris_fname = chris_dir + "ch4_lst_wetland_flux_merged_v2_Cs_050cm_025x025.nc"
#chris_fname = chris_dir + "ch4_lst_wetland_flux_merged_v2_Cs_100cm_025x025.nc"

#fname_out = chris_dir + "ch4_emis_eaf_025x03125_2010_2021.nc"
#fname_out = chris_dir + "ch4_emis_eaf_v2_025x03125_2010_2021.nc"
#fname_out = chris_dir + "ch4_emis_eaf_v2_Cs_025x03125_2010_2021.nc"

#fname_out = chris_dir + "ch4_emis_eaf_v2b Cs_025x03125_2010_2021.nc"

fname_out = chris_dir + "ch4_emis_eaf_v2c_Cs_030cm_025x03125_2010_2021.nc"
#fname_out = chris_dir + "ch4_emis_eaf_v2c_Cs_050cm_025x03125_2010_2021.nc"
#fname_out = chris_dir + "ch4_emis_eaf_v2c_Cs_100cm_025x03125_2010_2021.nc"

#%%
# Define bounds of regional EAF domain (can alsways cut down after)

latmin= -10
latmax= 20
lonmin = 20
lonmax = 55

start_date = "20100101"
end_date = "20211231"

fnames_wet = filenames(wet_dir, wet_str, start=start_date, end=end_date, freq='YS')
fnames_edgar = filenames(edgar_dir, edgar_str, start=start_date, end=end_date, freq='YS')
fnames_gfed = filenames(gfed_dir, gfed_str, start=start_date, end=end_date, freq='YS')

# Takes forever to do all years. GFED was daily not monthly which caused the issue

ds_wet = read_reduced_netcdfs(fnames_wet, lonmin,lonmax, latmin,latmax) # 2010-2020
ds_edgar = read_reduced_netcdfs(fnames_edgar, lonmin,lonmax, latmin,latmax)   # 2010-2018
ds_gfed = read_reduced_netcdfs(fnames_gfed, lonmin,lonmax, latmin,latmax)  # 2010-2021

ds_chris = open_ds(chris_fname)
#ds_wet = read_reduced_single_nc(fnames_wet[0], lonmin,lonmax, latmin,latmax)
#ds_edgar = read_reduced_single_nc(fnames_edgar[0], lonmin,lonmax, latmin,latmax)
#ds_gfed = read_reduced_single_nc(fnames_gfed[0], lonmin,lonmax, latmin,latmax)

emis_wet = ds_wet.emis_ch4.fillna(0.) /1.e6/24./60/60   # converted from mgch4/m2/day
emis_edgar = ds_edgar.emi_ch4     # kg/m2/s
emis_gfed = ds_gfed.ch4_emissions  # kg/m2/s
emis_chris0 = (ds_chris.flux0.sel(time=slice(start_date,end_date))).fillna(0.)
emis_chris2 = (ds_chris.flux2.sel(time=slice(start_date,end_date))).fillna(0.)
emis_chris4 = (ds_chris.flux4.sel(time=slice(start_date,end_date))).fillna(0.)
emis_chris6 = (ds_chris.flux6.sel(time=slice(start_date,end_date))).fillna(0.)

dlat=0.25
dlon=0.3125

lat_out = np.arange(-8,18 + dlat, dlat)
lon_out = np.arange(22.5, 52.5+dlon,dlon)

nlat_out = len(lat_out)
nlon_out=len(lon_out)
#%%
# Regrid each emis field to common grid. 

lat_edgar = emis_edgar.lat.values
lon_edgar = emis_edgar.lon.values
time_edgar = emis_edgar.time.values
ntime_edgar = len(time_edgar)

q_edgar_new = np.zeros((ntime_edgar,nlat_out,nlon_out))
for ti in range(ntime_edgar):
    q_edgar_new[ti,:,:],dum = acrg_regrid.regrid2d(emis_edgar[ti,:,:].values, 
             lat_edgar,lon_edgar,lat_out,lon_out)

lat_wet = emis_wet.lat.values
lon_wet = emis_wet.lon.values
time_wet = emis_wet.time.values
ntime_wet = len(time_wet)

q_wet_new = np.zeros((ntime_wet,nlat_out,nlon_out))
for ti in range(ntime_wet):
    q_wet_new[ti,:,:],dum = acrg_regrid.regrid2d(emis_wet[ti,:,:].values, 
             lat_wet,lon_wet,lat_out,lon_out)
    
lat_gfed = emis_gfed.lat.values
lon_gfed = emis_gfed.lon.values
time_gfed = emis_gfed.time.values
ntime_gfed = len(time_gfed)

q_gfed_new = np.zeros((ntime_gfed,nlat_out,nlon_out))
for ti in range(ntime_gfed):
    q_gfed_new[ti,:,:],dum = acrg_regrid.regrid2d(emis_gfed[ti,:,:].values, 
             lat_gfed,lon_gfed,lat_out,lon_out)

lat_chris = emis_chris0.lat.values
lon_chris = emis_chris0.lon.values
time_chris = emis_chris0.time.values
ntime_chris = len(time_chris)

q_chris_new0 = np.zeros((ntime_chris,nlat_out,nlon_out))
q_chris_new2 = np.zeros((ntime_chris,nlat_out,nlon_out))
q_chris_new4 = np.zeros((ntime_chris,nlat_out,nlon_out))
q_chris_new6 = np.zeros((ntime_chris,nlat_out,nlon_out))
for ti in range(ntime_chris):
    q_chris_new0[ti,:,:],dum = acrg_regrid.regrid2d(emis_chris0[ti,:,:].values, 
             lat_chris,lon_chris,lat_out,lon_out)
    
    q_chris_new2[ti,:,:],dum = acrg_regrid.regrid2d(emis_chris2[ti,:,:].values, 
             lat_chris,lon_chris,lat_out,lon_out)
    
    q_chris_new4[ti,:,:],dum = acrg_regrid.regrid2d(emis_chris4[ti,:,:].values, 
             lat_chris,lon_chris,lat_out,lon_out)
    
    q_chris_new6[ti,:,:],dum = acrg_regrid.regrid2d(emis_chris6[ti,:,:].values, 
             lat_chris,lon_chris,lat_out,lon_out)

#%%
# Combine emissions. All fields monthly except EDGAR. 
# But need to account for those that don't cover later years
    
# Firstly put everythin onto common time grid - i.e. EDGAR monthly and wetcharts thru 2021
dates_out = pd.date_range(start=start_date, end=end_date, freq="MS")
ntime_out = len(dates_out)    
    
emis_edgar_new2 = np.zeros((ntime_out, nlat_out,nlon_out))

for ti in range(ntime_out):
    if dates_out[ti].year < 2019:
        ti_ed = ti//12
        emis_edgar_new2[ti,:,:] = q_edgar_new[ti_ed,:,:]
        
    else:
        emis_edgar_new2[ti,:,:] = q_edgar_new[-1,:,:]

emis_wet_new2 = np.zeros((ntime_out, nlat_out,nlon_out))

base=0
for ti in range(ntime_out):
    if dates_out[ti].year < 2021:
        emis_wet_new2[ti,:,:] = q_wet_new[ti,:,:]
    else:
        q_wet_clim = np.mean(q_wet_new[base::12,:,:], axis=0)  
        emis_wet_new2[ti,:,:] = q_wet_clim
        base+=1


emis_total = emis_edgar_new2 + emis_wet_new2 + q_gfed_new

da_emis = xarray.DataArray(data= emis_total, dims=["time", "lat", "lon"], 
                           coords = {"time": dates_out,
                                     "lat": lat_out,
                                     "lon": lon_out})
    
da_chris_out = xarray.DataArray(data= q_chris_new0, dims=["time", "lat", "lon"], 
                           coords = {"time": dates_out,
                                     "lat": lat_out,
                                     "lon": lon_out})
    
# Replace wetcharts with CHris's field between coords
    
lonmin_chris = lon_chris.min()
lonmax_chris = lon_chris.max()
latmin_chris = lat_chris.min()
latmax_chris = lat_chris.max()

wh_lon0 = np.where(np.abs(lon_out-lonmin_chris) <= dlon/2)[0][0]
wh_lon1 = np.where(np.abs(lon_out-lonmax_chris) <= dlon/2)[0][0]

wh_lat0 = np.where(np.abs(lat_out-latmin_chris) <= dlon/2)[0][0]
wh_lat1 = np.where(np.abs(lat_out-latmax_chris) <= dlon/2)[0][0]

da_emis2 = da_emis.copy()
da_emis2[:,wh_lat0:wh_lat1+1, wh_lon0:wh_lon1+1] = (
        da_emis[:,wh_lat0:wh_lat1+1, wh_lon0:wh_lon1+1] 
        - emis_wet_new2[:,wh_lat0:wh_lat1+1, wh_lon0:wh_lon1+1] 
        + q_chris_new0[:,wh_lat0:wh_lat1+1, wh_lon0:wh_lon1+1] )
        
# Create dataset with totals and individual sectors in. Save to file.
# Write a function within inversion script to calculate prior emissions.

ds_out= xarray.Dataset()
ds_out["ch4_total"] = da_emis
ds_out["ch4_total_lst0"] = da_emis2
ds_out["ch4_total_lst2"] = da_emis2 - q_chris_new0 + q_chris_new2
ds_out["ch4_total_lst4"] = da_emis2 - q_chris_new0 + q_chris_new4
ds_out["ch4_total_lst6"] = da_emis2 - q_chris_new0 + q_chris_new6

ds_out["ch4_wetcharts"] = (["time", "lat", "lon"], emis_wet_new2)
ds_out["ch4_edgar"] = (["time", "lat", "lon"],emis_edgar_new2)
ds_out["ch4_gfed"] = (["time", "lat", "lon"],q_gfed_new)
ds_out["ch4_lst0"] = (["time", "lat", "lon"],q_chris_new0)
ds_out["ch4_lst2"] = (["time", "lat", "lon"],q_chris_new2)
ds_out["ch4_lst4"] = (["time", "lat", "lon"],q_chris_new4)
ds_out["ch4_lst6"] = (["time", "lat", "lon"],q_chris_new6)

for key in list(ds_out.keys()):
    
    ds_out[key].attrs["units"] = "kg/m2/s"
    
ds_out.attrs["comments"] = "CH4 emissions from WetCHARTs v1.3.2, EDGAR v6, GFEDv4 and Chris Taylors LST approach"
 
ds_out.to_netcdf(path=fname_out)

#%%
cmin=0
cmax=1.e-9
proj = ccrs.PlateCarree()
fig,axs = plt.subplots(2, figsize=(8,8), subplot_kw=dict(projection=proj))

p2=axs[0].pcolormesh(lon_out-dlon/2, lat_out-dlat/2, da_emis.mean(dim="time"), 
            transform=ccrs.PlateCarree(), cmap='viridis', vmin=cmin, vmax=cmax)

axs[0].coastlines()
axs[0].add_feature(cfeature.BORDERS, edgecolor='black', alpha=0.8)

p2=axs[1].pcolormesh(lon_out-dlon/2, lat_out-dlat/2, da_emis2.mean(dim="time")-da_emis.mean(dim="time"), 
            transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=-1.e-9, vmax=1.e-9)

axs[1].coastlines()
axs[1].add_feature(cfeature.BORDERS, edgecolor='black', alpha=0.8)

cbaxes2 = fig.add_axes([0.1, 0.1, 0.8, 0.04]) 
##[left, bottom, width, height],
cb = plt.colorbar(p2, cax = cbaxes2, orientation='horizontal', extend='both', label = 'Emissions')  


sudd_total = da_emis.sel(lon=slice(lonmin_chris,lonmax_chris), lat=slice(latmin_chris,latmax_chris)).mean(dim=["lat","lon"])
sudd_total2 = da_emis2.sel(lon=slice(lonmin_chris,lonmax_chris), lat=slice(latmin_chris,latmax_chris)).mean(dim=["lat","lon"])

sudd_total_wet = ds_out.ch4_wetcharts.sel(lon=slice(lonmin_chris,lonmax_chris), lat=slice(latmin_chris,latmax_chris)).mean(dim=["lat","lon"])
sudd_total_edg = ds_out.ch4_edgar.sel(lon=slice(lonmin_chris,lonmax_chris), lat=slice(latmin_chris,latmax_chris)).mean(dim=["lat","lon"])
sudd_total_gfed = ds_out.ch4_gfed.sel(lon=slice(lonmin_chris,lonmax_chris), lat=slice(latmin_chris,latmax_chris)).mean(dim=["lat","lon"])

#lat_sudd = sudd_total.lat.values
#lon_sudd = sudd_total.lon.values

#area = areagrid.areagrid(lat_sudd,lon_sudd)
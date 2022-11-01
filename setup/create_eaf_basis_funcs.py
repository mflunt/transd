#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 17 09:12:15 2021

Create E. Africa basis functions for inversion

Will just run with a flat prior so I have the exact dy/dx and prior dist doesn't matter.

In S. Sudan need 0.5x0.625
In E. Africa (Uganda, Kenya, w. Ethiopia) 1x1.25 
For rest use 2x2.5

Maybe just do 1x1.25 over regions of interest in S. Sudan, Uganda, Ethiopia.
Outside of that use 2x2.5.

Could run as 2 separate runs to save space. 120 tracers for inner region.
Think 2x2.5 might leave me with 156-30 so 126 for outer region. But can ignore any not on land.
Which might save another 1/5th. So maybe 100 2x2.5 regions? 

@author: mlunt
"""

import xarray
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import regionmask
import pandas as pd

latmin = -8
latmax= 18
lonmin = 22.5
lonmax = 52.5
#dlon = 0.625
#dlat = 0.5

dlon = 0.3125
dlat = 0.25

lon = np.arange(lonmin,lonmax+dlon,dlon)
lat = np.arange(latmin,latmax+dlat,dlat)
nlon=len(lon)
nlat=len(lat)



lonmin_sud = 27.5
lonmax_sud = 40

latmin_sud = 0
latmax_sud = 12

#nlon_sud = int((lonmax_sud-lonmin_sud)/dlon/2)
#nlat_sud = int((latmax_sud-latmin_sud)/dlat/2)

nlon_sud = int((lonmax_sud-lonmin_sud)/dlon/4)
nlat_sud = int((latmax_sud-latmin_sud)/dlat/4)

mask = regionmask.defined_regions.natural_earth.countries_110.mask(lon, lat, xarray=True)

lonmin_eaf = 30
lonmax_eaf = 45

#latmin_eaf = 

emis_regions_inner = np.zeros((nlat,nlon))
emis_regions_outer = np.zeros((nlat,nlon))

#nlon25 = int((nlon-1)/4)
#nlat2 = int((nlat-1)/4)

nlon25 = int((nlon-1)/8)
nlat2 = int((nlat-1)/8)

wh_sud_lon0 = np.where(lon == lonmin_sud)[0][0]
wh_sud_lat0 = np.where(lat == latmin_sud)[0][0]

count_inner = 1
for lai in range(nlat_sud):
    for loi in range(nlon_sud):
        #emis_regions_inner[wh_sud_lat0+lai*2: wh_sud_lat0+(lai+1)*2, wh_sud_lon0+loi*2: wh_sud_lon0+(loi+1)*2] = count_inner
        emis_regions_inner[wh_sud_lat0+lai*4: wh_sud_lat0+(lai+1)*4, wh_sud_lon0+loi*4: wh_sud_lon0+(loi+1)*4] = count_inner
        count_inner+=1

count_outer=1
for lai in range(nlat2):
    for loi in range(nlon25):
        
        wh_fin = np.where(np.ravel(np.isfinite(mask[lai*8:(lai+1)*8, loi*8:(loi+1)*8].values)))
        
        if (lon[loi*8] >= lonmin_sud and lon[loi*8] < lonmax_sud
            and lat[lai*8] >= latmin_sud and lat[lai*8] < latmax_sud): 
        
            a=1
        #elif np.isfinite(mask[lai*8, loi*8].values) == False: 
        elif len(wh_fin[0]) <4:
        #elif np.sum(mask[lai*8:(lai+1)*8, loi*8:(loi+1)*8]) <1: 
            a=2
        else:
            emis_regions_outer[lai*8:(lai+1)*8, loi*8:(loi+1)*8] = count_outer
            count_outer+=1
            
#        if (lon[loi*4] >= lonmin_sud and lon[loi*4] < lonmax_sud
#            and lat[lai*4] >= latmin_sud and lat[lai*4] < latmax_sud): 
#        
#            a=1
#        elif np.isfinite(mask[lai*4, loi*4].values) == False:     
#            a=2
#        else:
#            emis_regions_outer[lai*4:(lai+1)*4, loi*4:(loi+1)*4] = count_outer
#            count_outer+=1
#%%
emis_inner2 = emis_regions_inner.copy()
wh_inner = np.where(emis_inner2 >0)
emis_inner2[wh_inner] = emis_inner2[wh_inner]+count_outer-1
cmin=0
cmax=218     
proj = ccrs.PlateCarree()
fig2,ax2 = plt.subplots(figsize=(8,6), subplot_kw=dict(projection=proj))

p2=ax2.pcolormesh(lon, lat, emis_regions_outer+ emis_inner2,
            transform=ccrs.PlateCarree(), cmap='viridis', vmin=cmin, vmax=cmax)
ax2.coastlines()
ax2.add_feature(cfeature.BORDERS, edgecolor='black', alpha=0.8)
cbaxes2 = fig2.add_axes([0.1, 0.1, 0.8, 0.04]) 
##[left, bottom, width, height],
cb = plt.colorbar(p2, cax = cbaxes2, orientation='horizontal', extend='both', label = 'Region count')  

# Add some gridlines

#%%
# Create a different mask file for both sets. 
# Also create a flat emission field (at 0.25x0.3125) resolution
# Run monthly model runs x2. Should run quicket than EU runs due to smaller domain: 26*30 vs 32*50
# So should be 2x quicker. Set up 12 run scripts at once. Might get 1 year done by end of weekend?

ds_out_inner = xarray.Dataset()
for zi in range(count_inner-1):
    ds_out_inner["bf_" + str(zi+1).zfill(2)] = (("time", "lat", "lon"), np.zeros((1,nlat,nlon)))
    wh = np.where(emis_regions_inner == zi+1)
    basis_dum = emis_regions_inner.copy()*0.
    basis_dum[wh]= 1.
    ds_out_inner["bf_" + str(zi+1).zfill(2)][0,:,:] =  basis_dum

keys_inner = list(ds_out_inner.keys())
for reg in keys_inner:
    # Make sure edges of domain are 0.
    ds_out_inner[reg][0,:3,:]=0.
    ds_out_inner[reg][0,-3:,:] = 0.
    ds_out_inner[reg][0,:,:3]=0.
    ds_out_inner[reg][0,:,-3:]=0.
    ds_out_inner[reg].attrs["units"] = "1"
    ds_out_inner[reg].attrs["long_name"] = "Mask_for_" + reg + "_region"

ds_out_inner.coords["lat"] = lat
ds_out_inner.coords["lon"] = lon
ds_out_inner.coords["time"] = [pd.to_datetime("20190101")]

ds_out_inner["lat"].attrs["units"] = "degrees_north"
ds_out_inner["lat"].attrs["standard_name"] = "latitude"
ds_out_inner["lat"].attrs["long_name"] = "latitude" 
ds_out_inner["lat"].attrs["axis"] = "Y" 

ds_out_inner["lon"].attrs["units"] = "degrees_east"
ds_out_inner["lon"].attrs["standard_name"] = "longitude"
ds_out_inner["lon"].attrs["long_name"] = "longitude" 
ds_out_inner["lon"].attrs["axis"] = "X" 

out_dir = "/home/mlunt/ceph/verify/model_settings/DARE_runs/CH4/eaf_runs/masks/"
fname_out_inner = out_dir + "eaf_bfs_inner.nc"
#ds_out_inner.to_netcdf(path=fname_out_inner, mode="w")

#%%
ds_out_outer = xarray.Dataset()
for zi in range(count_outer-1):
    ds_out_outer["bf_" + str(zi+1).zfill(2)] = (("time", "lat", "lon"), np.zeros((1,nlat,nlon)))
    wh = np.where(emis_regions_outer == zi+1)
    basis_dum = emis_regions_outer.copy()*0.
    basis_dum[wh]= 1.
    ds_out_outer["bf_" + str(zi+1).zfill(2)][0,:,:] =  basis_dum

keys_outer = list(ds_out_outer.keys())
for reg in keys_outer:
    # Make sure edges of domain are 0.
    ds_out_outer[reg][0,:3,:]=0.
    ds_out_outer[reg][0,-3:,:] = 0.
    ds_out_outer[reg][0,:,:3]=0.
    ds_out_outer[reg][0,:,-3:]=0.
    ds_out_outer[reg].attrs["units"] = "1"
    ds_out_outer[reg].attrs["long_name"] = "Mask_for_" + reg + "_region"

ds_out_outer.coords["lat"] = lat
ds_out_outer.coords["lon"] = lon
ds_out_outer.coords["time"] = [pd.to_datetime("20190101")]

ds_out_outer["lat"].attrs["units"] = "degrees_north"
ds_out_outer["lat"].attrs["standard_name"] = "latitude"
ds_out_outer["lat"].attrs["long_name"] = "latitude" 
ds_out_outer["lat"].attrs["axis"] = "Y" 

ds_out_outer["lon"].attrs["units"] = "degrees_east"
ds_out_outer["lon"].attrs["standard_name"] = "longitude"
ds_out_outer["lon"].attrs["long_name"] = "longitude" 
ds_out_outer["lon"].attrs["axis"] = "X" 

fname_out_outer = out_dir + "eaf_bfs_outer.nc"
#ds_out_outer.to_netcdf(path=fname_out_outer, mode="w")

#%%
# Create flat emission file 90.25x0.3215)

lat_flux = np.arange(latmin,latmax+0.25,0.25)
lon_flux = np.arange(lonmin,lonmax+0.3125,0.3125)

cnt_hires = regionmask.defined_regions.natural_earth.countries_110.mask(lon_flux, lat_flux, xarray=True)
flux_out = (cnt_hires*0.+1.).fillna(0.)
flux_out2 = flux_out.expand_dims(dim="time")

# Need to put flux out into realistic units. Each grid cell say 0.01 Tg?
# Or assume a rate of 2 mg/m2/hour.  = 5.555e-10 kg/m2/s - so maybe use 5e-10 or 1e-9.
flux_out2 = flux_out2*1.e-9 # Now in kg/m2/s.

ds_out2 = xarray.Dataset()

#ds_out2["emis_ch4"] = (["time", "lat", "lon"], flux_out2)
ds_out2["emis_ch4"] =  flux_out2

ds_out2["emis_ch4"].attrs["units"] = "kg/m2/s"
ds_out2["emis_ch4"].attrs["standard_name"] = "emissions"
#ds_out2[varname].attrs["comments"] = "Flux rate of 2 mgCH4/m2/hour for outer regions, 4mgCH4/m2/hour for mid and 10 mgCH4/ms/hour for inner."
ds_out2["emis_ch4"].attrs["comments"] = "Flux rate of 1e-9 kg/m2/s for all land grid cells."

#ds_out2.coords["lat"] = lat_emi
#ds_out2.coords["lon"] = lon_emi
ds_out2.coords["time"] = [pd.to_datetime("20190101")]

ds_out2["lat"].attrs["units"] = "degrees_north"
ds_out2["lat"].attrs["standard_name"] = "latitude"
ds_out2["lat"].attrs["long_name"] = "latitude" 
ds_out2["lat"].attrs["axis"] = "Y" 

ds_out2["lon"].attrs["units"] = "degrees_east"
ds_out2["lon"].attrs["standard_name"] = "longitude"
ds_out2["lon"].attrs["long_name"] = "longitude" 
ds_out2["lon"].attrs["axis"] = "X" 

out_dir_emis = "/home/mlunt/ceph/verify/model_settings/DARE_runs/CH4/eaf_runs/emissions/"
fname_out_emis = out_dir_emis + "emis_ch4_flat_eaf_025x03125_2019.nc"
#ds_out2.to_netcdf(path=fname_out_emis, mode="w")


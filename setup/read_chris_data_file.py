#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 09:37:43 2022

Read in Chris' Grads binary file

@author: mlunt
"""


import numpy as np
import xarray
import matplotlib.pyplot as plt
import pandas as pd

data_dir = "/home/mlunt/datastore/EAF/emissions/"
#fname = data_dir + "CH4_flux_from_wetland_025_vn1_h20h21v08.gra"
#fname = data_dir + "CH4_flux_from_wetland_025_vn2_h20h21v08.gra"
#fname = data_dir + "CH4_flux_from_wetland_025_vn2_h20h21v08_Cs_var_010cm.gra"

#fname = data_dir + "CH4_flux_from_wetland_025_vn2_EAfrica_Cs_var_010cm.gra"

#fname = data_dir + "CH4_flux_from_wetland_025_vn2_merged_Cs_varP_030cm.gra"
#fname = data_dir + "CH4_flux_from_wetland_025_vn2_merged_Cs_varP_050cm.gra"
fname = data_dir + "CH4_flux_from_wetland_025_vn2_merged_Cs_varP_100cm.gra"

#file = open(fname, "rb")


#dtype = np.dtype('B')
dtype = np.dtype('<f')  # Little endian float type
try:
    with open(fname, "rb") as f:
        np_data = np.fromfile(f,dtype=dtype)
    print(np_data)
except IOError:
    print('Error While Opening the file!')    
    
    
cats = ["pc_lst", "flux0", "flux2", "flux4", "flux6" ]
ncat = 5

#nlon = 40
#nlat = 40
#ntime = 228
#
#lat = np.arange(0,10, 0.25)
#lon = np.arange(27.5, 37.5, 0.25)

nlon = 110
nlat = 80
ntime = 228

ntime = 252

lat = np.arange(-5.,15, 0.25)
lon = np.arange(25., 52.5, 0.25)

#dates = pd.date_range(start = "20030101", end="20211231", freq= "MS")
dates = pd.date_range(start = "20020101", end="20221231", freq= "MS")

ndata = len(np_data)

#%%
# Order  (fastest to slowest) = lon,lat, cat, time

data_list={}

for cat in cats:
    
    data_list[cat] = []

for ti in range(ntime):
    
    for ct in range(ncat):
    
        data_v = np_data[ti*ncat*nlon*nlat + ct*nlon*nlat: ti*ncat*nlon*nlat + (ct+1)*nlon*nlat]
        
        data_v[data_v <-900] = np.nan
        
        data_2d = np.reshape(data_v, (nlat,nlon))
        
        da_data_2d = xarray.DataArray(data= data_2d, dims=["lat", "lon"], coords={"lat": lat, "lon":lon})
        
        data_list[cats[ct]].append(da_data_2d)
        
# Concat along time dimensions
ds_out=xarray.Dataset()
for cat in cats:
    ds_out[cat] = xarray.concat(data_list[cat], dim="time")
    #ds_out[cat].coords["time"] = dates
    
ds_out.coords["time"] = dates
#%%
# Looks sensible. Now need to save to file. Then combine with
# other wetland fields (Wetcharts?) and EDGAR and GFED for prior. 
# Get everything on common 0.25 x0.3125 grid.
threshold={"flux0": "0 K", 
            "flux2": "-2 K",
            "flux4": "-4 K",
            "flux6": "-6 K",}
# Save as native file to begin with so I don't ruin anything.
    
for cat in cats[1:]:
    ds_out[cat].attrs["unit"] = "kg/m2/s"
    ds_out[cat].attrs["long_name"] = "ch4_wetland_flux" 
    ds_out[cat].attrs["comment"] = "ch4 wetland flux calculated with LST threshold of " + threshold[cat]
    
ds_out["pc_lst"].attrs["comment"] = "Fraction of pixels with valid LST in each month"

#fname_out = data_dir + "ch4_lst_wetland_flux_sudd_v2_025x025.nc"
fname_out = data_dir + "ch4_lst_wetland_flux_merged_v2_Cs_100cm_025x025.nc"
#ds_out.to_netcdf(path=fname_out)
    
    
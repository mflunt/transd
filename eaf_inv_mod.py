#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 12:45:20 2021

Module file containing functions for EnKF

Probably need to turn this inot a class one day

@author: mlunt
"""

import numpy as np
import xarray
#import scipy.linalg as sp_linalg
#from acrg_grid import haversine
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
#import matplotlib.colorbar as cbar
#import scipy.stats as spstats
#from dateutil.relativedelta import relativedelta
import glob
import pandas as pd
#from numba import jit


def open_ds(fname):
    """
    Open netcdf file as xarray dataset. 
    Loads to memor then closes the file.
    """
    with xarray.open_dataset(fname) as ds:
        ds.load()
    return ds

def filenames(start, end, file_dir, file_string):
    """
    Output a list of available file names,
    for given directory and date range.
    Assumes monthly files
    """
    
    # Convert into time format
    #days = pd.DatetimeIndex(start = start, end = end, freq = "1D").to_pydatetime()
    
    days = pd.date_range(start = start, end = end, freq = "1D")
    yrmnd = [str(d.year) + str(d.month).zfill(2) + str(d.day).zfill(2) for d in days]

    files = []
    for ymd in yrmnd:
        f=glob.glob(file_dir  + "/" + file_string +  
                    ymd + "*.nc")
        if len(f) > 0:
            files += f
    files.sort()

    if len(files) == 0:
        print("Can't find file: " + file_dir + "/" + file_string + ymd +  "*.nc")
                        
    return files

def read_netcdfs(files, dim = "time"):
    '''
    Use xray to open sequential netCDF files. 
    Makes sure that file is closed after open_dataset call.
    '''
    datasets = [open_ds(p) for p in sorted(files)]
    
    # Need to add time dimension to lat and lon
    #for ds in datasets:
    #    ds["lon"] = ds["lon"].expand_dims(dim="time")
    #    ds["lat"] = ds["lat"].expand_dims(dim="time")
    
    combined = xarray.concat(datasets, dim)
    return combined   

###@jit(nopython=True)
def numba_calc_R(deltaspace, deltatime, rho, tau, y_uncert):
        
#        Qinv_temp = np.zeros((nobs,nobs))
#        for ti in range(ndays):
#            
#            
#            Q_block = np.exp((-1.)*deltaspace[cum_nmeas2[ti]:cum_nmeas2[ti+1],
#                             cum_nmeas2[ti]:cum_nmeas2[ti+1]] / rho)
#            Q_block_inv = np.linalg.inv(Q_block)
#            Qinv_temp[cum_nmeas2[ti]:cum_nmeas2[ti+1],
#                      cum_nmeas2[ti]:cum_nmeas2[ti+1]] = Q_block_inv.copy()
#        
#        dum1 = np.dot(np.diag(1./y_uncert),Qinv_temp)
#        R_inv = np.dot(dum1,np.diag(1./y_uncert))
        
        #Qinv_temp = np.zeros((nobs,nobs))
       
    Q_block = np.exp((-1.)*deltaspace/rho) * np.exp((-1.)*deltatime/tau)
                       
    Q_block_inv = np.linalg.inv(Q_block)
       
    dum1 = np.dot(np.diag(1./y_uncert),Q_block_inv)
    R_inv = np.dot(dum1,np.diag(1./y_uncert))
        
    return R_inv
    

def custom_sqrtm(A):
    
    B,S = np.linalg.eigh(A)
    diag_B_sq = np.diag(np.sqrt(B))
    sqrt_A=np.dot(np.dot(S,diag_B_sq),np.linalg.inv(S))
    
    return sqrt_A

def read_sat_output(species, run_dir, run_str, run_date, 
                        start_date, end_date, keys, 
                        satellite = "GOSAT", read_obs=True,
                        bremen=False):
    y_mod_all_dict={}
    #for key in keys:
    #    y_mod_all_dict[key]=[]
    
        
        
    fnames = filenames(start_date, end_date, run_dir,  run_str)

    ds_site1 = read_netcdfs(fnames, dim="nobs")
    for key in keys:
        y_mod_all_dict[key] = np.ravel(ds_site1[key].values)
    


    y_mod_all_sites=[]
    for key in keys:
        if len(y_mod_all_dict[key])>0:
            y_mod_all_sites.append(np.hstack(y_mod_all_dict[key]))
        
    y_mod_all_out = np.stack(y_mod_all_sites)
    
    if read_obs == True:
        
        y_obs = np.ravel(ds_site1["XCH4_obs"].values)
        #y_std = np.ravel(ds_site1["XCH4_uncertainty"].values)
        y_lat = np.ravel(ds_site1["lat"].values)
        y_lon = np.ravel(ds_site1["lon"].values)
        y_date = np.ravel(ds_site1["obs_date"].values)
        
        if satellite == "GOSAT":
            y_std = np.ravel(ds_site1["XCH4_uncertainty"].values)
        elif satellite == "TROPOMI" and bremen == False:
            y_std = np.ravel(ds_site1["XCH4_precision"].values)
        elif satellite == "TROPOMI" and bremen == True:
            y_std = np.ravel(ds_site1["XCH4_uncertainty"].values)
            
        y_dict_out={"obs": y_obs,
                    "std": y_std,
                    "lat": y_lat,
                    "lon": y_lon,
                    "date":y_date}
        
        if satellite =="GOSAT":
            y_flag = np.ravel(ds_site1["retr_flag"].values)
            y_dict_out["flag"] = y_flag
        
        if satellite == "TROPOMI" and bremen ==True:
            
            albedo = np.ravel(ds_site1["albedo"].values)
            surf_rough = np.ravel(ds_site1["surface_roughness"].values)
            land_frac = np.ravel(ds_site1["land_fraction"].values)
            cloud_param = np.ravel(ds_site1["cloud_parameter"].values)
            
            y_dict_out["albedo"] = albedo
            y_dict_out["surf_rough"] = surf_rough
            y_dict_out["land_frac"] = land_frac
            y_dict_out["cloud_param"] = cloud_param
        elif satellite == "TROPOMI" and bremen == False:
            
            albedo = np.ravel(ds_site1["SWIR_albedo"].values)
            aot = np.ravel(ds_site1["SWIR_AOT"].values)
            asize = np.ravel(ds_site1["aerosol_size"].values)
            xco = np.ravel(ds_site1["XCO"].values)
            
            y_dict_out["albedo"] = albedo
            y_dict_out["AOT"] = aot
            y_dict_out["asize"] = asize
            y_dict_out["xco"] = xco
        
        return y_mod_all_out, y_dict_out
    else:
        return y_mod_all_out    

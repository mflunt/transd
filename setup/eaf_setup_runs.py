#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 16:40:23 2021

Write restart file for CO2 and CO inversions

@author: mlunt
"""

from enkf_code import GC_setup_mod as ensemble_mod
import matplotlib.pyplot as plt
import subprocess
import glob
import pandas as pd
from dateutil.relativedelta import relativedelta
import xarray
import numpy as np
import re

def write_input_file(start_date, end_date, 
                     run_dir, output_dir, template_dir, 
                      species, run_key, met_type, lag=False):
    """
    Write the input.geos file for each assimilation run. 
    Length of run needs to be for whole of lag period.
    
    Inputs:
        start_date (string): "YYYYMMDD" start of emissions window
        lag_start (string): "YYYYMMDD" Start of lag window/end of emissions window
        end_date (string): "YYYYMMDD" End of lag window
        N_ens (int): NUmber of ensemble members
        bc_dir (string): Directory where BC files are kept.
        spc_BC (list): list of boundary variables to include. - Not neededI don't think
        
    Requires:
        Template file - input_ensemble.template - to be setup correctly. 
    """
    if met_type == "MERRA2":
        template_file = template_dir + run_key + "_" + species + "_merra2_input.template"
    else:
        template_file = template_dir + run_key + "_" + species + "_input.template"
    # Read in template input file
    # Read in the file
    #with open(template_dir + "input.template", 'r') as file :
    #  filedata = file.read()
    ## Replace the target string
    
    if lag == True:
        input_file_name = run_dir + "input_files/lag_input.geos"
    else:
        input_file_name = run_dir + "input_files/window_input.geos"
     
    #with open(template_dir + "input_ensemble_CH4.template", "r") as in_file:
    with open(template_file, "r") as in_file:
        filedata = in_file.readlines()
    
    with open(input_file_name, "w") as out_file:
        for line in filedata:
            
            if "StartDate" in line: 
                line = line.replace('{StartDate}', start_date)
            if "EndDate" in line:
                line = line.replace('{EndDate}', end_date)          
           
            if "{outputdir}" in line:
                line = line.replace('{outputdir}', output_dir)     
                
            out_file.write(line)

    print("Successfully written file " + input_file_name)
    
    return

def write_mask_file(latmin,latmax,lonmin,lonmax, dlat,dlon, fname_out, start_date):
    """
    Write scale factor netcdf file for GC run. 
    
    2 timesteps - 1st is emissions period.
                - 2nd is lag period where scale factors should be 0.
                 -3rd is end of lag period 
    
    """
    ds_out = xarray.Dataset()
    
    lat_out = np.arange(latmin, latmax+dlat, dlat)
    lon_out = np.arange(lonmin, lonmax+dlon, dlon)
    
    pd_days = [pd.to_datetime(start_date)]
    
    nlat = len(lat_out)
    nlon = len(lon_out)
    
    mask_field = np.zeros((1, nlat, nlon))
    
    mask_field[0,3:-3,3:-3] = 1.
            
    ds_out["MASK"] = (("time", "lat","lon"), mask_field)
    
    ds_out.coords['lat'] = lat_out
    ds_out.coords['lon'] = lon_out

    ds_out.coords['time'] = pd_days
    
    ds_out.time.attrs["standard_name"] =  'time'
    ds_out.time.attrs["long_name"] =  'time'
    #ds_out.time.attrs["units"] =  'hours since 2000-01-01'
    
    ds_out["lat"].attrs["standard_name"] =  'latitude'
    ds_out["lat"].attrs["long_name"] =  'latitude'
    ds_out["lat"].attrs["units"] =  'degrees_north'
    ds_out["lat"].attrs["axis"] =  'Y'
    
    ds_out["lon"].attrs["standard_name"] =  'longitude'
    ds_out["lon"].attrs["long_name"] =  'longitude'
    ds_out["lon"].attrs["units"] =  'degrees_east'
    ds_out["lon"].attrs["axis"] =  'X'
    
    ds_out["MASK"].attrs["units"] = "unitless"
    ds_out["MASK"].attrs["long_name"] = "NAF mask"
    ds_out["MASK"].attrs["scale_factor"] = 1.
    
    ds_out.to_netcdf(path=fname_out, mode='w')
    
    print("Successfully written file: " + fname_out)

    return
    
def write_bc_ensemble_v12(start_date, start_lag, end_lag, ensemb_names, bc_dir,bc_str,
                            out_dir, species, ic_spc_names, 
                            latmin, latmax, lonmin, lonmax,
                            bc_split="None",
                            met_type = "GEOS-FP",
                            BC_ensemble = False):
    """
    Write BC input files for ensemble GEOS-Chem runs.
    
    Easiest way then is to read the CH4 field from global run. Then:
        1) Apply scaling factors of 1 for new fields
        2) apply oter scaling factos for older fields.
    
    bc scalings =  dictionary of arrays of size (nwindows, nBC) e.g. if north-south = (nwindows, 2)
    
    v12.4/v12.5
    Read in netcdf files from global run with CH4 field
    Copy this field for each ensemble member and apply scale factors from x_ensemble
    
    Only apply scale factors for length of assim window.
    But - these may require I use a longer lag window. Is there a way of avoiding this?
    
    For years not in range 2018-2019 I sholud copy one of those years. Maybe 2019. 
    
    Need a way to deal with leap years...
    
    Added a key for different met types - i.e. need different BC for MERRA2 (needs to extend 2 boxes into domain)
    
    """
    
    # Find list of file names that span date range of interest
    if start_date[:4] in ["2018", "2019"]:
        files = ensemble_mod.filenames(bc_dir, bc_str, start=start_date, end=end_lag, freq="D")
        
        start_lag2=start_lag
        end_lag2=end_lag
        
    else:
        
        start_date2 = "2018" + start_date[4:]
        
        # Need to deal with leap years
        if start_date[:4] in (["2012", "2016", "2020"]) and start_date[4:6]=="02":
            start_date2 = "20200201"
        
        if start_date[4:6] in ["11", "12"]:
            end_lag2 = "2019" + end_lag[4:]
        elif start_date[:4] in (["2012", "2016", "2020"]) and start_date[4:6]=="02":
            end_lag2 = "2020" + end_lag[4:]
        else:
            end_lag2 = "2018" + end_lag[4:]
        if start_date[4:6] == "12":
            start_lag2= "2019" + start_lag[4:]
        elif start_date[:4] in (["2012", "2016", "2020"]) and start_date[4:6]=="02":
            start_lag2= "2020" + start_lag[4:]
        else:
            start_lag2= "2018" + start_lag[4:]
        
        files = ensemble_mod.filenames(bc_dir, bc_str, start=start_date2, end=end_lag2, freq="D")
        
        date_range = pd.date_range(start=start_date, end=end_lag, freq="1D")
        
    #fname = bc_dir + "BC.20140102"
    # Now loop through all files and create new daily files in turn:

    # New version for 12.4 - Need to edit this.

    for ti, bc_file in enumerate(files):
        
        date_str = re.search("([0-9]{8})", bc_file)[0]
        ds = ensemble_mod.open_ds(bc_file)
        
        # Cut BC ds down to domain size 
        ds_out = ds.sel(lon=slice(lonmin-5., lonmax+5.), lat=slice(latmin-4, latmax+4))
        #var_field = ds_out["SpeciesBC_" + species].copy()
        
        var_field = ds_out["SpeciesBC_" + species + "BC"].copy()
        
        # Need to drop all species out of input ds
        keys = list(ds_out.keys())
        
        #print(keys)
        ds_out = ds_out.drop(keys)
             
        if int(date_str) >= int(start_lag2):   
            var_field[:,:,:,:]=0.            # Set BC to 0 from the start hour of lag period.
            
        var_field_out0 = var_field*0.
        var_field_out0.attrs = var_field.attrs
        
        ds_out["SpeciesBC_" + species] = var_field
        
        # Attributes not carried over after opearting on data array. 
        # HEMCO reads in units, so need to include attrs in copied fields.
      
        # But I do still want to include IC species?
        if len(ic_spc_names) > 0:
            
            if bc_split == "NESW":
                
                #var_field of shape (ntime,nlev,nlat,nlon)
                
                lon  = var_field.lon.values
                lat  = var_field.lat.values
                
                wh_lat0 = np.where(lat == latmin)[0][0]
                wh_lat1 = np.where(lat == latmax)[0][0]
                
                wh_lon0 = np.where(lon == lonmin)[0][0]
                wh_lon1 = np.where(lon == lonmax)[0][0]
                
                var_N = var_field_out0.copy()
                var_E = var_field_out0.copy()
                var_S = var_field_out0.copy()
                var_W = var_field_out0.copy()
                
                if met_type == "MERRA2":
                    var_N[:,:,wh_lat1-2:,wh_lon0+2:] = var_field[:,:,wh_lat1-2:,wh_lon0+2:]*1.   
                    var_E[:,:,:wh_lat1-2,wh_lon1-2:] = var_field[:,:,:wh_lat1-2,wh_lon1-2:]*1.              
                    var_S[:,:,:wh_lat0+2,:wh_lon1-2:] = var_field[:,:,:wh_lat0+2,:wh_lon1-2:]*1.
                    var_W[:,:,wh_lat0+2:,:wh_lon0+2] = var_field[:,:,wh_lat0+2:,:wh_lon0+2]
                else:
                    var_N[:,:,wh_lat1-1:,wh_lon0+1:wh_lon1-1] = var_field[:,:,wh_lat1-1:,wh_lon0+1:wh_lon1-1]*1.               
                    var_E[:,:,:,wh_lon1-1:] = var_field[:,:,:,wh_lon1-1:]*1.               
                    var_S[:,:,:wh_lat0+1,wh_lon0+1:wh_lon1-1] = var_field[:,:,:wh_lat0+1:,wh_lon0+1:wh_lon1-1]*1.
                    var_W[:,:,:,:wh_lon0+1] = var_field[:,:,:,:wh_lon0+1]
                              
                ds_out["SpeciesBC_" +species + "_BC1"] = var_N 
                ds_out["SpeciesBC_" +species + "_BC2"] = var_E 
                ds_out["SpeciesBC_" +species + "_BC3"] = var_S 
                ds_out["SpeciesBC_" +species + "_BC4"] = var_W 
                           
                for spc_str in ["BC1", "BC2", "BC3", "BC4"]:
                    ds_out["SpeciesBC_" + species + "_" + spc_str].attrs["long_name"] = "Dry mixing ratio of species " + species + "_" + spc_str
#                if bc_file==files[0]:
#                    var_N_out = var_N.copy()
#                    var_E_out = var_E.copy()
            
            for name2 in ic_spc_names:
                if name2 == species + "BC":
                    ds_out["SpeciesBC_" + name2] = var_field 
                    ds_out["SpeciesBC_" + name2].attrs["long_name"] = "Dry mixing ratio of species " + name2
                elif name2 == species + "IC":
                    ds_out["SpeciesBC_" + name2] = var_field_out0
                    ds_out["SpeciesBC_" + name2].attrs["long_name"] = "Dry mixing ratio of species " + name2
                    #ds["SpeciesBC_" + name2].attrs = ch4_field.attrs 
                #ds["SpeciesBC_" + name2].attrs["long_name"] = "Dry mixing ratio of species " + name2
                      
        # Overwrite time index if year is not in 2018 or 2019    
        if start_date[:4] in ["2018", "2019"]:
            #print("BC file exists")
            fname_out = out_dir + bc_str + date_str + "_0000z.nc4"
        else:
            # Overwrite time index
            date_out = date_range[ti].strftime("%Y%m%d")            
            ds_times = pd.to_datetime(ds_out.time.values)
            
            #ds_out.index["time"] = times_out
            
            ds_time_strs = ds_times.strftime("%Y%m%d %H")
            times_out=[]
            for time_str in ds_time_strs:
                times_out.append(pd.to_datetime(date_out[:4] + time_str[4:]) )

            ds_out.coords["time"] = times_out
            #date_str = date_out
            fname_out = out_dir + bc_str + date_out + "_0000z.nc4"
        #return ds_out, times_out      
        
        #fname_out = out_dir + bc_str + date_str + "_0000z.nc4"
        ds_out.to_netcdf(path = fname_out)
     
    print("Written BC files between " + start_date + " and " + end_lag + " to disk." )
    
    return ds_out

#%%

run_dates = ["20130101", "20130201", "20130301",
             "20130401", "20130501", "20130601",
             "20130701", "20130801", "20130901",
             "20131001", "20131101", "20131201"]
    
#run_dates = ["20140101", "20140201", "20140301",
#             "20140401", "20140501", "20140601",
#             "20140701", "20140801", "20140901",
#             "20141001", "20141101", "20141201"]

#run_dates = ["20100101"]

run_key = "outer"

species="CH4"
bc_split="NESW"
assim_window=1
window_unit = "MS"
lag_period = 1
lag_unit="MS"

write_mask = False

lonmin = 22.5                      
lonmax = 52.5                      
latmin = -8                      
latmax= 18

dlat=0.25
dlon=0.3125

met_type = "MERRA2"
#met_type = "GEOS-FP"

#if run_key == "inner":
nbasis = 120
#elif run_key  == "outer":
#    nbasis = 93

for run_date in run_dates:

    run_name = run_key + "_" + species + "_" + run_date
    
    run_root = "/home/mlunt/ceph/verify/model_settings/DARE_runs/CH4/eaf_runs/"
    run_dir = run_root + run_name + "/"
    
    #output_root = "/home/mlunt/ceph/verify/model_outputs/DARE_runs/CH4/eaf_runs/"
    #output_dir = output_root + run_name + "/"
    output_dir = run_dir + "OutputDir/satellite/"
    
    gc_code_dir = "/geos/u73/mlunt/GC_code/Code.12.5.0/"
    template_dir = run_root + "templates/"
    gc_make_dir = template_dir + "rundir_template/"
    
    restart_out_dir = run_dir + "restarts/"
    fname_restart_out = restart_out_dir + "GEOSChem.Restart." + run_date + "_0000z.nc4"
    
    bc_output_dir = run_dir + "BC/"
      
    ref_conc={"CH4": 0}
    
    restart_temp_dir = "/home/mlunt/ceph/verify/model_settings/DARE_runs/CH4/eaf_runs/restart_inputs/" 
    restart_temp_file = restart_temp_dir + "GEOSChem.Restart.20190101_0000z.nc4"
    
    #bc_input_root = "/home/mlunt/ceph/verify/model_settings/DARE_runs/BC_inputs/"
    #bc_input_dir = bc_input_root + species + "/"
    bc_input_dir = "/geos/d21/mlunt/GC_output/BC_fields_scaled/SSA/"
    
    mask_root  = "/home/mlunt/ceph/verify/model_settings/DARE_runs/CH4/eaf_runs/masks/"
    fname_mask = mask_root + "NAF_mask.nc"
    
    #%%
    # Define timings of end of emis window and end of lag run
    
    pd_run_date = pd.to_datetime(run_date)
    if window_unit == "MS":
        pd_lag_start = pd_run_date + relativedelta(months=assim_window)
    
    lag_start_date = pd_lag_start.strftime('%Y%m%d')
    if lag_unit =="d":
        pd_lag_end = pd_lag_start + pd.Timedelta(lag_period, unit="d")
    elif lag_unit =="MS":    
        pd_lag_end = pd_lag_start + relativedelta(months=lag_period)
    lag_end_date = pd_lag_end.strftime('%Y%m%d')
    
    #%%
    
    basis_names=[]
    for xi in range(nbasis):    
        basis_names.append(species + "_E" + str(xi+1))
        
    spc_IC = [species+ "IC"]
    nBC=4
    for xi in range(nBC):
        spc_IC.append(species + "_BC" + str(xi+1))
    
    # Add reference species concentration    
    #spc_IC.append(species + "REF")
#    bc_str = "GEOSChem.BoundaryConditions."
#    write_bc_ensemble_v12(run_date, lag_start_date, lag_end_date, basis_names, bc_input_dir, bc_str,
#                            bc_output_dir, species, spc_IC,
#                            latmin,latmax,lonmin,lonmax,
#                            bc_split = bc_split, BC_ensemble=False)
    
#    a = ds_time.strftime("%Y%m%d %H")
#    b=[]
#    for string in a:
#        b.append(pd.to_datetime(dates_out[:4] + string[4:]) )
    #%%
    
    # Create directory structure
    subprocess.call(["mkdir", '-p', run_dir])
    subprocess.call(["mkdir", '-p', output_dir])
    
    subprocess.call(["mkdir", '-p', run_dir + 'input_files'])
    subprocess.call(["mkdir", '-p', run_dir + 'OutputDir'])
    
    subprocess.call(["ln", '-s', gc_code_dir, run_dir + "CodeDir"])
    
    make_files = glob.glob(gc_make_dir + "/*")
    
    for mfile in make_files:
        subprocess.call(["cp", mfile, run_dir]) # Can't use * wildcard without shell interpretation. 
    
    #subprocess.call(["mkdir", '-p', bc_run_dir])
    subprocess.call(["mkdir", '-p', bc_output_dir])
    subprocess.call(["mkdir", '-p', restart_out_dir])
    
    #%%
    # Write restart file
    ds_restart = ensemble_mod.write_restart_file_v12(restart_temp_file, 
                                                                   fname_restart_out, species, basis_names,spc_IC, 
                                                                   run_date, ref_conc = ref_conc[species], 
                                                                   write_IC = True,
                                                                   write_ens=True)
       
    #%%
    # Write BC files
    bc_str = "GEOSChem.BoundaryConditions."
    check_BC=False
    write_bc_files = True
    if check_BC:
        # Check whether BC files from global run exist or not. Raise an error if not.
        bc_check = ensemble_mod.filenames(bc_input_dir, bc_str, start=run_date, end=lag_end_date, freq="D")
        if len(bc_check) == 0:
            raise ValueError("No BC files available. Make sure you run global model first before attempting nested ensemble run. Or turn check_BC to False") 
       
    if write_bc_files == True:
        
        # Only write BCs for CH4 species
        # All ensemble fields have BC=0. Only interested in emissions component.
        ds_bc_out = write_bc_ensemble_v12(run_date, lag_start_date, lag_end_date, basis_names, bc_input_dir, bc_str,
                            bc_output_dir, species, spc_IC,
                            latmin,latmax,lonmin,lonmax,
                            bc_split = bc_split, 
                            met_type=met_type, BC_ensemble=False)
    
    #%%
    
    # Write input and HEMCO window and lag files and history and HEMCO diagn (just need to copy these ones)
    write_input_file(run_date, lag_start_date, 
                                          run_dir, output_dir, template_dir,  
                                          species, run_key, met_type, lag=False)
    
    # Write the input file for the lag period
    write_input_file(lag_start_date,lag_end_date, 
                                  run_dir, output_dir, template_dir,  
                                  species, run_key, met_type, lag=True)
    
    #%%     
    # Write the HEMCO_Config file for GEOS-Chem -  just need to copy the HEMCO file. Emis files already writen
           
    # 1. Pick correct HEMCO template file
    # 2. Copy to run_dir/input_files
    if met_type == "MERRA2":
        hemco_template = template_dir + run_key + "_" + species + "_merra2_HEMCO_Config.template"
    else:
        hemco_template = template_dir + run_key + "_" + species + "_HEMCO_Config.template"
    hemco_window_file = run_dir + "input_files/" + "window_HEMCO_Config.rc"
    
    subprocess.call(["cp", hemco_template, hemco_window_file])
           
    # Write hemco_lag
    # 1. Pick correct species HEMCO template file
    # 2. Copy to run_dir/input_files
    if met_type == "MERRA2":
        hemco_lag_template = template_dir + "lag" + "_" + species + "_merra2_HEMCO_Config.template"
    else:
        hemco_lag_template = template_dir + "lag" + "_" + species + "_HEMCO_Config.template"
    hemco_lag_file = run_dir + "input_files/" + "lag_HEMCO_Config.rc"
    subprocess.call(["cp", hemco_lag_template, hemco_lag_file])
    
    # Write the HISTORY.rc file for GEOS-Chem
    # Again just copy
    
    history_template = template_dir + "HISTORY.template"
    history_file = run_dir + "HISTORY.rc"
    subprocess.call(["cp", history_template, history_file])
    
    # Write HEMCO diagnostics file
    # Again just copy the template file
    hemco_diagn_template = template_dir + run_key + "_" + species + "_HEMCO_Diagn.template"
    diagn_file = run_dir + "HEMCO_Diagn.rc"
    subprocess.call(["cp", hemco_diagn_template, diagn_file])
    
    #if write_mask == True:
    #    
    #    write_mask_file(latmin,latmax,lonmin,lonmax,dlat,dlon, fname_mask, run_date)
    
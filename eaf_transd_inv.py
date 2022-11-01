#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 11:08:38 2020

Code for CO2-only inversion

Read in site files and combine 

Try firstly solving for bulk CO2 (CO2 + CO2 ff together)

Then try solving separately, but make sure it works first.- i.e with pseudo data

@author: mlunt
"""
import numpy as np
import xarray
import matplotlib.pyplot as plt
import glob
#import re
import pandas as pd
import json
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colorbar as cbar
from areagrid import areagrid
from dateutil.relativedelta import relativedelta
import eaf_inv_mod
import eaf_hbmcmc_mod
import subprocess
import argparse

def open_ds(fname):
    with xarray.open_dataset(fname) as ds:
        ds.load()
    return ds

#%%
parser = argparse.ArgumentParser(description='This is a demo script by Mark.')              
parser.add_argument("start", help="Start date string yyyymmdd") 
parser.add_argument("end", help="End date string yyyymmdd") 
parser.add_argument("version", help="Run version name")
args = parser.parse_args()

start_date = args.start
end_date = args.end
version = args.version


#start_date = "20200201"
#end_date = "20200229"
#version = "eaf_transd_brem_lst2_050cm"

#satellite = "GOSAT"
satellite = "TROPOMI"

bremen= True

if start_date == "20140101":
    start0 = "20140101"
else:
    start0 = "20100101"

data_root = "/home/mlunt/ceph/verify/model_settings/DARE_runs/CH4/eaf_runs/"
run_groups = ["outer", "inner"]


nbasis={}
nbasis["outer"] = 98
nbasis["inner"] = 120
species= "CH4"
nBC=4

soil_depth = "050cm"

emis_key = "ch4_total_lst2"  #_lst2"
#emis_key = "ch4_total"  # Wetcharts version
scale_prior =False
prior_scaling = 1.  #*1.5   # Make flat prior 50% larger

if satellite == "TROPOMI" and bremen == True:
    sat_fstr = "XCH4_model_" + satellite + "_bremen_"
else: 
    sat_fstr = "XCH4_model_" + satellite + "_"

#Thresholds for TROPOMI filtering
albedo_thresh0=0.05    #0.06
albedo_thresh1=0.35    #0.2
aot_thresh=0.1
asize_thresh=3.4
rough_thresh=200
#%%
##########################
# Set up transd parameters here:
transd=True
kmin=10             # Minimum number of regions
kmax=200           # Maximum number of regions
#kmax=200
k_ap = 50         # Starting number of regions

# Setup MCMC parameters here
#################
###############################
# Define MCMC parameters
para_temp=True   # 0 = off 1 = on

nIt = 20000   # 100000
burn_in = 20000   #100000
nsub = 100    # 100

# Paramters for mcmc inversion
#################################
# Define nbeta

nbeta = 8   # Number of parallel chains - needs to be defined even if no Parallel tempering
series = np.linspace(0.,1.,num=nbeta)
beta_1 = np.exp(0 + (np.log(250)-0)*series)
beta= 1./beta_1

####################################################
# Define stepsizes
stepsize_sig_y_0 = 0.2  # units of ppb
stepsize_0 = 0.5
stepsize_pdf_p1_0 = 0.1
stepsize_pdf_p2_0 = 0.1
#####################################################
# Set up hyper-parameters

x_pdf = 3   # 3  = lognormal
pdf_param1_pdf = 2
pdf_param2_pdf=2

pdf_param1_0 = 1.
pdf_param2_0 = 0.5

pdf_p1_hparam1_0 = pdf_param1_0*1.
pdf_p1_hparam2_0 = 0.1
pdf_p2_hparam1_0 = pdf_param2_0*1.
pdf_p2_hparam2_0 = 0.2

#########################################################
# Sey up sigma_y hyper-parameters

#sigma_model_ap = 25.    # Was 5.
if satellite == "GOSAT":
    sigma_model_ap = 10.    # Was 5.
elif satellite == "TROPOMI":
    sigma_model_ap= 20.

sigma_y_period = 31   # Length of each estimation period fro sigma_y in days. 
                     # The longer it is the quicker inversion will run
sigma_model_pdf = 2   # normal
#sigma_model_pdf = 3   # lognormal

sigma_model_hparam1_0 = sigma_model_ap*1.   # Units == ppb
sigma_model_hparam2_0 = sigma_model_ap/5.   # units = ppb 

#sigma_model_hparam1_0 = 0.5   # Units == ppb
#sigma_model_hparam2_0 = 40.   # units = ppb   # Was 30

stepsize_clon = 6.    # Stepsize for longitude for move
stepsize_clat = 4.    # Stepsize for latitude for move
stepsize_bd=1.5        # Stepsize for change in x during birth step

####################################################

#%%

#sites= ["TAC", "RGL", "BSD", "HFD"]
#sites=["NOR"]
inv_out_dir = "/home/mlunt/datastore/EAF/inv_outputs/paper/"  + version + "/"

BC_names = []
for xi in range(nBC):
    BC_names.append(species + "_BC" + str(xi+1))    

spc_names={}
for run_str in run_groups:
    spc_names[run_str]=[]
    for xi in range(nbasis[run_str]):
        #spc_names.append("CO2R" + str(xi+1).zfill(2) )
        spc_names[run_str].append(species + "_E" + str(xi+1))
 

spc_mod_dict1 ={}
spc_mod_dict2 ={}
for spc in spc_names:
    spc_mod_dict1[spc]=[]
    spc_mod_dict2[spc]=[]
 
# Set up run dates
#This assumes each assimilation window is 1 month
# Set up which dates to read in. Max 1 assim window and 1 lag windows

if start_date == start0:
    run_dates = [start_date]
else:    
    pd_run_date0  = pd.to_datetime(start_date) - relativedelta(months=1)
    run_date0 = pd_run_date0.strftime("%Y%m%d")
    run_dates   = [run_date0, start_date]
        
nruns = len(run_dates)

if start0 in run_dates:
    nIC = 1
else:
    nIC=0
    
ngroups = len(run_groups)

#%%

H_basis_rdate={}
H_bc_rdate = {}
x_group_rdate={}
for run_date in run_dates:
   
    H_basis_rdate[run_date] = []
    
    if run_date == start0:
        x_group_rdate[run_date]=["IC"]
    else:
        x_group_rdate[run_date]=[]
    #H_bc_run_date[run_date] = []
    
    
    for xi in range(nBC):
        x_group_rdate[run_date].append("BC")
        
    for run_str in run_groups:
        run_name=run_str + "_" + species + "_" + run_date   
        site_out_dir = data_root + run_name + "/OutputDir/satellite/" +satellite + "/"  
        
        if run_str == "outer":
            if run_date == start_date:
                y_bc_stack, y_data_dict = eaf_inv_mod.read_sat_output(species, site_out_dir, sat_fstr, 
                                                                      run_date, start_date, end_date, BC_names, 
                                                                      read_obs=True, satellite = satellite, bremen=bremen)
                
            else:
                y_bc_stack = eaf_inv_mod.read_sat_output(species, site_out_dir, sat_fstr,
                                                         run_date, start_date, end_date, BC_names, 
                                                         read_obs=False, satellite = satellite, bremen=bremen)
                
            H_bc_rdate[run_date] = y_bc_stack
            
        # Read in fixed Basis_function names
        y_basis_stack = eaf_inv_mod.read_sat_output(species,  site_out_dir, sat_fstr,
                                                    run_date, start_date, end_date, spc_names[run_str], 
                                                    read_obs=False, satellite = satellite, bremen=bremen)
        
        # Read in IC term if needed.
        if run_date == start0:
            y_ic_stack = eaf_inv_mod.read_sat_output(species, site_out_dir, sat_fstr,
                                                     run_date, start_date, end_date, [species + "IC"], 
                                                     read_obs=False, satellite=satellite, bremen=bremen)
            
            y_IC = y_ic_stack.squeeze()
        
    
        H_basis_rdate[run_date].append(y_basis_stack)
        
        
    
        for xi in range(nbasis[run_str]):
            x_group_rdate[run_date].append(run_str)
        
    
#%%    

# Aggregate model outputs into one Jacobian matrix 
nobs_unfilt = len(y_data_dict["obs"])
nbasis_all = nbasis["inner"]+nbasis["outer"]
nstate = (nBC + nbasis_all)*nruns + nIC
nstate0 = nBC+nbasis_all

H_unfilt = np.zeros((nstate,nobs_unfilt))


if nIC == 1:
    H_unfilt[0,:] = y_IC
    
for ri, run_date in enumerate(run_dates):
    H_unfilt[nIC + ri*(nstate0): nIC + ri*(nstate0)+nBC] = H_bc_rdate[run_date]
    
    cum_nbasis=0
    for xi, run_str in enumerate(run_groups):
           
            H_unfilt[nIC + ri*(nstate0)+nBC + cum_nbasis: nIC+ ri*(nstate0)+nBC+cum_nbasis + nbasis[run_str]] = H_basis_rdate[run_date][xi]
            cum_nbasis+= nbasis[run_str] 
            
x_group_list=[]
for run_date in run_dates:
    x_group_list.append(x_group_rdate[run_date])
   
x_group = np.hstack(x_group_list)


#%%

if satellite == "TROPOMI":
   
    albedo = y_data_dict["albedo"]
    y_lat_temp = y_data_dict["lat"]
    y_lon_temp = y_data_dict["lon"]
    y_date_temp = y_data_dict["date"]
    
    if bremen == True:
        surf_rough = y_data_dict["surf_rough"]
    
#        wh_data = np.where(( albedo >= albedo_thresh0) & 
#                       (albedo < albedo_thresh1)  &
#                       (surf_rough < rough_thresh) & 
#                       (y_lat_temp <= 14.5) & 
#                       (y_lon_temp < 42.) )[0]
#        
        wh_data = np.where(( albedo >= albedo_thresh0) & 
                       (albedo < albedo_thresh1)  &
                       (surf_rough < rough_thresh) )[0]
    else:
        asize = y_data_dict["asize"]
        aot = y_data_dict["AOT"]
        xco = y_data_dict["xco"]
        
#        wh_data = np.where(( albedo >= albedo_thresh0) & 
#                           (albedo < albedo_thresh1)  &
#                           (aot < aot_thresh) & 
#                           (asize > asize_thresh) &
#                           (y_lat_temp <= 14.5) & 
#                           (y_lon_temp < 42.) )[0]
#        
        # Less striict filtering. Added 21/7. Need to test impact and compare to GOSAT
        wh_data = np.where(( albedo >= albedo_thresh0) & 
                           (albedo < albedo_thresh1)  &
                           (aot < aot_thresh) 
                            )[0]
    
    # MAybe more selective over lats and lons
    
    # What if I time average the TROPOMI data? 
    # Might reduce 1. noise, 2. data dimension 3. help long-lived signals stand out
    
    y_obs_temp = y_data_dict["obs"][wh_data]
    y_std_temp = y_data_dict["std"][wh_data]*4.
    y_lat_temp = y_data_dict["lat"][wh_data]
    y_lon_temp = y_data_dict["lon"][wh_data]
    y_albedo_temp = y_data_dict["albedo"][wh_data]
    
    # Loop through unique y_lat and unique_y_lon
    y_lat_unq = np.unique(y_lat_temp)
    y_lon_unq = np.unique(y_lon_temp)
 
    H_unfilt_temp = H_unfilt[:,wh_data]
    
    y_obs_av_list=[]
    y_std_av_list=[]
    y_lat_av_list=[]
    y_lon_av_list=[]
    H_av_list=[]
    y_IC_av_list=[]
    y_albedo_list=[]
    y_date_list=[]
    for lati in y_lat_unq:
        for loni in y_lon_unq:
            
            wh_unq = np.where((y_lat_temp == lati) & (y_lon_temp == loni))[0]

            if len(wh_unq)> 0:
                y_obs_av_list.append(y_obs_temp[wh_unq].mean())
                y_std_av_list.append(y_std_temp[wh_unq].mean())
                y_lat_av_list.append(lati)
                y_lon_av_list.append(loni)
                y_albedo_list.append(y_albedo_temp[wh_unq].mean())
                y_date_list.append(y_date_temp[wh_unq].min())
                
                H_av_list.append(H_unfilt_temp[:,wh_unq].mean(axis=1))
                
                if start_date == start0:
                    y_IC_av_list.append(y_ic_stack.squeeze()[wh_unq].mean())
                    y_IC = np.hstack(y_IC_av_list)
                    
    H_temp = np.vstack(H_av_list)
    y_obs = np.hstack(y_obs_av_list)
    y_std = np.hstack(y_std_av_list)
    y_lat = np.hstack(y_lat_av_list)
    y_lon = np.hstack(y_lon_av_list)
    y_albedo = np.hstack(y_albedo_list)
    y_date = np.hstack(y_date_list)
    
    # Need to better filter out data over Yemen.
    
    print(species, len(y_data_dict["obs"]))
    print("Well mixed:")
    print(len(wh_data))
    print ("Every nth:")
    print(len(y_obs))

    #print ("Every 4th:")
    #print(len(wh_data[::4]))
    
    # Set ~3000 data points as a max
#    if len(wh_data) > 4200:
#        nskip = int(np.round(len(wh_data)/4000.))   # Wsa 3200 and 3000
#    else:
#        nskip=1
#    
#    base=1
#    print ("Every nth:")
#    print(len(wh_data[base::nskip]))
#    #if len(wh_data)/3000 > 20000:
#    
#    y_obs = y_data_dict["obs"][wh_data[base::nskip]]
#    y_std = y_data_dict["std"][wh_data[base::nskip]]*2.
#    y_date = y_data_dict["date"][wh_data[base::nskip]]
#    y_lat = y_data_dict["lat"][wh_data[base::nskip]]
#    y_lon = y_data_dict["lon"][wh_data[base::nskip]]
#    y_xco = y_data_dict["xco"][wh_data[base::nskip]]
#    y_albedo = y_data_dict["albedo"][wh_data[base::nskip]]
#    y_aot = y_data_dict["AOT"][wh_data[base::nskip]]
#    if start_date == start0:
#        y_IC = y_ic_stack.squeeze()[wh_data[base::nskip]]
#      
#    H_temp = np.transpose(H_unfilt[:,wh_data[base::nskip]])
    
else:
    flag = y_data_dict["flag"]
    wh_flag = np.where(flag <2)[0]    #0 = land; 1 = glint
       
    y_obs = y_data_dict["obs"][wh_flag]
    y_std = y_data_dict["std"][wh_flag]*1.5
    y_date = y_data_dict["date"][wh_flag]
    y_lat = y_data_dict["lat"][wh_flag]
    y_lon = y_data_dict["lon"][wh_flag]
    
    if start_date == start0:
        y_IC = y_ic_stack.squeeze()[wh_flag]
      
    H_temp = np.transpose(H_unfilt[:,wh_flag])
    
    #H_temp= np.transpose(H_unfilt)
    
nobs = len(y_obs)  


#dates = np.unique(y_date)
# End of obs read section    

#%%
# Read in emissions file and apply emissions field to H
#emis_file = "/home/mlunt/datastore/EAF/emissions/ch4_emis_eaf_025x03125_2010_2021.nc"
#emis_file = "/home/mlunt/datastore/EAF/emissions/ch4_emis_eaf_v2_Cs_025x03125_2010_2021.nc"
#emis_file = "/home/mlunt/datastore/EAF/emissions/ch4_emis_eaf_v2b_Cs_025x03125_2010_2021.nc"

emis_file = "/home/mlunt/datastore/EAF/emissions/ch4_emis_eaf_v2c_Cs_" + soil_depth + "_025x03125_2010_2021.nc"

ds_emis = open_ds(emis_file)

emis_map_ap={}
for run_date in run_dates:
    emis_map_ap[run_date] = ds_emis[emis_key].sel(time=run_date)

# Read in basis functions
# Map posterior onto basis functions
basis_dir = "/home/mlunt/ceph/verify/model_settings/DARE_runs/CH4/eaf_runs/masks/"

fname_inner = basis_dir + "eaf_bfs_inner.nc"
fname_outer = basis_dir + "eaf_bfs_outer.nc"

ds_inner = open_ds(fname_inner)
ds_outer = open_ds(fname_outer)

outer_keys = list(ds_outer.keys())
inner_keys = list(ds_inner.keys())

# Loop through x_post and bf inner
wh_outer = np.where(x_group == "outer")[0]
wh_inner = np.where(x_group == "inner")[0]

wh_outer_rdate={}
wh_inner_rdate={}
for xi, run_date in enumerate(run_dates):
    
    wh_outer_rdate[run_date] = wh_outer[nbasis["outer"]*xi:nbasis["outer"]*(xi+1) ]
    wh_inner_rdate[run_date] = wh_inner[nbasis["inner"]*xi:nbasis["inner"]*(xi+1) ]

x_post_map = ds_outer.bf_01[0,:,:]*0.
x_post_map2 = ds_outer.bf_01[0,:,:]*0.

q_basis = np.zeros((nstate))+1

for xi,basis in enumerate(outer_keys):
    for run_date in run_dates:
        dum1 = ds_outer[basis][0,:,:]*emis_map_ap[run_date]
        q_basis[wh_outer_rdate[run_date][xi]] = (dum1.where(ds_outer[basis][0,:,:] ==1) ).mean()/1.e-9 
        # Divide by 1.e-9 to account for flat emissions rate used in prior model runs
    
for xi,basis in enumerate(inner_keys):
    for run_date in run_dates:
        dum1 = ds_inner[basis][0,:,:]*emis_map_ap[run_date]
        q_basis[wh_inner_rdate[run_date][xi]] = (dum1.where(ds_inner[basis][0,:,:] ==1) ).mean()/1.e-9
        # Divide by 1.e-9 to account for flat emissions rate used in prior model runs


# Take mean emissions in each basis function region and multiply by H.
# Then prior always starts from 1. Or have x_ap as the mean emissions (in units of 1e-9?) in each grid cell
H_temp2 = H_temp*q_basis

#H_temp2=H_temp.copy()


#%%
# Aggregate run_dates so I only have one set of estimates
H = np.zeros((nobs, nstate0 +nIC))

if nIC >0:
    H[:,0] = H_temp2[:,0]

for ri in range(nruns):
    H[:,nIC:] = H[:,nIC:] + H_temp2[:,ri*nstate0 +nIC:nstate0*(ri+1)+nIC]

Ngrid = nbasis_all
x_group2 = x_group[:nstate0+nIC]

#%%
if satellite == "TROPOMI":
    sigma_model = 12.
else:
    sigma_model= 4.
sigma_y = np.sqrt(y_std**2 + sigma_model**2)

wh_BC = np.where((x_group2 == "IC") | (x_group2 == "BC"))

x_ap = np.zeros((nstate0+nIC))+1   #+0.1
x_ap[wh_BC] = 1.
x_uncert = np.zeros((nstate0+nIC))+0.4


#%%
# Need to set up for transd inversion
# With multiple run dates need to work out how I do that.
# Obvious option is to apply same H_agg to both dates. 
# Other option is to add rdates together. So only have 1 set of terms. Think this migth be easiest approach.
# No need to recode fortran then. 
H_bc = H[:,:nIC+nBC]
h_v = H[:,nBC+nIC:]

x_agg = np.ones((k_ap+nIC+nBC))
x_v = np.ones((Ngrid + nIC + nBC))
y_mod_ap = np.dot(H,x_v)
n0_ap = y_mod_ap - y_obs  # Need to make sure the right way round.


#%%
# Define lats and lons of "raw grid"
# Need centre points of each basis function
# Not regular so need to loop through outer then inner. Find each basis mid point.
dlon = 0.3125/2.
dlat= 0.25/2.
lats=np.zeros((Ngrid))
lons = np.zeros((Ngrid))
for xi,basis in enumerate(outer_keys):
    basis = ds_outer[basis][0,:,:]
    basis2 = basis.where(basis == 1, drop=True)
    basis_lat = (basis2.lat.values-dlat).mean()
    basis_lon = (basis2.lon.values-dlon).mean()
    lats[xi] = basis_lat
    lons[xi] = basis_lon
    
for xi,basis in enumerate(inner_keys):
    basis = ds_inner[basis][0,:,:]
    basis2 = basis.where(basis == 1, drop=True)
    basis_lat = (basis2.lat.values-dlat).mean()
    basis_lon = (basis2.lon.values-dlon).mean()
    lats[xi+nbasis["outer"]] = basis_lat
    lons[xi+nbasis["outer"]] = basis_lon



#%%

###################################
# Call MCMC routine here
#pdf_param1_ap = x_ap * pdf_param1_0
#pdf_param2_ap = x_ap*0. + pdf_param2_0
#pdf_p1_hparam1 = x_ap * pdf_p1_hparam1_0
#pdf_p1_hparam2 = x_ap*0 + pdf_p1_hparam2_0
#pdf_p2_hparam1 = pdf_param2_ap * 1.
#pdf_p2_hparam2 = pdf_param2_ap * pdf_p2_hparam2_0

#stepsize = (x_ap *0. +1)* stepsize_0
#stepsize_pdf_p1 = pdf_param1_ap * stepsize_pdf_p1_0
#stepsize_pdf_p2 = pdf_param2_ap * stepsize_pdf_p2_0

# Need to pass in lon and lat ibins for sigma_y  -do later


mcmc_out_dict = eaf_hbmcmc_mod.run_transd_mcmc(nIC,nBC, k_ap, kmin, kmax, beta, x_agg,
                    y_obs, n0_ap, y_std, y_lat, y_lon,
                    H_bc, h_v, lats, lons,
                    nIt, nsub, burn_in,
                    pdf_param1_0, pdf_param2_0,
                    pdf_p1_hparam1_0, pdf_p1_hparam2_0,
                    pdf_p2_hparam1_0, pdf_p2_hparam2_0,
                    x_pdf, pdf_param1_pdf, pdf_param2_pdf,    
             stepsize_0, stepsize_pdf_p1_0, stepsize_pdf_p2_0,
             sigma_model_pdf, sigma_model_ap, sigma_model_hparam1_0, sigma_model_hparam2_0,
             stepsize_sig_y_0, 
             stepsize_clat, stepsize_clon, stepsize_bd,
             para_temp=para_temp, rjmcmc=transd
             )

#mcmc_out_dict = eaf_hbmcmc_mod.run_mcmc(start_date, end_date, nIC, nBC, ngroups, nstate0, nruns,
#             sigma_y_period, dates, nIt, nsub, burn_in,
#             beta, x_ap, H, y_obs, y_std,y_date,y_mod_ap, y_lat,y_lon,
#             x_pdf_all, pdf_param1_pdf, pdf_param2_pdf, pdf_param1_ap, pdf_param2_ap, 
#             pdf_p1_hparam1, pdf_p1_hparam2,
#             pdf_p2_hparam1, pdf_p2_hparam2,
#             stepsize, stepsize_pdf_p1, stepsize_pdf_p2,
#             sigma_model_pdf, sigma_model_ap, sigma_model_hparam1_0, sigma_model_hparam2_0,
#             stepsize_sig_y_0, rho,
#             para_temp=False
#             )

#%%

x_it= mcmc_out_dict["x_it"]
y_post_it = mcmc_out_dict["y_it"]
sigma_y_it = mcmc_out_dict["sigma_y_it"]
sigma_model_it = mcmc_out_dict["sigma_model_it"]
x_post_mcmc = np.mean(x_it,axis=1)
y_post_mcmc = np.mean(y_post_it,axis=1)
#pdf_param2_mean = np.mean(mcmc_out_dict["pdf_param2_it"], axis=1)
ratio_out = mcmc_out_dict["ratio"]

ratio_sig_y_out = mcmc_out_dict["ratio_sigma_y"]
ratio_birth = mcmc_out_dict["ratio_birth"]
ratio_death = mcmc_out_dict["ratio_death"]
ratio_move = mcmc_out_dict["ratio_move"]


#fig,ax = plt.subplots()

#ax.plot(y_mod_ap, label="Model prior")
#ax.plot(y_obs, 'o', markersize=2, label="ICOS and DECC obs")
#ax.plot(y_post, label = "Model posterior")

#ax.plot(y_post_mcmc,  label = "MCMC posterior")
#ax.plot(y_obs, 'o', markersize=2, label="ICOS and DECC obs")

#ax.set_ylabel("CO2 (ppm)")
#ax.set_xlabel("Observation count")

#leg = ax.legend()

#%%
## Map posterior onto basis functions

# First need to map x onto regions for each iteration

# i.e firstly turn x_it back to x_v_it
k_it = mcmc_out_dict["k_it"]
regions_it = mcmc_out_dict["regions_it"]

nIt_out  = len(k_it)
x_post_v_it = np.zeros((Ngrid,nIt_out))
for it in range(nIt_out):
    
    k_i = k_it[it]
    for reg in range(k_i):
        wh = np.where(regions_it[:,it] == reg+1)[0]  # Need to add 1 to account for Fortran indexing
        x_post_v_it[wh,it] = x_it[reg+nIC+nBC,it]
        

## Then average x_v_it(dim="nIt")
x_post_v_mean = x_post_v_it.mean(axis=1)
# Then use process below to map

if len(wh_outer) > nbasis["outer"]:
    wh_outer2 = wh_outer[nbasis["outer"]:]
    wh_inner2 = wh_inner[nbasis["inner"]:]
else:
    wh_outer2=wh_outer.copy()
    wh_inner2=wh_inner.copy()

x_post_map = ds_outer.bf_01[0,:,:]*0.
x_post_map2 = ds_outer.bf_01[0,:,:]*0.
q_post_map = ds_outer.bf_01[0,:,:]*0.
q_post_map2 = ds_outer.bf_01[0,:,:]*0.
q_ap_map = ds_outer.bf_01[0,:,:]*0.
H_map = ds_outer.bf_01[0,:,:]*0.
for xi,basis in enumerate(outer_keys):
    x_post_map = x_post_map + ds_outer[basis][0,:,:]*x_post_v_mean[wh_outer[xi]-nBC-nIC]  #*q_basis[wh_outer2[xi]]
    q_post_map = q_post_map + ds_outer[basis][0,:,:]*x_post_v_mean[wh_outer[xi]-nBC-nIC]*q_basis[wh_outer2[xi]]*1.e-9
    q_ap_map = q_ap_map + ds_outer[basis][0,:,:]*q_basis[wh_outer2[xi]]*1.e-9
    
    #H_map = H_map + ds_outer[basis][0,:,:]*H[800,wh_outer[xi]] 
for xi,basis in enumerate(inner_keys):
    x_post_map = x_post_map + ds_inner[basis][0,:,:]*x_post_v_mean[wh_inner[xi]-nBC-nIC]  #*q_basis[wh_inner2[xi]]
    q_post_map = q_post_map + ds_inner[basis][0,:,:]*x_post_v_mean[wh_inner[xi]-nBC-nIC]*q_basis[wh_inner2[xi]]*1.e-9
    q_ap_map = q_ap_map + ds_inner[basis][0,:,:]*q_basis[wh_inner2[xi]]*1.e-9
    
    #H_map = H_map + ds_inner[basis][0,:,:]*H[800,wh_inner[xi]] 

#%%
lon=ds_outer.lon.values
lat=ds_outer.lat.values
dlon=lon[1]-lon[0]
dlat=lat[1]-lat[0]
cmin=0
cmax=2.
cmin2=-1.e-9
cmax2=1.e-9     
#proj = ccrs.PlateCarree()
#fig2,ax2 = plt.subplots(2, figsize=(8,6), subplot_kw=dict(projection=proj))
#
#p2=ax2[0].pcolormesh(lon-dlon/2, lat-dlat/2, x_post_map, 
#            transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=cmin, vmax=cmax)
#
#ax2[0].coastlines()
#ax2[0].add_feature(cfeature.BORDERS, edgecolor='black', alpha=0.8)
#
#p2b=ax2[1].pcolormesh(lon-dlon/2, lat-dlat/2, q_post_map-q_ap_map, 
#            transform=ccrs.PlateCarree(), cmap='RdBu_r', vmin=cmin2, vmax=cmax2)
#ax2[1].coastlines()
#ax2[1].add_feature(cfeature.BORDERS, edgecolor='black', alpha=0.8)
#
#
#cbaxes2 = fig2.add_axes([0.1, 0.1, 0.8, 0.04]) 
##[left, bottom, width, height],
#cb = plt.colorbar(p2, cax = cbaxes2, orientation='horizontal', extend='both', label = 'Emissions')  

#%%
# calculate emissions sum.

area = areagrid(lat,lon)
q_sum = (q_post_map*area).sum()*60*60*24*365./1.e9    #x_post is in units of 1.e-9
#q_sum2 = (q_post_map2*area).sum()*60*60*24*365./1.e9

q_sum_ap = (q_ap_map*area).sum()*60*60*24*365./1.e9 
#%%

wh_ch4bc = np.where((x_group2 == "BC") | (x_group2 =="IC"))[0]
y_ap_ch4bc = np.dot(H[:,:nBC+nIC],x_ap[:nBC+nIC])
y_post_ch4bc = np.dot(H[:,:nBC+nIC],x_post_mcmc[:nBC+nIC])

# Save outputs
# x_post
# x_uncert
# y_post

# Split into different groups

ds_out = xarray.Dataset({"y_post_it":(["nobs", "nIt"], y_post_it),
                         "y_mod_ap":(["nobs"], y_mod_ap),
                         "y_lat": (["nobs"], y_lat), 
                         "y_lon": (["nobs"], y_lon), 
                         "y_obs":(["nobs"], y_obs),
                         "y_date": (["nobs"], y_date),   
                         "y_obs_std": (["nobs"], y_std),
                         "sigma_y_it":(["nobs", "nIt"], sigma_y_it),
                         "y_post_bc": (["nobs"], y_post_ch4bc),
                         "y_ap_bc": (["nobs"], y_ap_ch4bc),
                         "k_it": (["nIt"], k_it),
                         "x_bc_it":(["nBC", "nIt"], x_it[:nBC+nIC,:]),
                         #"x_it":(["basis", "nIt"], x_post_v_it),
                         #"x_ap_uncert":(["regions"], pdf_param2_mean),
                         "x_group": (["regions"], x_group2),
                         
                         "groups":(["ngroups"], run_groups),
                         "run_dates":(["nruns"], run_dates),
                         
                         "acc_ratio": (["nBC1"], ratio_out),
                         "sigma_y_ratio":  ratio_sig_y_out,
                         "ratio_birth": ratio_birth,
                         "ratio_death": ratio_death,
                         "ratio_move": ratio_move,
                         
                         "nIC":nIC,
                         #"nBC": nBC,
                         "burn_in": burn_in,
                         "nsub": nsub
                         },
                         )
#ds_out.attrs["comment"] = "Output from 31/08/21 with lognormal sigma_model and prior 20ppb to eb consistent with NAME rjmcmc"

cum_nbasis=0
for ri, run_str in enumerate(run_groups):
    
    #x_group = x_post[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*ri: 
    #    (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*(ri+1)]
    #x_uncert_group = x_post_uncert[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*ri: 
    #    (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*(ri+1)]
    
#    H_unfilt[nIC + ri*(nstate0)+nBC + cum_nbasis: nIC+ ri*(nstate0)+nBC+cum_nbasis + nbasis[run_str]] = H_basis_rdate[run_date][xi]
#            cum_nbasis+= nbasis[run_str] 

    x_ap_group = x_ap[nIC +nBC + cum_nbasis:nIC +nBC + cum_nbasis + nbasis[run_str]]

    x_it_group = x_post_v_it[cum_nbasis:cum_nbasis + nbasis[run_str],:]
        
    #pdf_p2_group = pdf_param2_mean[nIC +nBC + cum_nbasis: 
    #     nIC +nBC + cum_nbasis + nbasis[run_str]]
        
    cum_nbasis+= nbasis[run_str] 
#    x_ap_group = x_ap[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*ri: 
#        (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*(ri+1)]
#    x_it_group = x_it[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*ri: 
#        (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*(ri+1),:]
#        
#    pdf_p2_group = pdf_param2_mean[(nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*ri: 
#        (nBC + nbasis*ngroups)*(nruns-1) + nIC +nBC + nbasis*(ri+1)]
        
    ds_out["x_it_"  + run_str] = (("basis_" + run_str , "nIt"), x_it_group)
    ds_out["x_ap_"  + run_str] = (("basis_" + run_str), x_ap_group)
    #ds_out["x_ap_uncert_"  + run_str] = (("basis_"+ run_str), pdf_p2_group)
    
    ds_out["x_it_"  + run_str].encoding['least_significant_digit'] = 3 
    ds_out["x_ap_"  + run_str].encoding['least_significant_digit'] = 2 
    #ds_out["x_ap_uncert_"  + run_str].encoding['least_significant_digit'] = 2 
    #ds_out["x_post_" + run_str] = (("basis"), x_group)
    #ds_out["x_uncert_" + run_str] = (("basis"), x_uncert_group)

    # Also save out y_mod_ap_group to file
    # And y_post_group  - want to see how much each term is contributing to measurements at different sites


fname_out = inv_out_dir + "inv_out_"  + version + "_" +  start_date + ".nc"

if len(glob.glob(inv_out_dir)) ==0:
    subprocess.call(["mkdir", '-p', inv_out_dir])

for key in list(ds_out.keys()):    
    ds_out[key].encoding['zlib'] = True 

for key in ["y_post_it", "x_bc_it", "y_mod_ap", "y_obs",
            "sigma_y_it",  "acc_ratio", "sigma_y_ratio"]:
    ds_out[key].encoding['least_significant_digit'] = 3 

ds_out.to_netcdf(path=fname_out) 

# Write a separate script to calculate grid cell and country emissions and uncertainties
# Then plot up national emissions from that. and compare to EnKF results


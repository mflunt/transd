#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 09:48:42 2021

MOdule file to set up MCMC code

Basically to keep line size down from actual run script

@author: mlunt
"""
import numpy as np
import xarray
import matplotlib.pyplot as plt
import pandas as pd
import json
import time as run_time
import test_mod as mcmc_mod
import datetime as dt
from dateutil.relativedelta import relativedelta

def fn_multipledistances(origin, lat, lon):
    # Calculate the distance between the point of interest and 
    #  ALL points in the grid
    radius = 6371 #km
    lat0 = origin[0]
    lon0 = origin[1]
    
    distances = np.zeros(len(lat))    
            
    #lat2 = lat
    #lon2 = lon[j

    dlat = np.radians(lat-lat0)
    dlon = np.radians(lon-lon0)
    
    a = np.sin(dlat/2) * np.sin(dlat/2) + np.cos(np.radians(lat0)) \
        * np.cos(np.radians(lat)) * np.sin(dlon/2) * np.sin(dlon/2)
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    
    d_j = radius * c           
    distances = d_j
              
    return distances

def get_nsigma_y_space(lon_bin, lat_bin, lon, lat,
                  nmeasure, sigma_values):      
    """
    Defines the indices of the measurement vector which are described by different 
    uncertainty hyperparameters (sigma_model)
    By default all sigma_models are split by spatital position.
    Bins into which to divide measurements must be defined and passe dinto function.     
    """     
    y_bl=np.zeros((nmeasure))
    nsigma=0
    ydim1=0
    sigma_models=[]
    
    nlon_bins = len(lon_bin)-1
    nlat_bins = len(lat_bin)-1
   
    for loi in range(nlon_bins):
        for lai in range(nlat_bins):
            
            wh=np.where(np.logical_and(np.logical_and(lon>=lon_bin[loi],
                                lon<lon_bin[loi+1]),
                                np.logical_and(lat>=lat_bin[lai],
                                   lat<lat_bin[lai+1])))
        
            if len(wh[0]) > 0:
                y_bl[wh]=nsigma
                sigma_models.append(sigma_values)
                nsigma+=1
                
            n_obs = len(wh[0])
            if n_obs > ydim1:
                ydim1 = n_obs*1
                                 
    # INDEX R
    print (ydim1, nsigma)
    R_indices = np.zeros((ydim1,nsigma), dtype=np.uint16)
    for ii in range(nsigma):      
        wh_bl=np.where(y_bl == ii)
        nwh=len(wh_bl[0])
        R_indices[:nwh,ii]=wh_bl[0]+1
        if nwh < ydim1:
            R_indices[nwh:,ii]=np.max(wh_bl)+1
    
    ydim2=nsigma*1
    
    return R_indices, ydim1, ydim2   #, np.asarray(sigma_models)

def get_nsigma_y_time_space(lon_bin, lat_bin, lon, lat,
                  nmeasure, sigma_values,
                  y_time, start_date, end_date,
                  bl_period=1):      
    """
    Defines the indices of the measurement vector which are described by different 
    uncertainty hyperparameters (sigma_model)
    By default all sigma_models are split by spatital position.
    Bins into which to divide measurements must be defined and passe dinto function.     
    """  
    d0=start_date   # Need to be in pandas datetime format
    d1=end_date    # Needs to be in pandas datetime format
    delta = d1 - d0
    nmonths = int(np.round(delta.days/30.415))
    
    ntime_periods = np.int(np.ceil(nmonths/np.float(bl_period)))
    
    y_bl=np.zeros((nmeasure))
    nsigma=0
    ydim1=0
    sigma_models=[]
    
    nlon_bins = len(lon_bin)-1
    nlat_bins = len(lat_bin)-1
    #mf_time_temp = y_time
    mf_time_temp2=pd.to_datetime(y_time)
    bl_start=d0
    for ti in range(ntime_periods):
        bl_end=bl_start+relativedelta(months=bl_period)
        for loi in range(nlon_bins):
            for lai in range(nlat_bins):
            
                wh=np.where(np.logical_and(np.logical_and(np.logical_and(lon>=lon_bin[loi],
                                    lon<lon_bin[loi+1]),
                                    np.logical_and(lat>=lat_bin[lai],
                                       lat<lat_bin[lai+1])),
                                    np.logical_and(mf_time_temp2>=bl_start,
                                       mf_time_temp2<bl_end)))
            
                if len(wh[0]) > 0:
                    y_bl[wh]=nsigma
                    sigma_models.append(sigma_values)
                    nsigma+=1
                    
                #bl_start=bl_start+relativedelta(months=bl_period)
                n_obs = len(wh[0])
                if n_obs > ydim1:
                    ydim1 = n_obs*1
                    
        bl_start=bl_start+relativedelta(months=bl_period)
    # INDEX R
    print (ydim1, nsigma)
    #return y_bl, ydim1, nsigma
    R_indices = np.zeros((ydim1,nsigma), dtype=np.uint16)
    for ii in range(nsigma):      
        wh_bl=np.where(y_bl == ii)
        nwh=len(wh_bl[0])
        R_indices[:nwh,ii]=wh_bl[0]+1
        if nwh < ydim1:
            R_indices[nwh:,ii]=np.max(wh_bl)+1
    
    ydim2=nsigma*1
    
    return R_indices, ydim1, ydim2, np.asarray(sigma_models)

def get_nsigma_y(y_time, y_site, y_spc,start_date, end_date, nmeasure, 
                 bl_period,sites, species, nsites_spc):
      
    
    d0=pd.to_datetime(start_date)
    d1=pd.to_datetime(end_date)
    delta = d1 - d0
    ndays = delta.days + 1
    
    #bl_period = ndays/1.
    y_bl=np.zeros((nmeasure))
    
    nsigma=0
    nsigma_max = np.int(np.ceil(ndays/np.float(bl_period)))
    ntime_stn=np.zeros((nsites_spc))
    
    ydim1=0
    
    for spc in species:
    
        for si,site in enumerate(sites):
            
            wh_site_spc = np.where( (y_site==site) & (y_spc == spc))[0]
            if len(wh_site_spc >0):
                
                #mf_time_temp = y_time[y_site==site]
                mf_time_temp = y_time[wh_site_spc]
                
                #fp_data_H3 = fp_data_H[sites[si]].dropna("time", how="all")        
                nsigma_stn=0
                bl_start=d0
                #mf_time_temp=fp_data_H3.time.values
                #mf_time_temp=fp_data_H[sites[si]].time.values
                mf_time_temp2=pd.to_datetime(mf_time_temp)
                
                ntime_stn[si]=len(mf_time_temp)
                
                for ti in range(nsigma_max):
                    bl_end=bl_start+dt.timedelta(days=bl_period)
                    
                    wh=np.where(np.logical_and(mf_time_temp2>=bl_start,
                                               mf_time_temp2<bl_end))
                                                  
                    if len(wh[0]) > 0:
                        # This line is a problem - it assumes everything is sequential
                        #y_bl[wh+np.sum(ntime_stn[:si],dtype=np.uint16)]=nsigma_stn+nsigma
                        #nsigma_stn+=1
                        
                        y_bl[wh_site_spc[wh]] = nsigma_stn+nsigma
                        nsigma_stn+=1
                        
                    bl_start=bl_start+dt.timedelta(days=bl_period)
                    
                    n_obs = len(wh[0])
                    if n_obs > ydim1:
                        ydim1 = n_obs*1
            
                nsigma+=nsigma_stn
    
    # INDEX R
    #print(ydim1, nsigma)
    #return y_bl, ydim1,nsigma
    #return y_bl, ydim1, nsigma
    R_indices = np.zeros((ydim1,nsigma), dtype=np.uint16)
    for ii in range(nsigma):      
        wh_bl=np.where(y_bl == ii)
        nwh=len(wh_bl[0])
        R_indices[:nwh,ii]=wh_bl[0]+1
        if nwh < ydim1:
            R_indices[nwh:,ii]=np.max(wh_bl)+1
    
    ydim2=nsigma*1
    
    return R_indices, ydim1, ydim2 

#%%
def run_mcmc(start_date, end_date, nIC, nBC, ngroups, nstate0, nruns,
             bl_period, dates, nIt, nsub, burn_in,
             beta, x_ap, H, y_obs, y_std,y_date,y_mod_ap,y_lat,y_lon, 
             x_pdf_all, pdf_param1_pdf, pdf_param2_pdf, pdf_param1_ap, pdf_param2_ap, 
             pdf_p1_hparam1, pdf_p1_hparam2,
             pdf_p2_hparam1, pdf_p2_hparam2,
             stepsize, stepsize_pdf_p1, stepsize_pdf_p2,
             sigma_model_pdf, sigma_model_ap, sigma_model_hparam1_0, sigma_model_hparam2_0,
             stepsize_sig_y_0, rho,
             para_temp=False
             ): 
    """
    Function to define all arrays needed for MCMC from the input scalars
    
    Also calls mcmc_mod to do MCMC
    
    Just need to figure out how to account for different run_dates 
    and different species groups CO, Co2gee, reco ff etc. 
    
    Use a y_spc variable to separate out where species are the same.
    
    Sites needs to be a list of both sites[CO2] and sites[CO]
    """
    ########################################################
    # Define full x hyperparam vecotrs for size k
    
    nobs= len(y_obs)
    nstate = len(x_ap)
    nbeta = len(beta)
    #nsites = len(sites)
    n0_ap = y_mod_ap - y_obs # This has to be this way round to be consistent with fortran script
    
    #%%
    # Set up off-diagonal correlation structure
    ###############################################################
    # But need to work out how to do this for multiple species
    
    nobs_day_list = []
    #for ss, si in enumerate(sites):
    #    wh_site = np.ravel(np.where(y_site == si))
    #    nobs_site[ss]=len(wh_site)
        
    deltaspace=np.zeros((nobs,nobs))+1.e12

    
    for di, date in enumerate(dates):
        wh_day = np.where(y_date == date )[0]
            
        if len(wh_day) >0:
            nobs_day_list.append(len(wh_day))
            
            for whi in wh_day:
                deltaspace[whi,whi] = 0.
#                deltaspace_site = fn_multipledistances([y_lat[whi],y_lon[whi]], 
#                  y_lat[wh_day],y_lon[wh_day])
#                ##tdelta = np.absolute(y_time[wh_site]-y_time[whi]).astype('timedelta64[m]')
#                ##deltaspace_site=tdelta/np.timedelta64(1, 'h')
#                
#                ##deltatime_site[wh_site] = np.absolute(y_time[wh_site]-y_time[whi])
#                deltaspace[whi,wh_day] = deltaspace_site
                #deltaspace[wh_day,whi] = deltaspace_site

    nobs_day = np.asarray(nobs_day_list)
    nday_max=np.max(nobs_day)    
    
    ndays = len(nobs_day) # This should be legth of sites[CO2] + sites[CO]
    
    cum_nobs = np.zeros((ndays), dtype=np.int16) # Fine up to 32768 (2**15)
    for ti in range(1,ndays):
        cum_nobs[ti] = np.sum(nobs_day[:ti])
    
    
#    out_dict={"deltaspace": deltaspace,
#              "nobs_day": nobs_day,
#              "cum_nobs": cum_nobs}
#    return out_dic
        
    
    lon_bin = np.arange(20.,65,10.)
    lat_bin = np.arange(-10.,24.,8.)
    
    #lon_bin = np.arange(20.,65,40.)
    #lat_bin = np.arange(-10.,24.,32.)


    R_indices, ydim1, ydim2 = get_nsigma_y_space(lon_bin, lat_bin, y_lon,y_lat,                                                     
                      nobs, sigma_model_ap)     
    
    

    #R_indices, ydim1, ydim2, sigma_model0= get_nsigma_y_time_space(lon_bin, lat_bin, lon,lat,                                                     
    #              nmeasure, sigma_model_ap,
    #              y_time, pandas.to_datetime(start_date),
    #              pandas.to_datetime(end_date),bl_period=sigma_y_period)  
    
    
    
    # Define R indices
    #R_indices, ydim1, ydim2 = get_nsigma_y(y_time,y_date, start_date, end_date, 
    #                                   nobs, bl_period,dates, ndays)

    #out_dict={"R_indices": R_indices,
    #          "ydim1": ydim1,
    #          "ydim2": ydim2}
    #return out_dict

    ################################################
    # Define sigma model params
    sigma_model0=np.zeros((ydim2))
    sigma_model0[:]=sigma_model_ap  
    
    sigma_measure = y_std.copy()
    error_structure = sigma_measure*0.+1.
    
    sigma_model_hparam1 = np.zeros((ydim2))+sigma_model_hparam1_0
    sigma_model_hparam2 = np.zeros((ydim2))+sigma_model_hparam2_0
    stepsize_sigma_y = np.zeros((ydim2))  + stepsize_sig_y_0
    
    ########################################################
    # Set-up all variables for parallel tempering:
    pdf_param1 = np.zeros((nstate,nbeta))
    pdf_param2 = np.zeros((nstate,nbeta))
    x=np.zeros((nstate,nbeta))
    sigma_model = np.zeros((ydim2,nbeta))
    n0=np.zeros((nobs,nbeta))  
    for ib in range(nbeta):  
        x[:,ib]=x_ap.copy()  
        n0[:,ib] = n0_ap.copy()
        sigma_model[:,ib]=sigma_model0.copy()
        pdf_param1[:,ib] = pdf_param1_ap
        pdf_param2[:,ib] = pdf_param2_ap
    
    
    #x_pdf_all = np.zeros((nstate), dtype=np.int8) + x_pdf
    
    sigma_model_pdf_all = np.zeros((ydim2), dtype=np.int8) + sigma_model_pdf
    
    # Call mcmc script
    nit_sub=nIt/nsub
    if para_temp == True:
        para_temp_in = 1
    else:
        para_temp_in = 0
        
        
#    out_dict={"R_indices": R_indices,
#              "nsites_spc": nsites_spc,
#              "x_ap": x_ap,
#              "pdf_param1_ap": pdf_param1_ap,
#              "pdf_param2_ap": pdf_param2_ap}
#    return out_dict
        
   
    #%%
    # MCMC version
    print ("Starting MCMC...")
    
    startt = run_time.time()
    
    # y_corr options
    x_it, y_post_it, sigma_model_it, sigma_y_it, \
        n0T_out, pdf_param1_out, pdf_param2_out, accept, reject, \
        accept_sigma_y, reject_sigma_y, accept_swap, reject_swap, \
        tot_acc_x, tot_acc_p1, tot_acc_p2, tot_acc_sigma_y, \
        accept_all, reject_all, accept_sigma_y_all, reject_sigma_y_all = mcmc_mod.mcmc_corr.hbtdmcmc(beta, x, H,y_obs, n0, 
            pdf_param1, pdf_param2, sigma_model, sigma_measure,error_structure, 
            R_indices, sigma_model_hparam1, sigma_model_hparam2, stepsize_sigma_y, sigma_model_pdf_all, 
            rho, deltaspace,  
            nobs_day, cum_nobs, nday_max, 
            para_temp_in, x_pdf_all, burn_in, 
            pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, pdf_param2_pdf, 
            stepsize, stepsize_pdf_p1,stepsize_pdf_p2, nIt, nsub, nit_sub, 
            nbeta, nstate, nobs, ydim1, ydim2, ndays)    
        
    
    endt = run_time.time()
    print ("MCMC complete in ", (endt-startt), " seconds")
    
    ratio = 1.*accept/(accept+reject)
    ratio_sigma_y = accept_sigma_y/(accept_sigma_y+reject_sigma_y)
    #ratio_all = accept_all/(accept_all+reject_all)
    
    #ds_out = xarray.Dataset(#"x_it": (["nstate", "nIt"], x_it),
    #                        "y_it": (["nobs","nIt"], y_post_it)
                                     
    out_dict = {"x_it": x_it,
                "y_it": y_post_it,
                "sigma_model_it": sigma_model_it,
                "sigma_y_it": sigma_y_it,
                "pdf_param2_it": pdf_param2_out,
                "ratio": ratio,
                "ratio_sigma_y": ratio_sigma_y,               
                }                                 
    
    return out_dict

#%%

def run_transd_mcmc(nIC,nBC, k_ap, kmin, kmax, beta, x_agg,
                    y_obs, n0_ap, y_std, y_lat, y_lon,
                    H_bc, h_v, lats, lons,
                    nIt, nsub, burn_in,
                    pdf_param10, pdf_param20,
                    pdf_p1_hparam10, pdf_p1_hparam20,
                    pdf_p2_hparam10, pdf_p2_hparam20,
                    x_pdf0, pdf_param1_pdf, pdf_param2_pdf,    
             stepsize, stepsize_pdf_p1, stepsize_pdf_p2,
             sigma_model_pdf, sigma_model_ap, sigma_model_hparam1_0, sigma_model_hparam2_0,
             stepsize_sig_y_0, 
             stepsize_clat, stepsize_clon, stepsize_bd,
             para_temp=False, rjmcmc=False
             ): 
    """
    Function to define all arrays needed for transd MCMC from the input scalars
    
    Also calls mcmc_mod to do MCMC
    
    Will use uncorrelated code, so no need for deltatimes or distances.
    
                    
             stepsize, stepsize_pdf_p1, stepsize_pdf_p2,
    
    
    Input (integer) scalars:
        nIC = Number of IC terms
        nBC = Number of boundary condition terms
        k_ap = Prior number of aggregated nodes
        kmin = Min number of nodes
        kmax = Max number of nodes
        nIt = NUmber of iterations
        nsub = nth iteration to store in output
        burn_in = Number of burn in iterations
        sigma_model_pdf = Integer definining sigma model PDF (1=Uniform, 2= Gaussian, 3=Lognormal)
        x_pdf0 = Prior pdf of scale factor values (1=Uniform, 2= Gaussian, 3=Lognormal)
        pdf_param1_pdf = Prior pdf of hparam 1 (1=Uniform, 2= Gaussian, 3=Lognormal)
        pdf_param2_pdf = Prior pdf of hparam 2 (1=Uniform, 2= Gaussian, 3=Lognormal)
        
    Input (float) scalars:
        pdf_param10 = Prior parameter 1 (i.e mean)
        pdf_param20 = Prior parameter 2 (std dev)
        pdf_p1_hparam10
        pdf_p1_hparam20,
        pdf_p2_hparam10
        pdf_p2_hparam20
        stepsize = 
        stepsize_pdf_p1 = 
        stepsize_pdf_p2 = 
        sigma_model_ap = Initital value for sigma model
        sigma_model_hparam1_0 = Hyperparameter 1 for sigma model (i.e. mean or lower limit)
        sigma_model_hparam2_0 = Hyperparameter 2 for sigma model (i.e. std dev or upper limit)
        stepsize_sig_y_0 = Jummp size for sigma model parameters
        stepsize_clat = latitude jump size for moves
        stepsize_clon = Longitude jump size for moves
        stepsize_bd = Parameter jump size for birth/death
        
    
    Input arrays:
        beta (nbeta) = Tempering values
        x_agg (k_ap) = Prior scale factor (usually 1)
        y_obs (nobs) = Observations
        n0_ap (nobs) = y_mod_ap - y_obs (Needs to be this way round)
        y_std (nobs) = Measurement uncertainty
        y_lat (nobs) = Satellite observation latitudes
        y_lon (nobs)  = Satellite obs longitudes
        H_bc (nobs, nBC+nIC) = Boundary and initial condition sensitvity
        h_v (nobs, Ngrid) = Grid level emissions sensitivity
        lats (Ngrid) = Lat centre points of sensitivity matrix
        lons (Ngrid) = Lon centre points of sensitivity matrix
        
    Input booleans:
        para_temp = Do parallel teempering or not
        rjmcmc = Use reversible jump or not. If not will just use a fixed grid
        
    """
    nobs= len(y_obs)
    nbeta = len(beta)
    Ngrid = len(lats)
    nIC1 = nIC+nBC+1
    nIBC = nIC + nBC
    
    kICmax = kmax + nIBC 
    
    latmin = np.min(lats)
    latmax = np.max(lats)
    lonmin = np.min(lons)
    lonmax = np.max(lons)
    
    
    #n0_ap = y_mod_ap - y_obs # This has to be this way round to be consistent with fortran script
    
    lon_bin = np.arange(20.,65,10.)
    lat_bin = np.arange(-10.,24.,8.)
    
    #lon_bin = np.arange(20.,65,40.)
    #lat_bin = np.arange(-10.,24.,32.)

    R_indices, ydim1, ydim2 = get_nsigma_y_space(lon_bin, lat_bin, y_lon,y_lat,                                                     
                      nobs, sigma_model_ap)     
    

    h_agg0 = np.zeros((nobs, k_ap+nIBC))
    h_agg0[:, :nIBC] = H_bc

    plat,plon, h_agg, regions_v = aggregate_H(kmax, kICmax, nIBC, nbeta, Ngrid, nobs,  k_ap,
                h_v, lats, lons, h_agg0)

    ################################################
    # Define sigma model params
    sigma_model0=np.zeros((ydim2))
    sigma_model0[:]=sigma_model_ap  
    
    sigma_measure = y_std.copy()
    
    sigma_model_hparam1 = np.zeros((ydim2))+sigma_model_hparam1_0
    sigma_model_hparam2 = np.zeros((ydim2))+sigma_model_hparam2_0
    #stepsize_sigma_y = np.zeros((ydim2))  + stepsize_sig_y_0
    stepsize_sigma_y_all=sigma_model0*stepsize_sig_y_0
    ########################################################
    # Set-up all variables for parallel tempering:
    

    sigma_model = np.zeros((ydim2,nbeta))
    n0=np.zeros((nobs,nbeta))  
    k=np.zeros((nbeta),dtype=np.int)+k_ap
    
#    for ib in range(nbeta):  
#        x[:,ib]=x_ap.copy()  
#        n0[:,ib] = n0_ap.copy()
#        sigma_model[:,ib]=sigma_model0.copy()
#        pdf_param1[:,ib] = pdf_param1_ap
#        pdf_param2[:,ib] = pdf_param2_ap
        
    x=np.zeros((kICmax,nbeta))
    #sigma_y = np.zeros((nmeasure,nbeta))
    sigma_model = np.zeros((ydim2,nbeta))

    for ib in range(nbeta):  
        x[:k_ap+nIBC,ib]=x_agg.copy()
        sigma_model[:,ib]=sigma_model0.copy()
        n0[:,ib] = n0_ap.copy()
    
    
    #x_pdf_all = np.zeros((nstate), dtype=np.int8) + x_pdf
    
    #sigma_model_pdf_all = np.zeros((ydim2), dtype=np.int8) + sigma_model_pdf
    #%%
    # Call mcmc script
    nit_sub=nIt/nsub
    if para_temp == True:
        para_temp_in = 1
    else:
        para_temp_in = 0
        
    if rjmcmc == True:
        rjmcmc_in = 1
    else:
        rjmcmc_in = 0
        
    pdf_param1 = np.zeros((kICmax,nbeta))
    pdf_param2 = np.zeros((kICmax,nbeta))

    #########################################
    sigma_clon = stepsize_clon*1.
    sigma_clat = stepsize_clat*1.
    sigma_bd=np.mean(x_agg[nIBC:])*stepsize_bd
    
    ################################################
    # TUNING OF INDIVIDUAL PARAMETER STEPSIZES AND UNCERTAINTIES
    
    stepsize_all=np.zeros((nIBC+1))+stepsize
    stepsize_pdf_p1_all=np.zeros((nIBC+1))+(stepsize_pdf_p1)
    stepsize_pdf_p2_all=np.zeros((nIBC+1))+(stepsize_pdf_p2)
    stepsize_all[:nIBC]=stepsize_all[:nIBC]/200.
    stepsize_pdf_p1_all[:nIBC]=stepsize_pdf_p1_all[:nIBC]/10.
    stepsize_pdf_p2_all[:nIBC]=stepsize_pdf_p2_all[:nIBC]/10.
    
    stepsize_all[-1]=stepsize_all[-1]*2.
    stepsize_pdf_p1_all[-1]=stepsize_pdf_p1_all[-1]/5.
    stepsize_pdf_p2_all[-1]=stepsize_pdf_p2_all[-1]
    
    stepsize_all[:2]=stepsize_all[:2]*10.
    stepsize_all[2]=stepsize_all[2]*100.
    
    
    pdf_param1[:,:]=pdf_param10
    pdf_param2[:,:]=pdf_param20
    
    pdf_p1_hparam1=np.zeros((nIC1))+pdf_p1_hparam10
    pdf_p1_hparam2=np.zeros((nIC1))+pdf_p1_hparam20
    
    pdf_p2_hparam1=np.zeros((nIC1))+pdf_p2_hparam10
    pdf_p2_hparam2=np.zeros((nIC1))+pdf_p2_hparam20
    
    x_pdf = np.zeros((nIC1), dtype=np.int)+x_pdf0
    #pdf_param1_pdf=np.zeros((nIC1), dtype=np.int)+pdf_param1_pdf0
    #pdf_param2_pdf=np.zeros((nIC1), dtype=np.int)+pdf_param2_pdf0
    
    #pdf_param1_pdf=pdf_param1_pdf0
    #pdf_param2_pdf=pdf_param2_pdf0
    
    
    #out_dict={"plat": plat,
    #          "plon": plon,
    #          "h_agg": h_agg,
    #          "regions_v": regions_v}
    #return out_dict
    
     # MCMC version
     
    import tdmcmc_mod
    print ("Starting TransD MCMC...")
    
    startt = run_time.time()
    
    k_it, x_out, y_out, regions_out, plon_out, plat_out, sigma_model_out,sigma_y_out, \
    n0T_out,pdf_param1_out,pdf_param2_out, accept, reject, \
    accept_birth, reject_birth, accept_death, reject_death, accept_move, reject_move, \
    accept_sigma_y, reject_sigma_y, accept_swap, reject_swap, \
    stepsize_x_out, stepsize_p1_out, stepsize_p2_out, \
    stepsize_sigma_y_out, accept_all, reject_all, accept_birth_all, reject_birth_all, \
    accept_death_all, reject_death_all, accept_move_all, reject_move_all, \
    accept_sigma_y_all, reject_sigma_y_all = tdmcmc_mod.mcmc_uncorr.hbtdmcmc(      
    beta,k, x, h_agg,y_obs,n0, plon, plat, regions_v, 
    pdf_param1, pdf_param2, lons,lats, h_v, sigma_model, sigma_measure, 
    R_indices, sigma_model_hparam1, sigma_model_hparam2,
    stepsize_sigma_y_all, sigma_model_pdf, 
    sigma_clon, sigma_clat, rjmcmc_in, para_temp_in,
    lonmin, lonmax, latmin,latmax, sigma_bd, kmin, x_pdf, burn_in, 
    pdf_p1_hparam1, pdf_p1_hparam2, pdf_p2_hparam1, pdf_p2_hparam2, pdf_param1_pdf, 
    pdf_param2_pdf,stepsize_all, stepsize_pdf_p1_all,stepsize_pdf_p2_all, 
    nIt, nsub, nit_sub, nIBC, 
    nbeta, kmax, kICmax, nobs, Ngrid, ydim1, ydim2, nIC1)
      
    endt = run_time.time()
    print ("TransD MCMC complete in ", (endt-startt), " seconds")
    
    ratio = 1.*accept/(accept+reject)
    ratio_sigma_y = accept_sigma_y/(accept_sigma_y+reject_sigma_y)
    
    if rjmcmc == True:
        ratio_birth = accept_birth/(accept_birth+reject_birth)
        ratio_death = accept_death/(accept_death+reject_death)
        ratio_move = accept_move/(accept_move+reject_move)
    else:
        ratio_birth=0.
        ratio_death=0.
        ratio_move=0.
        
    if para_temp == True:
        ratio_swap = accept_swap/(accept_swap+reject_swap)
    else:
        ratio_swap=0
    #ratio_all = accept_all/(accept_all+reject_all)
    
    #ds_out = xarray.Dataset(#"x_it": (["nstate", "nIt"], x_it),
    #                        "y_it": (["nobs","nIt"], y_post_it)
                                     
    out_dict = {"x_it": x_out,
                "k_it": k_it,
                "regions_it": regions_out,
                "y_it": y_out,   # Needt to add this
                "sigma_model_it": sigma_model_out,
                "sigma_y_it": sigma_y_out,
                "pdf_param2_it": pdf_param2_out,
                "ratio": ratio,
                "ratio_sigma_y": ratio_sigma_y, 
                "ratio_birth": ratio_birth,
                "ratio_death": ratio_death,
                "ratio_move": ratio_move,
                "ratio_swap": ratio_swap,
                }                                 
    
    return out_dict
    
    
      
    
def aggregate_H(kmax, kICmax, nIC, nbeta, Ngrid, nobs,  k_ap,
                h_v, lats, lons, h_agg0):
    # Define prior model and regions with uniform distribution
    #######################################
   
    # Set up different starting nuclei locations for each chain 
    plon=np.zeros((kmax,nbeta))
    plat=np.zeros((kmax,nbeta))
    regions_v=np.zeros((Ngrid,nbeta),dtype=np.uint16)
    h_agg=np.zeros((nobs, kICmax,nbeta))
    #n0=np.zeros((nobs,nbeta))    
    

    for ib in range(nbeta):
        zgrid=np.random.randint(0,Ngrid,k_ap)
        #zlat=np.random.randint(0,nlat,k_ap)
        plon[:k_ap,ib] = lons[zgrid] # Lon locs of nuclei
        plat[:k_ap,ib] = lats[zgrid] # Lat locs of nuclei
        
        
        region = np.zeros((Ngrid), dtype=np.uint16)
        regions_v0=closest_grid(region, lats, lons, plat[:k_ap,ib], plon[:k_ap,ib], \
                np.arange(0, k_ap, dtype=np.uint16))
        #regions_v0 = np.ravel(regions0)
        regions_v[:,ib]=regions_v0.copy()+1
    
        for ri in range(k_ap):
            wh_ri = np.where(regions_v0 == ri)
            for ti in range(nobs):
                h_agg0[ti,ri+nIC]=np.sum(h_v[ti,wh_ri])
      
        #y_model = np.dot(h_agg0,x_agg) 
        #n0_ap = y_model-y
       
        h_agg[:,:k_ap+nIC,ib] = h_agg0.copy()
        #n0[:,ib]=n0_ap.copy()

    return plat,plon, h_agg, regions_v

def closest_grid(region_v, lats, lons, plat, plon, pind):
    """
    Find closest lats and lons to a given lat lon coordinate.
    Solve as a 1D problem
    
    """
    ngrid = len(lats)
          
    lai=0
    for idx in range(ngrid):
   
        maxdist=1e6
        for pi in pind: 
            dist=(lats[idx] - plat[pi])*(lats[idx] - plat[pi]) + (lons[idx] - plon[pi])*(lons[idx] - plon[pi])  
           #! Can't work out which way round??????????????
            if dist < maxdist:
                region_v[lai]=pi
                maxdist=dist                  
        lai+=1
   
    return region_v


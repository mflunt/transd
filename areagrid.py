#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 14:21:26 2018

Areagrid: Calculate a grid of areas (m2) given arrays of 
    latitudes and longitudes

    Example: 
    
    import areagrid
    lat=np.arange(50., 60., 1.)
    lon=np.arange(0., 10., 1.)
    
    area=areagrid(lat, lon)
@author: mlunt
"""

import numpy as np

def areagrid(lat, lon):

  re=6367500.0	#radius of Earth in m
  
  dlon=abs(np.mean(lon[1:] - lon[0:-1]))*np.pi/180.
  dlat=abs(np.mean(lat[1:] - lat[0:-1]))*np.pi/180.
  theta=np.pi*(90.-lat)/180.
  
  area=np.zeros((len(lat), len(lon)))
  
  for latI in range(len(lat)):
    if theta[latI] == 0. or theta[latI] == np.pi:
      area[latI, :]=(re**2)*abs(np.cos(dlat/2.)-np.cos(0.))*dlon
    else:
      lat1=theta[latI] - dlat/2.
      lat2=theta[latI] + dlat/2.
      area[latI, :]=((re**2)*(np.cos(lat1)-np.cos(lat2))*dlon)

  return area
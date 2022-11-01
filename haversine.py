#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 16:03:08 2014
 
FUNCTIONS
___________ 
 
 
distance


Calculates the distance in km between two lat/lon pairs using the haversine formula.


Example:

distance = acrg_grid.haversine.distance(origin, destination, radius=radius)


Inputs:

origin: [lat_1, lon_1] two element array containing the lat and lon of your original point
destination: [lat_2, lon_2] two element array containing the lat and lon of your destination point
radius: radius in km. Defaults to 6371km.


Output:

distance between the two points in km



CLASSES
_________


multipledistances

Calculates the distance in km between a single lat/lon pair and a grid of lat/lons using the haversine formula.
Also identifies the minimum distance and it's location.


Example:

distances = haversine.multipledistances(origin, lat, lon, radius=radius)

Inputs:

origin: [lat_1, lon_1] two element array containing the lat and lon of your original point
lat: an array of latitudes of length n
lon: an array of longitudes of lenth m
radius: radius in km. Defaults to 6371km.

Outputs:

distances= an n by m array containing the distances from every point in the nxm grid to the origin
mindist = the minimum distance
mindist_index = two element array containing the index in n and m corresponding to the minium distance
e.g. if the minimum distance is at the 2nd latitude and 17th longitude then mindist_index = [2,17]
mindist_loc = two element array containing the lat and lon corresponding to the minimum distance


"""
 
import math
import numpy as np
#from numba import jit
#import pdb

radius = 6371 #km
 
def distance(origin, destination, radius=radius):
    lat1, lon1 = origin
    lat2, lon2 = destination

    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
 
    return d

def multiple_distances(origin, lats,lons):
    
    lat0 = origin[0]
    lon0 = origin[1]
    
    dlats = np.radians(lats-lat0)
    dlons = np.radians(lons-lon0)
    
    a = np.sin(dlats/2) * np.sin(dlats/2) + np.cos(np.radians(lat0)) \
            * np.cos(np.radians(lat0)) * np.sin(dlons/2) * np.sin(dlons/2)
            
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
        
    distances = radius * c        
    
    return distances

    

    


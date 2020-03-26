#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experiment.
"""

# setup
import os
import sys
alp = os.path.abspath('../LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import zfun

import netCDF4 as nc4
import numpy as np

def box_inds(lon0,lon1,lat0,lat1,grid_fn):

    # label output dict
    D = dict()

    # get Datasets
    ds = nc4.Dataset(grid_fn)

    # gather some fields, for convenience
    lon_rho = ds['lon_rho'][0,:].squeeze()
    lat_rho = ds['lat_rho'][:,0].squeeze()
    lon_u = ds['lon_u'][0,:].squeeze()
    lat_u = ds['lat_u'][:,0].squeeze()
    lon_v = ds['lon_v'][0,:].squeeze()
    lat_v = ds['lat_v'][:,0].squeeze()

    # find rho index rel. to reference point, find u and v indeces rel. to rho point
    # relation of u and v pts to rho pt determined from the definitions on the
    # ROMS grid here: https://www.myroms.org/wiki/Numerical_Solution_Technique
    lonr0 = zfun.find_nearest_ind(lon_rho,lon0)
    D['lonr0'] = lonr0
    D['lonu0'] = lonr0-1
    D['lonv0'] = lonr0

    lonr1 = zfun.find_nearest_ind(lon_rho,lon1)
    D['lonr1'] = lonr1
    D['lonu1'] = lonr1
    D['lonv1'] = lonr1

    latr0 = zfun.find_nearest_ind(lat_rho,lat0)
    D['latr0'] = latr0
    D['latu0'] = latr0
    D['latv0'] = latr0-1

    latr1 = zfun.find_nearest_ind(lat_rho,lat1)
    D['latr1'] = latr1
    D['latu1'] = latr1
    D['latv1'] = latr1

    ds.close()
    
    return D


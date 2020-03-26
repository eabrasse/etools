#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot results of a particle tracking experiment.
"""

# setup
import os
import sys
alp = os.path.abspath('/pmr4/eab32/LiveOcean/alpha')
if alp not in sys.path:
    sys.path.append(alp)
import Lfun
import zrfun
import zfun
import netCDF4 as nc
import pickle
import numpy as np

outdir = 'momentum_balance'
Lfun.make_dir(outdir,clean=True)

#get ocean background
dir0 = '/pmr4/eab32/LiveOcean_ROMS/output/aestus1_base_ae1/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

NT=len(f_list)

fn_list = os.listdir(dir0 + f_list[0])
fna = [x for x in fn_list if x[-11:-8]=='avg']
fnd = [x for x in fn_list if x[-11:-8]=='dia']

oceanfna = dir0 + f_list[0] + '/' + fna[0]

dsa = nc.Dataset(oceanfna)
salt = dsa['salt'][0,:,:,:]
NZ, NY, NX = np.shape(salt)
G,S,T = zrfun.get_basic_info(oceanfna)
H = G['h'][:]
# lonr = G['lon_rho'][:]
# lonu = G['lon_u'][:]
# lonv = G['lon_v'][:]
# latr = G['lat_rho'][:]
# latu = G['lat_u'][:]
# latv = G['lat_v'][:]
maskr = G['mask_rho'][:]
masku = G['mask_u'][:]
maskv = G['mask_v'][:]

count = 0

#initialize variable holders
AVG_vars = {'salt' : np.zeros((NT, NY, NX)),
            'z' : np.zeros((NT, NZ, NY, NX)),
            'lonr' : G['lon_rho'][:],
            'lonu' : G['lon_u'][:],
            'lonv' : G['lon_v'][:],
            'latr' : G['lat_rho'][:],
            'latu' : G['lat_u'][:],
            'latv' : G['lat_v'][:]}
DIA_vars = {'u_hadv' : np.zeros((NT, NZ, NY, NX-1)),
            'u_accel' : np.zeros((NT,NZ,NY,NX-1)),
            'u_vadv' : np.zeros((NT, NZ, NY, NX-1)),
            'u_vvisc' : np.zeros((NT, NZ, NY, NX-1)),
            'u_cor' : np.zeros((NT, NZ, NY, NX-1)),
            'u_prsgrd' : np.zeros((NT, NZ, NY, NX-1)),
            'v_hadv' : np.zeros((NT, NZ, NY-1, NX)),
            'v_accel' : np.zeros((NT,NZ,NY-1,NX)),
            'v_vadv' : np.zeros((NT, NZ, NY-1, NX)),
            'v_vvisc' : np.zeros((NT, NZ, NY-1, NX)),
            'v_cor' : np.zeros((NT, NZ, NY-1, NX)),
            'v_prsgrd' : np.zeros((NT, NZ, NY-1, NX))}

for i in range(NT):
    print('%i/%i' % (i,NT))
    #identify file names for that day
    fn_list = os.listdir(dir0 + f_list[i])
    
    fnd = [x for x in fn_list if x[-11:-8]=='dia']
    fna = [x for x in fn_list if x[-11:-8]=='avg']
    
    #load in data for each hour and find the average
    for h in np.arange(0,24):
    
        oceanfnd = dir0 + f_list[i] + '/' + fnd[h]
        oceanfna = dir0 + f_list[i] + '/' + fna[h]
    
        dsd = nc.Dataset(oceanfnd)
        dsa = nc.Dataset(oceanfna)
        
        G,S,T = zrfun.get_basic_info(oceanfna)
        zeta = dsa['zeta'][:]
        salt = dsa['salt'][0,-1,:,:]
        zr = zrfun.get_z(H,zeta,S,only_rho=True)

        D = dict()

        if h==0:
            for key in dsd.variables.keys():
                if key in DIA_vars.keys():
                    DIA_vars[key][i,:,:,:] = dsd[key][:]
            z = zr
            sl = np.ma.masked_where(maskr==0, salt)

        else:
            for key in dsd.variables.keys():
                if key in DIA_vars.keys():
                    DIA_vars[key][i,:,:,:] = DIA_vars[key][i,:,:,:]+dsd[key][:]
            z = z + zr
            sl = sl + np.ma.masked_where(maskr==0, salt)
        
        dsd.close()
        dsa.close()
    
    #take the average
    for key in DIA_vars:
        DIA_vars[key][i,:,:,:] = DIA_vars[key][i,:,:,:]/24
    AVG_vars['z'][i,:,:,:] = z/24
    AVG_vars['salt'][i,:,:] = sl/24

outfna = outdir+'AVG_vars.p'
outfnd = outdir+'DIA_vars.p'
pickle.dump(AVG_vars,open(outfna,'wb'))
pickle.dump(DIA_vars,open(outfnd,'wb'))

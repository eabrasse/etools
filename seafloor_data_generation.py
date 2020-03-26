#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a filtered time series of a lat-lon slice of ROMS output
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
import netCDF4 as nc4
import numpy as np

import pickle

#get ocean background
dir0 = '/pmr4/eab32/LiveOcean_ROMS/output/aestus1_base_ae1/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

f_list=f_list[:]
ndays = len(f_list)
nt = ndays*24

z_ind = 5
tt = 0
for t in range(ndays):
    print('day %i' % t)
    
    fn_list = os.listdir(dir0 + f_list[t])
    fn_list.sort()
    fnd_list = [x for x in fn_list if x[6:9]=='dia']
    fna_list = [x for x in fn_list if x[6:9]=='avg']
    
    hrs = len(fnd_list)
    for hr in range(hrs):
        
        oceanfna = dir0 + f_list[t] + '/' + fna_list[hr]
        oceanfnd = dir0 + f_list[t] + '/' + fnd_list[hr]
    
        dsa = nc4.Dataset(oceanfna)
        dsd = nc4.Dataset(oceanfnd)
        
        #read in variable slices from data
        rho = dsa['rho'][:]
        density = rho[0,z_ind,:,:].squeeze()

        u_prsgrd = dsd['u_prsgrd'][:]
        u_prsgrd = 0.5*(u_prsgrd[0,z_ind,1:-1,1:] + u_prsgrd[0,z_ind,1:-1,:-1])
        v_prsgrd = dsd['v_prsgrd'][:]
        v_prsgrd = 0.5*(v_prsgrd[0,z_ind,1:,1:-1] + v_prsgrd[0,z_ind,:-1,1:-1])
    
        u_cor = dsd['u_cor'][:]
        u_cor = 0.5*(u_cor[0,z_ind,1:-1,1:] + u_cor[0,z_ind,1:-1,:-1])
        v_cor = dsd['v_cor'][:]
        v_cor = 0.5*(v_cor[0,z_ind,1:,1:-1] + v_cor[0,z_ind,:-1,1:-1])
    
        u_ageo = u_prsgrd + u_cor
        v_ageo = v_prsgrd + v_cor
    
        if tt==0:
            NXr, NYr = np.shape(density)
            NXu, NYu = np.shape(u_prsgrd)
            NXv, NYv = np.shape(v_prsgrd)
            #initialize matrices
            density_mat = np.zeros((nt, NXr, NYr))
            u_prsgrd_mat = np.zeros((nt, NXu, NYu))
            u_ageo_mat = np.zeros((nt, NXu, NYu))
            v_prsgrd_mat = np.zeros((nt, NXv, NYv))
            v_ageo_mat = np.zeros((nt, NXv, NYv))
    
        #dump data slices into filterable matrices
        density_mat[tt,:,:] = density
        u_prsgrd_mat[tt,:,:] = u_prsgrd
        u_ageo_mat[tt,:,:] = u_ageo
        v_prsgrd_mat[tt,:,:] = v_prsgrd
        v_ageo_mat[tt,:,:] = v_ageo
    
        dsa.close()
        dsd.close()
        tt = tt+1
    
#filter time series
print('filtering...')
days = range(35,ndays*24-35,24)
Ndays = len(days)
density_gf = np.zeros((Ndays,NXr,NYr))
u_prsgrd_gf= np.zeros((Ndays,NXu,NYu))
u_ageo_gf = np.zeros((Ndays,NXu,NYu))
v_prsgrd_gf = np.zeros((Ndays,NXv,NYv))
v_ageo_gf = np.zeros((Ndays,NXv,NYv))

filt = zfun.godin_shape()
filt = filt / filt.sum()
filt = np.reshape(filt,(np.shape(filt)[0],1,1))
filt_r = np.tile(filt,(1,NXr,NYr))
filt_u = np.tile(filt,(1,NXu,NYu))
filt_v = np.tile(filt,(1,NXv,NYv))

for d in range(Ndays):
    t1 = days[d]
    t0 = t1-35
    t2 = t1+35+1
    density_gf[d,:,:] = np.sum(filt_r*density_mat[t0:t2,:,:],axis=0)
    u_prsgrd_gf[d,:,:] = np.sum(filt_u*u_prsgrd_mat[t0:t2,:,:],axis=0)
    u_ageo_gf[d,:,:] = np.sum(filt_u*u_ageo_mat[t0:t2,:,:],axis=0)
    v_prsgrd_gf[d,:,:] = np.sum(filt_v*v_prsgrd_mat[t0:t2,:,:],axis=0)
    v_ageo_gf[d,:,:] = np.sum(filt_v*v_ageo_mat[t0:t2,:,:],axis=0)
    
D = dict()
D['density_gf'] = density_gf
D['u_prsgrd_gf'] = u_prsgrd_gf
D['u_ageo_gf'] = u_ageo_gf
D['v_prsgrd_gf'] = v_prsgrd_gf
D['v_ageo_gf'] = v_ageo_gf

out_fn = '/pmr4/eab32/etools/seafloor_data.p'
pickle.dump(D, open(out_fn, 'wb'))
print('saving to %s' % out_fn)

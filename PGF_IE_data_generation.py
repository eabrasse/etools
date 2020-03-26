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
dir0 = '/pmr4/eab32/LiveOcean_ROMS/output/aestus2_unstrat_ae1/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

f_list=f_list[:90]
ndays = len(f_list)
nt = ndays*24

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
        density = rho[0,-1,:,:].squeeze()
        salt = dsa['salt'][0,-1,:,:].squeeze()

        u_prsgrd = dsd['u_prsgrd'][:]
        u_prsgrd_ext = 0.5*(u_prsgrd[0,-1,1:-1,1:] + u_prsgrd[0,-1,1:-1,:-1])
        u_prsgrd_int = 0.5*(u_prsgrd[0,0,1:-1,1:] + u_prsgrd[0,0,1:-1,:-1])
        v_prsgrd = dsd['v_prsgrd'][:]
        v_prsgrd_ext = 0.5*(v_prsgrd[0,-1,1:,1:-1] + v_prsgrd[0,-1,:-1,1:-1])
        v_prsgrd_int = 0.5*(v_prsgrd[0,0,1:,1:-1] + v_prsgrd[0,0,:-1,1:-1])
    
        if tt==0:
            NXr, NYr = np.shape(density)
            NXu, NYu = np.shape(u_prsgrd_ext)
            NXv, NYv = np.shape(v_prsgrd_ext)
            #initialize matrices
            density_mat = np.zeros((nt, NXr, NYr))
            salinity_mat = np.zeros((nt, NXr, NYr))
            u_prsgrd_ext_mat = np.zeros((nt, NXu, NYu))
            u_prsgrd_int_mat = np.zeros((nt, NXu, NYu))
            v_prsgrd_ext_mat = np.zeros((nt, NXv, NYv))
            v_prsgrd_int_mat = np.zeros((nt, NXv, NYv))
    
        #dump data slices into filterable matrices
        density_mat[tt,:,:] = density
        salinity_mat[tt,:,:] = salt
        u_prsgrd_ext_mat[tt,:,:] = u_prsgrd_ext
        u_prsgrd_int_mat[tt,:,:] = u_prsgrd_int
        v_prsgrd_ext_mat[tt,:,:] = v_prsgrd_ext
        v_prsgrd_int_mat[tt,:,:] = v_prsgrd_int
    
        dsa.close()
        dsd.close()
        tt = tt+1
    
#filter time series
print('filtering...')
days = range(35,ndays*24-35,24)
Ndays = len(days)
density_gf = np.zeros((Ndays,NXr,NYr))
salinity_gf = np.zeros((Ndays, NXr, NYr))
u_prsgrd_ext_gf= np.zeros((Ndays,NXu,NYu))
u_prsgrd_int_gf = np.zeros((Ndays,NXu,NYu))
v_prsgrd_ext_gf = np.zeros((Ndays,NXv,NYv))
v_prsgrd_int_gf = np.zeros((Ndays,NXv,NYv))

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
    salinity_gf[d,:,:] = np.sum(filt_r*salinity_mat[t0:t2,:,:],axis=0)
    u_prsgrd_ext_gf[d,:,:] = np.sum(filt_u*u_prsgrd_ext_mat[t0:t2,:,:],axis=0)
    u_prsgrd_int_gf[d,:,:] = np.sum(filt_u*u_prsgrd_int_mat[t0:t2,:,:],axis=0)
    v_prsgrd_ext_gf[d,:,:] = np.sum(filt_v*v_prsgrd_ext_mat[t0:t2,:,:],axis=0)
    v_prsgrd_int_gf[d,:,:] = np.sum(filt_v*v_prsgrd_int_mat[t0:t2,:,:],axis=0)

u_prsgrd_diff = u_prsgrd_int_gf - u_prsgrd_ext_gf
v_prsgrd_diff = v_prsgrd_int_gf - v_prsgrd_ext_gf
D = dict()
D['density_gf'] = density_gf
D['salinity_gf'] = salinity_gf
D['u_prsgrd_ext_gf'] = u_prsgrd_ext_gf
D['u_prsgrd_int_gf'] = u_prsgrd_int_gf
D['u_prsgrd_diff'] = u_prsgrd_diff
D['v_prsgrd_ext_gf'] = v_prsgrd_ext_gf
D['v_prsgrd_int_gf'] = v_prsgrd_int_gf
D['v_prsgrd_diff'] = v_prsgrd_diff

out_fn = '/pmr4/eab32/etools/PGF_IE_data_aestus2_unstrat.p'
pickle.dump(D, open(out_fn, 'wb'))
print('saving to %s' % out_fn)

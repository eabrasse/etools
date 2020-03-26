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
dir0 = '/pmr4/eab32/LiveOcean_ROMS/output/aestus2_base_ae1/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

f_list=f_list[:90]
ndays = len(f_list)
nt = ndays*24

z_ind = -1
tt = 0
for t in range(ndays):
    print('day %i' % t)
    
    fn_list = os.listdir(dir0 + f_list[t])
    fn_list.sort()
    fna_list = [x for x in fn_list if x[6:9]=='avg']
    
    hrs = len(fna_list)
    for hr in range(hrs):
        
        oceanfna = dir0 + f_list[t] + '/' + fna_list[hr]
    
        dsa = nc4.Dataset(oceanfna)
        
        #read in variable slices from data
        salt0 = dsa['salt'][:]
        salt = salt0[0,z_ind,:,:].squeeze()
    
        if tt==0:
            NXr, NYr = np.shape(salt)

            #initialize matrices
            salt_mat = np.zeros((nt, NXr, NYr))
    
        #dump data slices into filterable matrices
        salt_mat[tt,:,:] = salt

    
        dsa.close()
        tt = tt+1
    
#filter time series
print('filtering...')
days = range(35,ndays*24-35,24)
Ndays = len(days)
salt_gf = np.zeros((Ndays,NXr,NYr))

filt = zfun.godin_shape()
filt = filt / filt.sum()
filt = np.reshape(filt,(np.shape(filt)[0],1,1))
filt_r = np.tile(filt,(1,NXr,NYr))

for d in range(Ndays):
    t1 = days[d]
    t0 = t1-35
    t2 = t1+35+1
    salt_gf[d,:,:] = np.sum(filt_r*salt_mat[t0:t2,:,:],axis=0)
    
D = dict()
D['salt_gf'] = salt_gf

out_fn = '/pmr4/eab32/etools/salinity_data_aestus2.p'
pickle.dump(D, open(out_fn, 'wb'))
print('saving to %s' % out_fn)

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

import numpy as np
import netCDF4 as nc
import pickle

def get_delta(fn,pmpn,z,y,x):
    ds = nc.Dataset(fn)
    
    #delta is a small pain to calculate :/
    om = 1/ds[pmpn][y,x]
    # om0 = 0.5*(om00[:,:-1]+om00[:,1:])
    # om = np.tile(om0,(NZ,1,1))
    zeta = ds['zeta'][:]
    h = ds['h'][:]
    S = zrfun.get_basic_info(fn,only_S=True)
    zw = zrfun.get_z(h, zeta, S, only_w=True)
    Hz0 = zw[1:,:,:]-zw[:-1,:,:]
    
    if pmpn=='pn': #u grid
        Hz = 0.5*(Hz0[z,y,x+1]+Hz0[z,y,x])
    elif pmpn=='pm':
        Hz = 0.5*(Hz0[z,y,x]+Hz0[z,y+1,x])
    
    delta = Hz*om #delta0
    ds.close()
    return delta

# IMPORTANT CHOICES
varchoice = 'v'
latchoice = 45.25
lonchoice = -0.02


pmpn_d = {'u':'pn','v':'pm'}
pmpn = pmpn_d[varchoice]
lon_str = 'lon_' + varchoice
lat_str = 'lat_' + varchoice
trans_str = 'H' + varchoice + 'o' + pmpn[1]
accel_str = varchoice + '_accel'

dir0 = '/pmr4/eab32/LiveOcean_ROMS/output/aestus1_base_ae1/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

start_day = 30
end_day = 90

f_list = f_list[start_day:end_day+1]

NT=len(f_list)

series_len = NT*24

var_accel_old = np.zeros((series_len))
Ddvardt_oD = np.zeros((series_len))
vardDdt_oD = np.zeros((series_len))
var_accel_new = np.zeros((series_len))
ot = np.zeros((series_len))

t=0

for i in range(NT):
    print('%i/%i' % (i,NT))
    #identify file names for that day
    fn_list = os.listdir(dir0 + f_list[i])
    
    fnh = [x for x in fn_list if x[-11:-8]=='his']
    fnd = [x for x in fn_list if x[-11:-8]=='dia']
    fna = [x for x in fn_list if x[-11:-8]=='avg']
    
    #load in data for each hour and find the average
    for h in np.arange(0,24):
        
        oceanfnh0 = dir0 + f_list[i] + '/' + fnh[h]
        oceanfnh1 = dir0 + f_list[i] + '/' + fnh[h+1]   
        oceanfnd = dir0 + f_list[i] + '/' + fnd[h]
        oceanfna = dir0 + f_list[i] + '/' + fna[h]
        
        ds_his0 = nc.Dataset(oceanfnh0)
        ds_his1 = nc.Dataset(oceanfnh1)
        ds_dia = nc.Dataset(oceanfnd)
        ds_avg = nc.Dataset(oceanfna)
        
        if t==0:
            dt = ds_his1['ocean_time'][:]-ds_his0['ocean_time'][:] #dt
            lon_var = ds_his0[lon_str]
            lon_vec = lon_var[0,:]
            lat_var = ds_his0[lat_str]
            lat_vec = lat_var[:,0]
            x = zfun.find_nearest_ind(lon_vec,lonchoice)
            y = zfun.find_nearest_ind(lat_vec,latchoice)
            z = -1
        
        # variable to reproduce
        ot[t] = ds_avg['ocean_time'][0]
        var_accel_old[t] = ds_dia[accel_str][0,z,y,x]    
        dvar = ds_his1[varchoice][0,z,y,x]-ds_his0[varchoice][0,z,y,x] #du

        delta0 = get_delta(oceanfnh0,pmpn,z,y,x)
        delta1 = get_delta(oceanfnh1,pmpn,z,y,x)
        ddelta = delta1-delta0

        #load in averaged variables
        var = ds_avg[varchoice][0,z,y,x] #<u>
        vardelta = ds_avg[trans_str][0,z,y,x] #<udelta>
        delta = get_delta(oceanfna,pmpn,z,y,x)

        Ddvardt_oD[t] = dvar/dt
        vardDdt_oD[t] = var*ddelta/dt /delta
        var_accel_new[t] = Ddvardt_oD[t] + vardDdt_oD[t]
        
        t = t+1
        ds_his0.close()
        ds_his1.close()
        ds_dia.close()
        ds_avg.close()

# Save results for plotting
D = dict()
D_list = ['var_accel_old','Ddvardt_oD','vardDdt_oD',
          'var_accel_new','ot','varchoice',
          'lonchoice','latchoice','start_day']
for vn in D_list:
    D[vn] = locals()[vn]

out_dir = '/pmr4/eab32/etools/'
out_fn = out_dir + 'dvdt_series.p'
pickle.dump(D, open(out_fn, 'wb'))

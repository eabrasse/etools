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

dir = '/pmr4/eab32/LiveOcean_ROMS/output/aestus1_base_ae1/f2017.03.10/'

fn_his0 = dir + 'ocean_his_0001.nc'
ds_his0 = nc.Dataset(fn_his0)

fn_his1 = dir + 'ocean_his_0002.nc'
ds_his1 = nc.Dataset(fn_his1)

fn_avg = dir + 'ocean_avg_0001.nc'
ds_avg = nc.Dataset(fn_avg)

fn_dia = dir + 'ocean_dia_0001.nc'
ds_dia = nc.Dataset(fn_dia)

#first load in the variable I want to reproduce
u_accel = ds_dia['u_accel'][0,:]

#next, try to calculate it using u_accel = <delta*du/dt> / <delta> + <u*ddelta/dt> / <delta>
# need u, delta, udelta and dt, du, ddelta

dt = ds_his1['ocean_time'][:]-ds_his0['ocean_time'][:] #dt
du = ds_his1['u'][0,:]-ds_his0['u'][0,:] #du

NZ,NY,NX = du.shape

#delta is a small pain to calculate :/
on000 = 1/ds_his0['pn'][:]
on00 = 0.5*(on000[:,:-1]+on000[:,1:])
on0 = np.tile(on00,(NZ,1,1))
zeta0 = ds_his0['zeta'][:]
h = ds_his0['h'][:]
S0 = zrfun.get_basic_info(fn_his0,only_S=True)
zw0 = zrfun.get_z(h, zeta0, S0, only_w=True)
Hz00 = zw0[1:,:,:]-zw0[:-1,:,:]
Hz0 = 0.5*(Hz00[:,:,:-1]+Hz00[:,:,1:])
delta0 = Hz0*on0 #delta0

on100 = 1/ds_his1['pn'][:]
on10 = 0.5*(on100[:,:-1]+on100[:,1:])
on1 = np.tile(on10,(NZ,1,1))
zeta1 = ds_his1['zeta'][:]
S1 = zrfun.get_basic_info(fn_his1,only_S=True)
zw1 = zrfun.get_z(h, zeta1, S1, only_w=True)
Hz10 = zw1[1:,:,:]-zw1[:-1,:,:]
Hz1 = 0.5*(Hz10[:,:,:-1]+Hz10[:,:,1:])
delta1 = Hz1*on1 #delta1

ddelta = delta1-delta0

#load in averaged variables
u = ds_avg['u'][0,:] #<u>
udelta = ds_avg['Huon'][0,:] #<udelta>
delta = udelta/u # <udelta>/<delta> is this a silly way to get delta? see if it works


Ddudt_oD = du/dt
udDdt_oD = u*ddelta/dt /delta
u_accel_new = Ddudt_oD + udDdt_oD

#try to compare the different terms
print('u_accel = \n')
for i in range(5):
    print('%e \n' % u_accel[20,50,i])
    
print('u_accel_new = \n')
for i in range(5):
    print('%e \n' % u_accel_new[20,50,i])
    
print('<delta*du/dt> / <delta> = \n')
for i in range(5):
    print('%e \n' % Ddudt_oD[20,50,i])
    
print('<u*ddelta/dt> / <delta> = \n')
for i in range(5):
    print('%e \n' % udDdt_oD[20,50,i])

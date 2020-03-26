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

T0 = 60
T1 = 63
f_list=f_list[60:63]
ndays = len(f_list)
nt = ndays*24

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
        u0 = dsa['u'][:]
    
        if tt==0:
            NT, NZ, NY, NX = np.shape(u0)
            lonu = dsa['lon_u'][:]
            l0 = zfun.find_nearest_ind(lonu[0,:],0.3)
            
            #calculate A
            maskr = dsa['mask_rho']
            masku = dsa['mask_u']
            mu = np.tile(masku[:,l0],(NZ,1))
            h = dsa['h'][:,l0]
            h[maskr[:,l0]==0]=np.nan
            pn = dsa['pn'][:]
            dy = 1/pn[:,l0]
            dA = h*dy
            A = np.nansum(dA)
            H = np.nanmean(h)
            Cd = dsa['rdrg2'][:]
            
            #initialize matrix
            u_mat = np.zeros((nt, NZ, NY))

        u = u0[0,:,:,l0].squeeze()
        u[mu==0]=np.nan
        #dump data slices into filterable matrices
        u_mat[tt,:,:] = u

    
        dsa.close()
        tt = tt+1
    
#calculate Frf, M
Qr = 1500
ubar = Qr/A
g = 9.81
beta = 7.7e-4
s_ocean = 30.5
c = np.sqrt(g*beta*s_ocean*H)

Frf = ubar/c

umean = np.nanmean(u_mat)
uT = np.nanmean(np.abs(u_mat-umean))
omega = 2*np.pi/(12*60*60)
No = c/H

Msquared = (Cd*(uT**2))/(omega*No*(H**2))
M = np.sqrt(Msquared)

print('Frf = %.4f' % Frf)
print('M = %.4f' % M)

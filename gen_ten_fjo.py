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
import pickle
import numpy as np
import netCDF4 as nc

#LOAD IN VARIABLES
Ldir = Lfun.Lstart('aestus1', 'base')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'

Ldirroms = '/pmr4/eab32/LiveOcean_ROMS/'

dir0 = Ldirroms + 'output/' + Ldir['gtagex'] + '/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']
first_f = [i for i, s in enumerate(f_list) if '2017.02.14' in s]
f_list = f_list[first_f[0]:]

ocean_grid = dir0+f_list[0]+'/ocean_his_0001.nc'

trackfn = '/pmr4/eab32/LiveOcean_output/tracks/ae1_ndiv12_3d_turb/release_2017.02.14.nc'

list_choice = 'hadv','vadv','cor','prsgrd','accel','vvisc'

G, S, T = zrfun.get_basic_info(ocean_grid)
lonr = G['lon_rho'][:]
lonu = G['lon_u'][:]
lonv = G['lon_v'][:]
latr = G['lat_rho'][:]
latu = G['lat_u'][:]
latv = G['lat_v'][:]
maskr = G['mask_rho'][:]
masku = G['mask_u'][:]
maskv = G['mask_v'][:]
h = G['h'][:]


# GET PARTICLE PATHS
dst = nc.Dataset(trackfn)
tu0 = dst['u'][:]
tv0 = dst['v'][:]
tsalt0 = dst['salt'][:]
tlon0 = dst['lon'][:]
tlat0 = dst['lat'][:]
tz0 = dst['z'][:]
tot = dst['ot'][:]
tot = tot - tot[0]
tot = tot / 3600 / 24

NT,NP = tsalt0.shape
dia_list = []
avg_list = []
Ndays = int((NT-1)/24)
f_list = f_list[:Ndays+1]
for day in range(Ndays):
    f_dir = f_list[day]
    fn_list = os.listdir(dir0 + f_dir)
    fn_list.sort()
    #di0 = day*24
    #di1 = (day+1)*24
    d_list = [dir0+f_dir+'/'+x for x in fn_list if x[:9]=='ocean_dia']
    for d in d_list:
        dia_list.append(d)
    a_list = [dir0+f_dir+'/'+x for x in fn_list if x[:9]=='ocean_avg']
    for a in a_list:
        avg_list.append(a)
dia_list.append(dir0+f_list[-1]+'/ocean_dia_0001.nc')
avg_list.append(dir0+f_list[-1]+'/ocean_avg_0001.nc')

tma = np.zeros(NP,dtype=bool)
for i in range(NP):
    tma[i] = np.any(tsalt0[:,i]<20.0)

tsalt1 = tsalt0[:,tma]
tlon1 = tlon0[:,tma]
tlat1 = tlat0[:,tma]
tz1 = tz0[:,tma]
tu1 = tu0[:,tma]
tv1 = tv0[:,tma]

reflon=0.0
reflat=45.0

#CHOOSE INDIVIDUAL PARTICLE
particlelist = range(898,909)
plist_len = len(particlelist)
N = 80
for particlechoice in particlelist:
    print('Particle %i' % particlechoice)
    tlon = tlon1[:,particlechoice]
    tlat = tlat1[:,particlechoice]
    tz = tz1[:,particlechoice]
    tu = tu1[:,particlechoice]
    tv = tv1[:,particlechoice]
    
    if particlechoice==particlelist[0]:
        pdict = dict()
    
    dist_from_mouth = np.zeros(NT)
    for tt in range(NT):
        x,y =zfun.ll2xy(tlon[tt],tlat[tt],reflon,reflat)
        dist_from_mouth[tt] = np.sqrt(x**2 + y**2)
        dist_from_mouth[tt] = np.sign(tlon[tt]) * np.abs(dist_from_mouth[tt])
    
    #initialize dicts
    Du = dict()
    Dv = dict()

    #to build a list of momentum terms by position, go through each time step...
    print('reading in momentum budgets along track...')
    for i in range(NT):
        #open the relevant model output...
        ds_dia = nc.Dataset(dia_list[i])
        ds_avg = nc.Dataset(avg_list[i])
        
        zeta = ds_avg['zeta'][:]
        zr = zrfun.get_z(h,zeta,S,only_rho=True)
        
        #read in all of the momentum terms of interest...
        for key in sorted(list_choice):
            
            if i==0:
                Du[key] = np.zeros((NT))
                Dv[key] = np.zeros((NT))
                
            u_key = 'u_' + key
            v_key = 'v_' + key
            
            #for u and v at the particle's position for that time step
            indx_u = zfun.find_nearest_ind(lonu[0,:],tlon[i])
            indy_u = zfun.find_nearest_ind(latu[:,0],tlat[i])
            indz_u = zfun.find_nearest_ind(zr[:,indy_u,indx_u+1],tz[i]) #indx+1 because zr is on rho grid
            Du[key][i] = ds_dia[u_key][0,indz_u,indy_u,indx_u]
            
            indx_v = zfun.find_nearest_ind(lonv[0,:],tlon[i])
            indy_v = zfun.find_nearest_ind(latv[:,0],tlat[i])
            indz_v = zfun.find_nearest_ind(zr[:,indy_v+1,indx_v],tz[i]) #indy+1 because zr is on rho grid
            Dv[key][i] = ds_dia[v_key][0,indz_v,indy_v,indx_v]
        
        if np.mod(i,100)==0:
            print('time step %i' % i)
        #(don't forget to close the files!)
        ds_dia.close()
        ds_avg.close()

    theta = np.arctan(tv/tu)

    Da = dict()
    Dc = dict()

    print('rotating coordinates')
    for key in list_choice:
        # #rotate for the side plots
        alongtrack = np.zeros(NT)
        crosstrack = np.zeros(NT)
        for t in range(NT):
            alongtrack[t] = Du[key][t]*np.cos(theta[t]) + Dv[key][t]*np.sin(theta[t])
            crosstrack[t] = -Du[key][t]*np.sin(theta[t]) + Dv[key][t]*np.cos(theta[t])
        Da[key] = alongtrack
        Dc[key] = crosstrack
    
    # print('smoothing')
    # #should do some sort of smoothing
    # tlon_sm = zfun.filt_hanning(tlon,n=N)
    # tlon_sm = tlon_sm[~np.isnan(tlon_sm)]
    # tlat_sm = zfun.filt_hanning(tlat,n=N)
    # tlat_sm = tlat_sm[~np.isnan(tlat_sm)]
    # tz_sm = zfun.filt_hanning(tz,n=N)
    # tz_sm = tz_sm[~np.isnan(tz_sm)]
    # tot_sm = zfun.filt_hanning(tot,n=N)
    # tot_sm = tot_sm[~np.isnan(tot_sm)]
    # tu_sm = zfun.filt_hanning(tu,n=N)
    # tu_sm = tu_sm[~np.isnan(tu_sm)]
    # tv_sm = zfun.filt_hanning(tv,n=N)
    # tv_sm = tv_sm[~np.isnan(tv_sm)]
    # for key in Da.keys:
    #
    
    pdict[particlechoice] = {'tlon':tlon,'tlat':tlat,'tz':tz,
                            'tu':tu,'tv':tv,'tot':tot,'dist_from_mouth':dist_from_mouth,
                            'Da':Da,'Dc':Dc,}

dst.close()
outfn = 'pdict.p'
pickle.dump( pdict, open( outfn, 'wb'))

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

print('initializing...')

#LOAD IN VARIABLES
Ldir = Lfun.Lstart('aestus2', 'base')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'

Ldirroms = '/pmr4/eab32/LiveOcean_ROMS/'

dir0 = Ldirroms + 'output/' + Ldir['gtagex'] + '/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']
first_f = [i for i, s in enumerate(f_list) if '2017.02.15' in s]
f_list = f_list[first_f[0]:]

ocean_grid = dir0+f_list[0]+'/ocean_his_0001.nc'

trackfn = '/pmr4/eab32/LiveOcean_output/tracks/eliz_ae1_ndiv12_3d_turb/release_2017.02.15.nc'

dia_keys = 'hadv','vadv','cor','prsgrd','accel','vvisc'

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
dt = tot[2] - tot[1]
tot = tot - tot[0]
tot = tot / 3600 / 24

NT,NP = tsalt0.shape

tma = np.zeros(NP,dtype=bool)
for i in range(NP):
    tma[i] = np.any(tsalt0[:,i]<20.0)

tsalt1 = tsalt0[:,tma]
tlon1 = tlon0[:,tma]
tlat1 = tlat0[:,tma]
tz1 = tz0[:,tma]
tu1 = tu0[:,tma]
tv1 = tv0[:,tma]

Ndays = int((NT-1)/24)
reflon=0.0
reflat=45.0

#CHOOSE INDIVIDUAL PARTICLE
particlelist = range(np.sum(tma))
#initialize pdict
pdict = dict()
for particle in particlelist:
    pdict[particle] = dict()
    
    #loat in particle position and velocity
    pdict[particle]['tlon'] = tlon1[:,particle]
    pdict[particle]['tlat'] = tlat1[:,particle]
    pdict[particle]['tz'] = tz1[:,particle]
    pdict[particle]['tu'] = tu1[:,particle]
    pdict[particle]['tv'] = tv1[:,particle]
    
    #create new particle variables
    pdict[particle]['dist_from_mouth'] = np.zeros((NT))
    pdict[particle]['hab'] = np.zeros((NT))
    
    #create empty dicts to fill up with diagnostic data
    pdict[particle]['Da'] = dict()
    pdict[particle]['Dc'] = dict()
    
    pdict[particle]['Da']['unsteady'] = np.zeros((NT))
    pdict[particle]['Dc']['unsteady'] = np.zeros((NT))
    for key in dia_keys:
        pdict[particle]['Da'][key] = np.zeros((NT))
        pdict[particle]['Dc'][key] = np.zeros((NT))

# Now fill in dicts, day by day
tt = 0
for day in range(Ndays):
    print('running day %i' % day)
    f_dir = f_list[day]
    fn_list = os.listdir(dir0 + f_dir)
    fn_list.sort()
    
    d_list = [dir0+f_dir+'/'+x for x in fn_list if x[:9]=='ocean_dia']
    a_list = [dir0+f_dir+'/'+x for x in fn_list if x[:9]=='ocean_avg']
    h_list = [dir0+f_dir+'/'+x for x in fn_list if x[:9]=='ocean_his']
    
    hrs = len(d_list)
    for hr in range(hrs):
        ds_dia = nc.Dataset(d_list[hr])
        ds_avg = nc.Dataset(a_list[hr])
        ds_h0 = nc.Dataset(h_list[hr])
        ds_h1 = nc.Dataset(h_list[hr+1])
        
        zeta = ds_avg['zeta'][:]
        zr = zrfun.get_z(h,zeta,S,only_rho=True)

        for particle in particlelist: 
            p = pdict[particle]
            
            #build dist_from_mouth      
            x,y =zfun.ll2xy(p['tlon'][tt],p['tlat'][tt],reflon,reflat)
            p['dist_from_mouth'][tt] = np.sqrt(x**2 + y**2)
            p['dist_from_mouth'][tt] = np.sign(p['tlon'][tt]) * np.abs(p['dist_from_mouth'][tt])
            
            #build hab
            indx_r = zfun.find_nearest_ind(lonr[0,:],p['tlon'][tt])
            indy_r = zfun.find_nearest_ind(latr[:,0],p['tlat'][tt])
            p['hab'][tt] = h[indy_r,indx_r] + p['tz'][tt] #h is (+), tz is (-), and |h|>|tz| so adding them should produce a (+) difference
            
            theta = np.arctan(p['tv'][tt]/p['tu'][tt])
            
            #calculate du/dt
            indx_u = zfun.find_nearest_ind(lonu[0,:],p['tlon'][tt])
            indy_u = zfun.find_nearest_ind(latu[:,0],p['tlat'][tt])
            indz_u = zfun.find_nearest_ind(zr[:,indy_u,indx_u+1],p['tz'][tt]) #indx+1 because zr is on rho grid
            h0_u = ds_h0['u'][0,indz_u,indy_u,indx_u]
            h1_u = ds_h1['u'][0,indz_u,indy_u,indx_u]
            dudt = (h1_u-h0_u)/dt
        
            indx_v = zfun.find_nearest_ind(lonv[0,:],p['tlon'][tt])
            indy_v = zfun.find_nearest_ind(latv[:,0],p['tlat'][tt])
            indz_v = zfun.find_nearest_ind(zr[:,indy_v+1,indx_v],p['tz'][tt]) #indy+1 because zr is on rho grid
            h0_v = ds_h0['v'][0,indz_v,indy_v,indx_v]
            h1_v = ds_h1['v'][0,indz_v,indy_v,indx_v]
            dvdt = (h1_v-h0_v)/dt
            
            p['Da']['unsteady'][tt] = dudt*np.cos(theta) + dvdt*np.sin(theta)
            p['Dc']['unsteady'][tt] = -dudt*np.sin(theta) + dvdt*np.cos(theta)
            
            #build up diagnostic terms
            for key in dia_keys:
                u_key = 'u_' + key
                v_key = 'v_' + key
                
                indx_u = zfun.find_nearest_ind(lonu[0,:],p['tlon'][tt])
                indy_u = zfun.find_nearest_ind(latu[:,0],p['tlat'][tt])
                indz_u = zfun.find_nearest_ind(zr[:,indy_u,indx_u+1],p['tz'][tt]) #indx+1 because zr is on rho grid
                dia_u = ds_dia[u_key][0,indz_u,indy_u,indx_u]
            
                indx_v = zfun.find_nearest_ind(lonv[0,:],p['tlon'][tt])
                indy_v = zfun.find_nearest_ind(latv[:,0],p['tlat'][tt])
                indz_v = zfun.find_nearest_ind(zr[:,indy_v+1,indx_v],p['tz'][tt]) #indy+1 because zr is on rho grid
                dia_v = ds_dia[v_key][0,indz_v,indy_v,indx_v]
                
                p['Da'][key][tt] = dia_u*np.cos(theta) + dia_v*np.sin(theta)
                p['Dc'][key][tt] = -dia_u*np.sin(theta) + dia_v*np.cos(theta)
        
        tt = tt + 1
        ds_dia.close()
        ds_avg.close()
        ds_h0.close()
        ds_h1.close()
print('NT = %i' % NT)
print('tt = %i' % tt)
print('Saving!')
dst.close()
outfn = 'pdict0_aestus2_unstrat.p'
pickle.dump( pdict, open( outfn, 'wb'))

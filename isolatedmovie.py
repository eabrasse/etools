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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cmocean as cmo

plp = os.path.abspath('/pmr4/eab32/LiveOcean/plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

import netCDF4 as nc4
import numpy as np

#get ocean background
dir0 = '/pmr4/eab32/LiveOcean_ROMS/output/aestus1_base_ae1/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

f_list=f_list[45:]

trackfn = '/pmr4/eab32/LiveOcean_output/tracks/ae1_ndiv12_3d_turb/release_2017.02.14.nc'
dst = nc4.Dataset(trackfn)

#get particle positions
salt = dst['salt'][:]
ot = dst['ot'][:]
lon0 = dst['lon'][:]
lat0 = dst['lat'][:]
z0 = dst['z'][:]

time = (ot-ot[0])/3600
time = time/24 + 45
ndays=len(f_list)
cmap = cmo.cm.haline
nt = ndays * 24
# develop a filter for only particles that start at PSU =>30 and end at PSU <30
start_as_ocean = salt[0,:] > 29.9
felloff = np.zeros(len(start_as_ocean),dtype=bool)
end_as_channel = np.zeros(len(start_as_ocean),dtype=bool)
end_as_estuary = np.zeros(len(start_as_ocean),dtype=bool)
for i in np.arange(len(end_as_channel)):
    end_as_estuary[i] = np.any(salt[:nt,i]<15.0)
    end_as_channel[i] = np.any(salt[:nt,i]<20.0) #stricter here

# color code: start_as_ocean * end_as_estuary * end_as_channel = 'goldenrod', #0
#            start_as_ocean * end_as_estuary * ~end_as_channel = 'seagreen', #1
#           start_as_ocean * ~end_as_estuary * ~end_as_channel = 'royalblue', #2
tma = start_as_ocean * end_as_estuary# * end_as_channel
colorname = 'goldenrod'

#create slimmed down variables: initial pos only, and mask out particles that don't meet criteria
lon = lon0[0,tma]
lat = lat0[0,tma]
z = z0[0,tma]

lpad = 0.1
zpad = 10

for t in range(ndays):
    
    fn = f_list[t] + '/ocean_his_0001.nc'
    
    oceanfn = dir0 + fn
    
    dso = nc4.Dataset(oceanfn)
    
    if t==0:
        #get grid info from ocean
        lonr = dso['lon_rho'][:]
        latr = dso['lat_rho'][:]
        h = dso['h'][:]
        zeta = dso['zeta'][:]
        maskr = dso['mask_rho'][:]
        S = zrfun.get_basic_info(oceanfn,only_S=True)
        zr = zrfun.get_z(h,zeta,S,only_rho=True)
        NZ,NY,NX = np.shape(zr)
        midlat = zfun.find_nearest_ind(latr[:,0],45.0)
        
    bgsalt = dso['salt'][:]
    surfsalt = bgsalt[0,-1,:,:].squeeze()
    crosssalt = bgsalt[0,:,midlat,:].squeeze()
    nt = t * 24
    
    lon1 = lon0[nt,tma]
    lat1 = lat0[nt,tma]
    z1 = z0[nt,tma]
    
    fig = plt.figure(figsize=(12,4))
    
    ax1 = fig.add_subplot(1,3,2)
    a_ll = [lonr.min(), 1.2, latr.min(), latr.max()]
    bgs = np.ma.masked_where(maskr==0, surfsalt)
    ax1.pcolormesh(lonr, latr, bgs, vmin=0, vmax=30,
        cmap=cmap, alpha=1.0)
    ax1.axis(a_ll)
    ax1.set_aspect(1)
    ax1.set_xlabel('Longitude',fontweight='bold')
    ax1.set_ylabel('Latitude',fontweight='bold')
    #ax1.plot(lon,lat,'og',alpha=0.3)
    ax1.plot(lon1,lat1,marker='o',color='black',markerfacecolor='magenta',alpha=1.0,linestyle='None')
    #ax1.plot(lon0[:nt:24,tma],lat0[:nt:24,tma],'-k',alpha=0.2,linewidth=0.5)
    ax1.set_title('Plan view',fontweight='bold')

    ax2 = fig.add_subplot(1,3,3)
    a_lz = [lonr.min()+lpad, 1.2, -60, 5]
    p = ax2.pcolormesh(lonr[:NZ,:], zr[:,midlat,:].squeeze(), crosssalt, vmin=0, vmax=30,
        cmap=cmap, alpha=1.0)
    plt.colorbar(p)
    ax2.axis(a_lz)
    pfun.dar(ax2)
    ax2.set_aspect(1/80)
    ax2.set_xlabel('Longitude',fontweight='bold')
    ax2.set_ylabel('z (m)',fontweight='bold')
    #plot particles
    #ax2.plot(lon,z,'og',alpha=0.3)
    ax2.plot(lon1,z1,marker='o',color='black',markerfacecolor='magenta',alpha=1.0,linestyle='None')
    #ax2.plot(lon0[:nt:24,tma],z0[:nt:24,tma],'-k',alpha=0.2,linewidth=0.5)
    ax2.set_title('Cross shelf section',fontweight='bold')

    ax4 = fig.add_subplot(1,3,1)
    ax4.set_xlabel('Time (day)',fontweight='bold')
    ax4.set_ylabel('Salt (PSU)',fontweight='bold')
    ax4.plot(time[:nt:24],salt[:nt:24,tma],alpha=0.5,color=colorname)
    ax4.set_ylim([0,34])
    ax4.set_xlim([43,122])

    fig.subplots_adjust(wspace=0.3,hspace=0.2)
    dt = t + 45
    dts = str(dt).zfill(3)
    titlestr = 'Day %s' % dts
    outname = 'movie/pic_%s.png' % dts
    plt.suptitle(titlestr, fontsize = 14, fontweight = "bold")
    plt.savefig(outname)
    
    plt.close()
    dso.close()

dst.close()

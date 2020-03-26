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
ep = os.path.abspath('/pmr4/eab32/etools')
if ep not in sys.path:
    sys.path.append(ep)
import efun
plp = os.path.abspath('/pmr4/eab32/LiveOcean/plotting')
if plp not in sys.path:
    sys.path.append(plp)
import pfun

import netCDF4 as nc
import numpy as np

outdir = 'bottompressure'
Lfun.make_dir(outdir,clean=True)

#get ocean background
dir0 = '/pmr4/eab32/LiveOcean_ROMS/output/aestus1_base_ae1/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

#f_list=f_list[0::10]
f_list = f_list[:32]

ndays=len(f_list)
cmap1 = 'PuOr'
cmap2 = 'bwr'

fn_list = os.listdir(dir0 + f_list[0])
fnd = [x for x in fn_list if x[-11:-8]=='dia']
fna = [x for x in fn_list if x[-11:-8]=='avg']

oceanfnd = dir0 + f_list[0] + '/' + fnd[0]
oceanfna = dir0 + f_list[0] + '/' + fna[0]

dsd = nc.Dataset(oceanfnd)
dsa = nc.Dataset(oceanfna)

#load in remaining variable
u_prsgrd = dsd['u_prsgrd'][0,0,:,:]
v = dsa['v'][0,0,:,:]
salt = dsa['salt'][0,-1,:,:]

#get grid info from ocean
G, S, T = zrfun.get_basic_info(oceanfnd)
lonu = G['lon_u'][:]
lonv = G['lon_v'][:]
lonr = G['lon_rho'][:]
latu = G['lat_u'][:]
latv = G['lat_v'][:]
latr = G['lat_rho'][:]
masku = G['mask_u'][:]
maskr = G['mask_rho'][:]
maskv = G['mask_v'][:]

a_ll = [-1.0, 0.02, latr.min(), latr.max()]

bb = efun.box_inds(-0.65,-0.1,44.2,45.8,oceanfna)
bu = np.array([[bb['lonu0'],bb['lonu1']],[bb['latu0'],bb['latu1']]])
bv = np.array([[bb['lonv0'],bb['lonv1']],[bb['latv0'],bb['latv1']]])

upmin = u_prsgrd[bu[0,0]:bu[0,1],bu[1,0]:bu[1,1]].min()
upmax = u_prsgrd[bu[0,0]:bu[0,1],bu[1,0]:bu[1,1]].max()
upm = np.max([np.abs(upmin),np.abs(upmax)])*10
vmin = v[bv[0,0]:bv[0,1],bv[1,0]:bv[1,1]].min()
vmax = v[bv[0,0]:bv[0,1],bv[1,0]:bv[1,1]].max()
vm = np.max([np.abs(vmin),np.abs(vmax)])*50

for t in range(ndays):
    #identify file names for that day
    fn_list = os.listdir(dir0 + f_list[t])
    
    fnd = [x for x in fn_list if x[-11:-8]=='dia']
    fna = [x for x in fn_list if x[-11:-8]=='avg']
    
    #load in data for each hour and find the average
    for h in np.arange(0,24):
    
        oceanfnd = dir0 + f_list[t] + '/' + fnd[h]
        oceanfna = dir0 + f_list[t] + '/' + fna[h]
    
        dsd = nc.Dataset(oceanfnd)
        dsa = nc.Dataset(oceanfna)

        u_prsgrd = dsd['u_prsgrd'][0,0,:,:]
        v = dsa['v'][0,0,:,:]
        salt = dsa['salt'][0,-1,:,:]

        if h==0:
            upl = np.ma.masked_where(masku==0, u_prsgrd)
            vl = np.ma.masked_where(maskv==0, v)
            sl = np.ma.masked_where(maskr==0, salt)

        else:
            upl = upl + np.ma.masked_where(masku==0, u_prsgrd)
            vl = vl + np.ma.masked_where(maskv==0, v)
            sl = sl + np.ma.masked_where(maskr==0, salt)
        
        dsd.close()
        dsa.close()
    #take the average
    upl = upl/24
    vl = vl/24
    sl = sl/24
    
    #plot averaged data
    fig=plt.figure(figsize=(14,8))

    #plot surface plan view
    ax1 = fig.add_subplot(1,2,1)
    p = ax1.pcolormesh(lonu, latu, upl, cmap=cmap1,vmin=-upm, vmax = upm,alpha=1.0)
    plt.colorbar(p)
    ax1.contour(lonr,latr,sl,[27],colors='k',linewidths=1.0)
    ax1.axis(a_ll)
    ax1.set_xlabel('Longitude',fontweight='bold')
    ax1.set_ylabel('Latitude',fontweight='bold')
    ax1.set_title('X direction pressure gradient at seafloor')
    pfun.dar(ax1)

    ax2 = fig.add_subplot(1,2,2)
    p = ax2.pcolormesh(lonv, latv, vl, cmap=cmap2,vmin=-vm, vmax = vm,alpha=1.0)
    plt.colorbar(p)
    ax2.contour(lonr,latr,sl,[27],colors='k',linewidths=1.0)
    ax2.axis(a_ll)
    ax2.set_xlabel('Longitude',fontweight='bold')
    ax2.set_ylabel('Latitude',fontweight='bold')
    ax2.set_title('N-S velocity at seafloor')
    pfun.dar(ax2)


    fig.subplots_adjust(hspace=0.3)
    titlestr = '%s' % f_list[t]
    plt.suptitle(titlestr)
    tstr = str(t).zfill(3)
    outname = outdir + '/pic_%s.png' % tstr
    plt.savefig(outname)

    plt.close()

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

outdir = 'density'
Lfun.make_dir(outdir,clean=True)

#get ocean background
dir0 = '/pmr4/eab32/LiveOcean_ROMS/output/aestus1_base_ae1/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

#f_list=f_list[0::10]
#f_list = f_list[:50]

ndays=len(f_list)
cmap = 'rainbow'

fn_list = os.listdir(dir0 + f_list[0])
fna = [x for x in fn_list if x[-11:-8]=='avg']

oceanfna = dir0 + f_list[0] + '/' + fna[0]

dsa = nc.Dataset(oceanfna)

#load in remaining variable
rho = dsa['rho'][0,:,:,:]
NZ, NY, NX = np.shape(rho)


#get grid info from ocean
G, S, T = zrfun.get_basic_info(oceanfna)
lonr = G['lon_rho'][:]
latr = G['lat_rho'][:]
maskr0 = G['mask_rho'][:]
H = G['h'][:]

maskr = np.tile(maskr0,(NZ,1,1))

# create cross sectional axes
latvec = latr[:,0]
l_list = [45.6, 45.3, 45.0, 44.4]
llen = len(l_list)
l = np.zeros(len(l_list),dtype='int')
lonvec = np.zeros([len(l_list),np.shape(lonr)[1]])
for lat in range(len(l_list)):
      l0 = zfun.find_nearest_ind(latvec,l_list[lat])
      l[lat] = l0
      lonvec[lat,:] = lonr[l0,:]

a_ll = [-1.0, 0.02, latr.min(), latr.max()]

vmin = 22.7
vmax = 23.9

for t in range(ndays):
    #identify file names for that day
    fn_list = os.listdir(dir0 + f_list[t])
    
    fna = [x for x in fn_list if x[-11:-8]=='avg']
    
    #load in data for each hour and find the average
    for h in np.arange(0,24):
    
        oceanfna = dir0 + f_list[t] + '/' + fna[h]
    
        dsa = nc.Dataset(oceanfna)

        rho = dsa['rho'][0,:,:,:]
        zeta = dsa['zeta'][0,:,:]

        if h==0:
            rl = np.ma.masked_where(maskr==0, rho)
            zl = np.ma.masked_where(maskr0==0, zeta)
        else:
            rl = rl + np.ma.masked_where(maskr==0, rho)
            zl = zl + np.ma.masked_where(maskr0==0, zeta)
        dsa.close()
    #take the average
    rl = rl/24
    zl = zl/24

    zl = np.reshape(zl,(1,NY,NX))

    # fill out axis definitions with time varying z coordinate
    S = zrfun.get_basic_info(oceanfna,only_S=True)
    zr, zw = zrfun.get_z(H,zl,S)
    zvec = np.zeros([np.shape(zr)[0],len(l_list),np.shape(zr)[2]])
    for lat in range(len(l_list)):
        zvec[:,lat,:] = zr[:,l[lat],:]
    a_lz = [-1.0, 0.5, zvec.min()-5, zvec.max()+5]

    #plot averaged data
    #begin PLOT
    fig=plt.figure(figsize=(10,7))

    #plot surface plan view
    ax1 = plt.subplot2grid((llen,2),(0,0),colspan=1,rowspan=llen)
    rl = np.ma.masked_where(maskr0==0, rho[-1,:,:])
    p = ax1.pcolormesh(lonr, latr, rl, cmap='rainbow',vmin=vmin, vmax = vmax,alpha=0.7)
    plt.colorbar(p)

    ax1.axis(a_ll)
    ax1.set_xlabel('Longitude',fontweight='bold')
    ax1.set_ylabel('Latitude',fontweight='bold')
    ax1.text(0.8,0.9,'Surface',transform = ax1.transAxes,horizontalalignment = 'center',verticalalignment='center')

    for lat in range(len(l)):
        if l_list[lat]==45.0:
            bb = efun.box_inds(-0.65,1.5,44.0,46.0,oceanfna)
        else:
            bb = efun.box_inds(-0.65,0.0,44.0,46.0,oceanfna)
        ax1.plot(lonr[l[lat],:bb['lonr1']],latr[l[lat],:bb['lonr1']],'k--',linewidth=3.0)
    
        #plot density crossection
        ax = plt.subplot2grid((llen,2),(lat,1),colspan=1,rowspan=1)
        rz = rho[:,l[lat],:]
        xx,zz = np.meshgrid(lonvec[lat,:],zvec[:,lat,0])

        Nlist = np.linspace(vmin,vmax,15)
        ax.pcolormesh(xx, zvec[:,lat,:], rz, vmin=vmin, vmax = vmax,cmap='rainbow',alpha=0.7)
        ax.contour(xx,zvec[:,lat,:],rz,Nlist,colors='k',linewidths=0.5)
        ax.axis(a_lz)
        ax.set_xlabel('Longitude',fontweight='bold')
        ax.set_ylabel('Z (m)',fontweight='bold')
        axtext = 'Lat = %.1f' % l_list[lat]
    
        ax.text(0.7,0.1,axtext,transform = ax.transAxes,horizontalalignment = 'center',verticalalignment='center')
    
        if lat<llen-1:
            ax.set_xlabel('')


    fig.subplots_adjust(wspace=0.3,hspace=0.5)
    titlestr = 'Density contours %s' % f_list[t]
    plt.suptitle(titlestr)
    tstr = str(t).zfill(3)
    outname = outdir + '/pic_%s.png' % tstr
    plt.savefig(outname)

    plt.close()

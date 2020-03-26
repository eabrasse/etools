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
ep = os.path.abspath('/pmr4/eab32/etools')
if ep not in sys.path:
    sys.path.append(ep)
import efun
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

outdir = 'seafloor_prsgrd'
Lfun.make_dir(outdir, clean=True)

# f_list=f_list[45:]
ndays = len(f_list)
nt = ndays*24

lpad = 0.1
zpad = 10
ss=4
dia_scale = 100*3600*24 #m/s^2 to (cm/s)/day

for t in range(ndays):
    print('day %i' % t)
    
    fnd = f_list[t] + '/ocean_dia_0001.nc'
    fna = f_list[t] + '/ocean_avg_0001.nc'
    
    oceanfna = dir0 + fna
    oceanfnd = dir0 + fnd
    
    dsa = nc4.Dataset(oceanfna)
    dsd = nc4.Dataset(oceanfnd)
    
    if t==0:
        #get grid info from ocean
        lonr = dsa['lon_rho'][:]
        latr = dsa['lat_rho'][:]
        lonr_s = lonr[1:-1,1:-1]
        latr_s = latr[1:-1,1:-1]
        maskr = dsa['mask_rho'][:]
        b0 = efun.box_inds(-0.65,0.0,44.0,46.0,oceanfna)
        # lonu = dso['lon_u'][:]
        # latu = dso['lat_u'][:]
        # lonv = dso['lon_v'][:]
        # latv = dso['lat_v'][:]
        
    #plot background density
    rho = dsa['rho'][:]
    seafloor_rho = rho[0,0,:,:].squeeze()
    
    fig = plt.figure(figsize=(6,6))
    
    ax = plt.gca()
    a_ll = [lonr.min(), 1.2, latr.min(), latr.max()]
    rs = np.ma.masked_where(maskr==0, seafloor_rho)
    ax.pcolormesh(lonr, latr, rs, vmin=-30, vmax=30,
        cmap='RdBu', alpha=.5)
    ax.axis(a_ll)
    pfun.dar(ax)
    ax.set_xlabel('Longitude',fontweight='bold')
    ax.set_ylabel('Latitude',fontweight='bold')


    # add prsgrd quiver
    ax.set_title('Prsgrd + cor vectors',fontweight='bold')
    u_prsgrd = dsd['u_prsgrd'][:]
    v_prsgrd = dsd['v_prsgrd'][:]
    u_cor = dsd['u_cor'][:]
    v_cor = dsd['v_cor'][:]
    u = dia_scale*0.5*(u_prsgrd[0,0,1:-1,1:] + u_prsgrd[0,0,1:-1,:-1]+u_cor[0,0,1:-1,1:] + u_cor[0,0,1:-1,:-1])
    v = dia_scale*0.5*(v_prsgrd[0,0,1:,1:-1] + v_prsgrd[0,0,:-1,1:-1]+v_cor[0,0,1:,1:-1] + v_cor[0,0,:-1,1:-1])
    
    ax.quiver(lonr_s[b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
                latr_s[b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
                u[b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
                v[b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
                color='black',units='y',scale=500,scale_units='y',width=0.01)
    #quiver legend
    ax.quiver([0.7],[44.3],[100.0],[100.0],
                color='black',units='y',scale=500,scale_units='y',width=0.01)
    ax.text(0.7,44.5,'100 cm/s per day')

    dts = str(t).zfill(3)
    titlestr = 'Day %s' % dts
    outname = outdir + '/pic_%s.png' % dts
    plt.suptitle(titlestr, fontsize = 14, fontweight = "bold")
    plt.savefig(outname)
    
    plt.close()
    dsa.close()
    dsd.close()

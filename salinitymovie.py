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

f0 = 0
f_list=f_list[f0:f0+11]

outdir = 'salinity_movie'
Lfun.make_dir(outdir,clean=True)

#get particle positions
ndays=len(f_list)
cmap = 'rainbow'
nt = ndays * 24

lpad = 0.1
zpad = 10
tt=0
for t in range(ndays):
    print('Day %i' % t)
    fn_list = os.listdir(dir0+f_list[t])
    fn_list.sort()
    fna_list = [x for x in fn_list if x[6:9]=='avg']

    hrs = len(fna_list)
    for hr in range(hrs):
    
        oceanfn = dir0 + f_list[t] + '/' + fna_list[hr]
    
        dso = nc4.Dataset(oceanfn)
    
        if tt==0:
            #get grid info from ocean
            lonr = dso['lon_rho'][:]
            latr = dso['lat_rho'][:]
            maskr = dso['mask_rho'][:]

        salt = dso['salt'][:]
        surfsalt = salt[0,-1,:,:].squeeze()

    
        fig = plt.figure(figsize=(6,4))
    
        ax1 = fig.gca()
        a_ll = [lonr[0,6], 1.2, latr[6,0], latr[-6,0]]
        ss = np.ma.masked_where(maskr==0, surfsalt)
        p=ax1.pcolormesh(lonr, latr, ss, vmin=25, vmax=30,
        cmap=cmap, alpha=1.0)
        cbar = plt.colorbar(p)
        cbar.ax.set_ylabel('PSU',rotation=270)
        cbar.ax.get_yaxis().labelpad = 15
        ax1.axis(a_ll)
        pfun.dar(ax1)
        ax1.set_xlabel('Longitude',fontweight='bold')
        ax1.set_ylabel('Latitude',fontweight='bold')
        ax1.set_title('Surface salinity',fontweight='bold')

        dts = str(tt).zfill(3)
        titlestr = 'Day %i hour %i' % (f0+t,hr)
        outname = outdir + '/pic_%s.png' % dts
        plt.suptitle(titlestr, fontsize = 14, fontweight = "bold")
        plt.savefig(outname)
    
        plt.close()
        dso.close()
        tt = tt+1

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

outdir = 'seafloor_density'
Lfun.make_dir(outdir, clean=True)

# f_list=f_list[45:]
ndays = len(f_list)
nt = ndays*24

lpad = 0.1
zpad = 10

for t in range(ndays):
    print('day %i' % t)
    
    fn = f_list[t] + '/ocean_avg_0001.nc'
    
    oceanfn = dir0 + fn
    
    dso = nc4.Dataset(oceanfn)
    
    if t==0:
        #get grid info from ocean
        r0 = dso['rho'][:]
        seafloor_rho0 = r0[0,0,:,:].squeeze()
        
        lonr = dso['lon_rho'][:]
        latr = dso['lat_rho'][:]
        maskr = dso['mask_rho'][:]
        # lonu = dso['lon_u'][:]
        # latu = dso['lat_u'][:]
        # lonv = dso['lon_v'][:]
        # latv = dso['lat_v'][:]
        
    rho = dso['rho'][:]
    # u_prsgrd =
    seafloor_rho = rho[0,0,:,:].squeeze()
    seafloor_rho = seafloor_rho - seafloor_rho0
    
    fig = plt.figure(figsize=(6,6))
    
    ax = plt.gca()
    a_ll = [lonr.min(), 1.2, latr.min(), latr.max()]
    rs = np.ma.masked_where(maskr==0, seafloor_rho)
    p=ax.pcolormesh(lonr, latr, rs, vmin=-0.8, vmax=0.8,
        cmap='RdBu', alpha=1.0)
    ax.axis(a_ll)
    pfun.dar(ax)
    ax.set_xlabel('Longitude',fontweight='bold')
    ax.set_ylabel('Latitude',fontweight='bold')
    ax.set_title(r'Density anomaly ($kg/m^{3}$)',fontweight='bold')
    plt.colorbar(p)

    dts = str(t).zfill(3)
    titlestr = 'Day %s' % dts
    outname = outdir + '/pic_%s.png' % dts
    plt.suptitle(titlestr, fontsize = 14, fontweight = "bold")
    plt.savefig(outname)
    
    plt.close()
    dso.close()

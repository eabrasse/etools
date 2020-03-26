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
import pickle
import scipy.interpolate as intp


seafloor = True
if seafloor:
    data_fn = '/pmr4/eab32/etools/seafloor_data.p'
    zind = 5
    outdir = 'seafloor_movie'
    ax0_title = 'Density anomaly near seafloor'
    vmax= 0.2
    qscale = 250
    quiverleg = 50
else:
    data_fn = '/pmr4/eab32/etools/seasurface_data.p'
    zind = -5
    outdir = 'seasurface_movie'
    ax0_title = 'Density anomaly near sea surface'
    vmax= 5.0
    qscale = 250
    quiverleg = 50

#get ocean background
gridfn = '/pmr4/eab32/LiveOcean_ROMS/output/aestus1_base_ae1/f2017.01.01/ocean_his_0001.nc'

ds = nc4.Dataset(gridfn)
density0 = ds['rho'][0,zind,:,:].squeeze()
ds.close()

G = zrfun.get_basic_info(gridfn, only_G=True)
lonr = G['lon_rho'][:]
latr = G['lat_rho'][:]
maskr = G['mask_rho'][:]
maskuu = G['mask_u'][:]
maskvv = G['mask_v'][:]
lonr_s = lonr[1:-1,1:-1]
latr_s = latr[1:-1,1:-1]
b0 = efun.box_inds(-0.65,0.0,44.0,46.0,gridfn)
b1 = efun.box_inds(-0.2,1.5,44.8,45.2,gridfn)

D = pickle.load(open(data_fn,'rb'))

density = D['density_gf'][:]
upd = D['u_prsgrd_gf'][:]
vpd = D['v_prsgrd_gf'][:]
uad = D['u_ageo_gf'][:]
vad = D['v_ageo_gf'][:]

nt, ny, nx = np.shape(density)
maskr_s = maskr[1:-1,1:-1]
maskrr = np.reshape(maskr_s,(1,np.shape(maskr_s)[0],np.shape(maskr_s)[1]))
mask0 = np.tile(maskrr,(nt,1,1))

nngrid = 30
# set masked values to 0
upd[~mask0]=0
vpd[~mask0]=0
uad[~mask0]=0
vad[~mask0]=0

# mask out mouth
upd[:,b1['latu0']:b1['latu1'], b1['lonu0']:b1['lonu1']] = 0
vpd[:,b1['latv0']:b1['latv1'], b1['lonv0']:b1['lonv1']] = 0
uad[:,b1['latu0']:b1['latu1'], b1['lonu0']:b1['lonu1']] = 0
vad[:,b1['latv0']:b1['latv1'], b1['lonv0']:b1['lonv1']] = 0


Lfun.make_dir(outdir, clean=True)

ss=4
dia_scale = 100*3600*24 #m/s^2 to (cm/s)/day
a_ll = [lonr[0,3], 0.3, latr[3,0], latr[-3,0]]

for t in range(nt):
    print('day %i' % t)
    
    fig = plt.figure(figsize=(10,6))
    
    dens = density[t,:,:]-density0#density anomaly
    ax0 = fig.add_subplot(1,3,1)
    
    # create interpolant
    upi = intp.interp2d(lonr_s[0, :], latr_s[:, 0], upd[t,:])
    vpi = intp.interp2d(lonr_s[0, :], latr_s[:, 0], vpd[t,:])
    uai = intp.interp2d(lonr_s[0, :], latr_s[:, 0], uad[t,:])
    vai = intp.interp2d(lonr_s[0, :], latr_s[:, 0], vad[t,:])
        
    rs = np.ma.masked_where(maskr==0, dens)
    p = ax0.pcolormesh(lonr, latr, rs, vmin=-vmax, vmax=vmax,
        cmap='RdYlBu_r', alpha=1.0)
    ax0.axis(a_ll)
    pfun.dar(ax0)
    ax0.set_xlabel('Longitude',fontweight='bold')
    ax0.set_ylabel('Latitude',fontweight='bold')
    ax0.set_title(ax0_title, fontweight='bold')
    cbar = plt.colorbar(p,orientation='horizontal',pad=0.1)
    cbar.ax.set_xlabel(r'$kg/m^{3}$')
    
    if t==0:
        # create regular grid
        aaa = ax0.axis()
        daax = aaa[1] - aaa[0]
        daay = aaa[3] - aaa[2]
        axrat = np.cos(np.deg2rad(aaa[2])) * daax / daay
        x = np.linspace(aaa[0], aaa[1], int(round(nngrid * axrat)))
        y = np.linspace(aaa[2], aaa[3], nngrid)
        xx, yy = np.meshgrid(x, y)
    
    # interpolate to regular grid
    uup = upi(x, y)
    vvp = vpi(x, y)
    uua = uai(x, y)
    vva = vai(x, y)
    mask = uup != 0

    # add prsgrd quiver
    ax1 = fig.add_subplot(1,3,2)
    ax1.pcolormesh(lonr, latr, rs, vmin=-vmax, vmax=vmax,
        cmap='RdYlBu_r', alpha=.3)
    ax1.axis(a_ll)
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude',fontweight='bold')
    ax1.get_yaxis().set_visible(False)
    ax1.set_title('PG quiver',fontweight='bold')
    
    # ax1.quiver(lonr_s[b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
    #             latr_s[b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
    #             dia_scale*u_prsgrd[t,b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
    #             dia_scale*v_prsgrd[t, b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
    #             color='black',units='y',scale=qscale,scale_units='y',width=0.01)
    ax1.quiver(xx[mask], yy[mask], dia_scale*uup[mask], dia_scale*vvp[mask],color='k',units='y',scale=qscale,scale_units='y',width=0.01)
    #quiver legend
    ax1.quiver([0.02],[44.5],[quiverleg],[quiverleg],
                color='black',units='y',scale=qscale,scale_units='y',width=0.01)
    ax1.text(0.02,44.4,'%.0f cm/s\nper day' % quiverleg,fontsize=8)
    
    #add ageo pressure quiver
    ax2 = fig.add_subplot(1,3,3)
    ax2.pcolormesh(lonr, latr, rs, vmin=-vmax, vmax=vmax,
        cmap='RdYlBu_r', alpha=.3)
    ax2.axis(a_ll)
    pfun.dar(ax2)
    ax2.set_xlabel('Longitude',fontweight='bold')
    ax2.get_yaxis().set_visible(False)
    ax2.set_title('PG + Cor quiver',fontweight='bold')
    
    # ax2.quiver(lonr_s[b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
    #             latr_s[b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
    #             dia_scale*u_ageo[t, b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
    #             dia_scale*v_ageo[t, b0['latr0']:b0['latr1']:ss,b0['lonr0']:b0['lonr1']:ss],
    #             color='black',units='y',scale=qscale,scale_units='y',width=0.01)
    
    ax2.quiver(xx[mask], yy[mask], dia_scale*uua[mask], dia_scale*vva[mask],color='k',units='y',scale=qscale,scale_units='y',width=0.01)
    #quiver legend
    ax2.quiver([0.02],[44.5],[quiverleg],[quiverleg],
                color='black',units='y',scale=qscale,scale_units='y',width=0.01)
    ax2.text(0.02,44.4,'%.0f cm/s\nper day' % quiverleg,fontsize=8)

    dts = str(t+1).zfill(3)
    titlestr = 'Day %s' % str(t+1)
    outname = outdir + '/pic_%s.png' % dts
    plt.suptitle(titlestr, fontsize = 14, fontweight = "bold")
    plt.savefig(outname)
    
    plt.close()

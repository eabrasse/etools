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

import pickle
import numpy as np
import netCDF4 as nc

def daily_filt(vector):
    #assume daily data, as particle track data is
    LEN = len(vector)
    LEN_filt0 = np.floor(LEN/24)
    LEN_filt = int(LEN_filt0)
    vec_filtered = np.zeros(LEN_filt)
    for i in range(LEN_filt):
        initial = i * 24
        final = (i+1)*24
        value_filtered = np.mean(vector[initial:final])
        vec_filtered[i] = value_filtered
    return vec_filtered

#Establish lists
color_list = {'hadv':'cornflowerblue',
            'vadv':'seagreen',
            'cor':'goldenrod',
            'vvisc':'orchid',
            'accel':'indigo',
            'prsgrd':'firebrick',
            'cap':'coral',
            'ageo':'yellowgreen',
            'adv':'deepskyblue',
            'balance':'black'}
full_moment = 'hadv', 'vadv', 'cor', 'vvisc', 'prsgrd', 'accel', 'balance'
ageo_moment = 'adv', 'vvisc', 'ageo','accel'
bern_moment = 'vvisc', 'cap','accel'
bern_decomp = 'ageo', 'adv', 'cap'
Titles = {full_moment:'Full momentum balance',
        ageo_moment:'Momentum balance with geostrophic balance removed',
        bern_moment:'Momentum balance with CAP',
        bern_decomp:'Decomposition of CAP'}

list_choice = ageo_moment
    
AVGfn = 'momentum_balanceAVG_vars.p'
DIAfn = 'momentum_balanceDIA_vars.p'
trackfn = '/pmr4/eab32/LiveOcean_output/tracks/ae1_ndiv12_3d_turb/release_2017.02.14.nc'

AVG = pickle.load(open(AVGfn,'rb'))
DIA = pickle.load(open(DIAfn,'rb'))
for key in AVG.keys():
    locals()[key] = AVG[key]

# Use particle tracks to determine how many days we care about
# GET PARTICLE PATHS
dst = nc.Dataset(trackfn)
tsalt0 = dst['salt'][:]
tlon0 = dst['lon'][:]
tlat0 = dst['lat'][:]
tz0 = dst['z'][:]
tot = dst['ot'][:]
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

#CHOOSE INDIVIDUAL PARTICLE
particlechoice = 500
tlon = tlon1[:,particlechoice]
tlat = tlat1[:,particlechoice]
tz = tz1[:,particlechoice]

#should do some sort of smoothing
tlon_sm = daily_filt(tlon)
tlat_sm = daily_filt(tlat)
tz_sm = daily_filt(tz)
tot_sm = daily_filt(tot)
NT = len(tot_sm)

channelcross = np.where(tlon_sm>0)[0][0]

outdir = 'momentum_balance'
Lfun.make_dir(outdir,clean=True)

a_ll = [-0.65, 0.3, 44.7, 45.3]

for i in range(channelcross):
    
    fig=plt.figure(figsize=(9,7))
    ax = plt.gca()
    ax.axis(a_ll)
    
    #plot background
    sl = np.ma.masked_where(salt[45+i,:,:] > 40,salt[45+i,:,:])
    ax.pcolormesh(lonr,latr,sl,vmax=31,vmin=25,cmap='rainbow',alpha=0.35)
    
    #get and plot momentum forcings using quiver
    Du = dict()
    Dv = dict()
    for key in sorted(DIA):
        if key[0]=='u':
            #all keys are 3D, helpfully
            indx = zfun.find_nearest_ind(lonu[0,:],tlon_sm[i])
            indy = zfun.find_nearest_ind(latu[:,0],tlat_sm[i])
            indz = zfun.find_nearest_ind(z[45+i,:,indy,indx+1],tz_sm[i]) #indx+1 because zr is on rho grid
            Du[key[2:]] = DIA[key][45+i,indz,indy,indx]
        if key[0]=='v':
            indx = zfun.find_nearest_ind(lonv[0,:],tlon_sm[i])
            indy = zfun.find_nearest_ind(latv[:,0],tlat_sm[i])
            indz = zfun.find_nearest_ind(z[45+i,:,indy+1,indx],tz_sm[i]) #indy+1 because zr is on rho grid
            Dv[key[2:]] = DIA[key][45+i,indz,indy,indx]
    Du['ageo'] = Du['prsgrd']+Du['cor']
    Dv['ageo'] = Dv['prsgrd']+Dv['cor']
    Du['adv'] = Du['hadv']+Du['vadv']
    Dv['adv'] = Dv['hadv']+Dv['vadv']
    Du['cap'] = Du['ageo']+Du['adv']
    Dv['cap'] = Dv['ageo']+Dv['adv']
    Du['balance'] = Du['cap'] + Du['vvisc'] - Du['accel']
    Dv['balance'] = Dv['cap'] + Dv['vvisc'] - Dv['accel']

    scl = 0.000003
    count=0
    for key in list_choice:
        ax.quiver(tlon_sm[i],tlat_sm[i],Du[key],Dv[key],scale=scl,color=color_list[key],linewidth=0.5)
        ax.quiver(tlon_sm[i],tlat_sm[i],Du[key],Dv[key],scale=scl, edgecolor='k', facecolor='None', linewidth=.5)
        lattxt = 44.9 - 0.05 * count
        lontxt = 0.05
        ax.text(lontxt,lattxt,key,color=color_list[key],fontweight='bold')
        count=count+1

    #plot current particle position over quiver
    ax.plot(tlon_sm[i],tlat_sm[i],marker='o',color='black',markerfacecolor='green',markersize=10.0)


    titlestr = Titles[list_choice] + '\nParticle %i/%i' % (particlechoice, np.sum(tma))
    ax.set_title(titlestr)
    
    plotname = outdir + '/plot%04i.png' % i
    plt.savefig(plotname)
    plt.close('all')

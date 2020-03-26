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


trackfn1 = '/pmr4/eab32/LiveOcean_output/tracks/ae1_ndiv1_3d_turb/release_2017.02.14.nc'
trackfn12 = '/pmr4/eab32/LiveOcean_output/tracks/ae1_ndiv12_3d_turb/release_2017.02.14.nc'
dst1 = nc4.Dataset(trackfn1)
dst12 = nc4.Dataset(trackfn12)

#get particle positions
#cs0 = dst['cs'][:]
#z0 = dst['z'][:]
salt1 = dst1['salt'][:]
salt12 = dst12['salt'][:]
ot = dst1['ot'][:]

time = (ot-ot[0])/3600
time = time/24 + 45

cmap = cmo.cm.deep

nt = 75*24

# develop a filter for only particles that start at PSU =>30 and end at PSU <30
# start_as_ocean = salt[0,:] > 29.9
# felloff = np.zeros(len(start_as_ocean),dtype=bool)
# end_as_estuary = np.zeros(len(start_as_ocean),dtype=bool)
# end_as_channel = np.zeros(len(start_as_ocean),dtype=bool)
# for i in np.arange(len(end_as_channel)):
#     end_as_channel[i] = np.any(salt[:,i]<15.0) #stricter here
#     end_as_estuary[i] = np.any(salt[:,i]<25.0)
#
# tma0 = start_as_ocean * end_as_channel * end_as_estuary
# tma1 = start_as_ocean * ~end_as_channel * end_as_estuary
# tma2 = start_as_ocean * ~end_as_channel * ~end_as_estuary


fig = plt.figure(figsize=(12,5))
ax1 = fig.add_subplot(2,1,1)
#ax1.set_xlabel('Time (day)',fontweight='bold')
ax1.get_xaxis().set_visible(False)
ax1.set_ylabel('Salt (PSU)',fontweight='bold')
ax1.set_ylim([0,30])
# plt.plot(time[::24],salt[::24,tma2],alpha=0.5,color='royalblue')
# plt.plot(time[::24],salt[::24,tma1],alpha=0.5,color='seagreen')
# plt.plot(time[::24],salt[::24,tma0],alpha=0.5,color='goldenrod')
ax1.plot(time[::24],salt1[::24,:],alpha=0.5)
#plt.suptitle('Tracked from day 45 to day 90', fontsize = 14, fontweight = "bold")
ax1.set_title('ndiv=1')

ax2 = fig.add_subplot(2,1,2)
ax2.set_xlabel('Time (day)',fontweight='bold')
ax2.set_ylabel('Salt (PSU)',fontweight = 'bold')
ax2.set_ylim([0,30])
ax2.plot(time[::24],salt12[::24,:],alpha=0.5)
ax2.set_title('ndiv=12') 

plt.savefig('salinitytracks.png')
dst1.close()
dst12.close()

"""
Calculate TEF variables within a shelf box, maintaining temporal resolution
but averaging in space
"""

import os
import sys
import netCDF4 as nc
import numpy as np

pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
ep = os.path.abspath('../../etools')
if ep not in sys.path:
    sys.path.append(ep)
import efun
import zfun
import zrfun
import Lfun
import pickle

Ldir = Lfun.Lstart('aestus2', 'unstrat')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'

Ldirroms = '/pmr4/eab32/LiveOcean_ROMS/'

dir0 = Ldirroms + 'output/' + Ldir['gtagex'] + '/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

#testing
#f_list = f_list[:90] # control number of days

nt = len(f_list)*24

# initialize result arrays for history files
# (one more than the number of ave and dia files)
# these are all volume integrals with tidal variability
V_arr = np.nan * np.ones(nt+1)
Salt_arr = np.nan * np.ones(nt+1)
# initialize result arrays for average files
# 0 and 1 mean ocean and river ends
T_arr = np.nan + np.ones(nt) # ocean time

# initialize intermediate results arrays for TEF quantities
sedges = np.linspace(0, 32, 1001) # original was 35*20 + 1
sbins = sedges[:-1] + np.diff(sedges)/2
ns = len(sbins) # number of salinity bins

# loop over hours
tt = 0 # history
tta = 0 # averages
for f_dir in f_list:
    print(str(tta))
    
    # get volume-integrated quantities from the history files
    h_list = os.listdir(dir0 + f_dir)
    h_list.sort()
    h_list = [x for x in h_list if x[:9]=='ocean_his']
    if tt == 0:
        pass
    else:
        h_list = h_list[1:] # drop the zero hour for all but the first day
        
    for hi in h_list:
        fn = dir0 + f_dir + '/' + hi
        ds = nc.Dataset(fn)
        if tt == 0:
            
            Ind = efun.box_inds(-0.6, 1.5, 44.5, 45.5, fn)
            for vn in Ind.keys():
                globals()[vn] = Ind[vn]
        
            G, S, T = zrfun.get_basic_info(fn)
            h = G['h'][latr0:latr1+1, lonr0:lonr1+1]
            dx = G['DX'][latr0:latr1+1, lonr0:lonr1+1]
            dy = G['DY'][latr0:latr1+1, lonr0:lonr1+1]
            da = dx*dy

            ny, nx = h.shape
            
            tef_q_mouth0 = np.zeros((ns, nt, ny))
            tef_q_plume = np.zeros((ns, nt, nx))
            tef_q_shelf = np.zeros((ns, nt, ny))
            tef_q_south = np.zeros((ns, nt, nx))
            tef_qs_mouth0 = np.zeros((ns, nt, ny))
            tef_qs_plume = np.zeros((ns, nt, nx))
            tef_qs_shelf = np.zeros((ns, nt, ny))
            tef_qs_south = np.zeros((ns, nt, nx))

            lon_shelf = ds['lon_u'][latu0:latu1+1, lonu0]
            lat_shelf = ds['lat_u'][latu0:latu1+1, lonu0]            
            lon_mouth = ds['lon_u'][latu0:latu1+1, lonu1]
            lat_mouth = ds['lat_u'][latu0:latu1+1, lonu1]
            lon_south = ds['lon_v'][latv0, lonv0:lonv1+1]
            lat_south = ds['lat_v'][latv0, lonv0:lonv1+1]
            lon_plume = ds['lon_v'][latv1, lonv0:lonv1+1]
            lat_plume = ds['lat_v'][latv1, lonv0:lonv1+1]
                        
        zeta = ds['zeta'][0, latr0:latr1+1, lonr0:lonr1+1].squeeze()
        salt = ds['salt'][0, :, latr0:latr1+1, lonr0:lonr1+1].squeeze()
        zr, zw = zrfun.get_z(h, zeta, S)
        dzr = np.diff(zw, axis=0)
        
        V = np.sum(da.reshape((1,ny,nx))*dzr) # volume
        Salt = np.sum(da.reshape((1,ny,nx))*dzr*salt) # net salt
        
        # store results
        V_arr[tt] = V
        Salt_arr[tt] = Salt
        
        tt += 1
    
    # get arrays for flux claculations from the averages
    a_list = os.listdir(dir0 + f_dir)
    a_list.sort()
    a_list = [x for x in a_list if x[:9]=='ocean_avg']
    for ai in a_list:
        fn = dir0 + f_dir + '/' + ai
        ds = nc.Dataset(fn)
        T_arr[tta] = ds['ocean_time'][:].squeeze()

        # getting fluxes of volume and salt
        dq_shelf = ds['Huon'][0, : , latu0:latu1+1, lonu0].squeeze()
        dq_mouth = ds['Huon'][0, : , latu0:latu1+1, lonu1].squeeze()
        dq_south = ds['Hvom'][0, : , latv0, lonv0:lonv1+1].squeeze()
        dq_plume = ds['Hvom'][0, : , latv1, lonv0:lonv1+1].squeeze()
        dqs_shelf = ds['Huon_salt'][0, : , latu0:latu1+1, lonu0].squeeze()
        dqs_mouth = ds['Huon_salt'][0, : , latu0:latu1+1, lonu1].squeeze()
        dqs_south = ds['Hvom_salt'][0, : , latv0, lonv0:lonv1+1].squeeze()
        dqs_plume = ds['Hvom_salt'][0, : , latv1, lonv0:lonv1+1].squeeze()

        # then get the salinity averaged onto the u-grid on both open boundaries
        s_shelf = (ds['salt'][0, :, latr0:latr1+1, lonr0-1].squeeze() + ds['salt'][0, :, latr0:latr1+1, lonr0].squeeze())/2
        s_mouth = (ds['salt'][0, :, latr0:latr1+1, lonr1-1].squeeze() + ds['salt'][0, :, latr0:latr1+1, lonr1].squeeze())/2
        s_south = (ds['salt'][0, :, latr0-1, lonr0:lonr1+1].squeeze() + ds['salt'][0, :, latr0, lonr0:lonr1+1].squeeze())/2
        s_plume = (ds['salt'][0, :, latr1-1, lonr0:lonr1+1].squeeze() + ds['salt'][0, :, latr1, lonr0:lonr1+1].squeeze())/2

        # TEF variables
        # which are also area integrals at ocean and river ends
        
        for yy in range(ny):
            s_shelf0 = s_shelf[:,yy].squeeze()
            dq_shelf0 = dq_shelf[:,yy].squeeze()
            dqs_shelf0 = dqs_shelf[:,yy].squeeze()
            try:
                s_shelf00 = s_shelf0[s_shelf0.mask==False] # flattens the array
                dq_shelf00 = dq_shelf0[dq_shelf0.mask==False]
                dqs_shelf00 = dqs_shelf0[dqs_shelf0.mask==False]
            except:
                s_shelf00 = s_shelf0.flatten()
                dq_shelf00 = dq_shelf0.flatten()
                dqs_shelf00 = dqs_shelf0.flatten()
            inds = np.digitize(s_shelf00, sedges, right=True)
            dq_shelf00 = dq_shelf00.squeeze()
            dqs_shelf00 = dqs_shelf00.squeeze()
            counter = 0
            for ii in inds:
                tef_q_shelf[ii-1,tta,yy] += dq_shelf00[counter]
                tef_qs_shelf[ii-1,tta,yy] += dqs_shelf00[counter]
                counter += 1
            #
            s_mouth0 = s_mouth[:,yy].squeeze()
            dq_mouth0 = dq_mouth[:,yy].squeeze()
            dqs_mouth0 = dqs_mouth[:,yy].squeeze()
            
            try:
                s_mouth00 = s_mouth0[s_mouth0.mask==False]
                dq_mouth00 = dq_mouth0[dq_mouth0.mask==False]
                dqs_mouth00 = dqs_mouth0[dqs_mouth0.mask==False]

            except:
                s_mouth00 = s_mouth0.flatten()
                dq_mouth00 = dq_mouth0.flatten()
                dqs_mouth00 = dqs_mouth0.flatten()

            inds = np.digitize(s_mouth00, sedges, right=True)
            dq_mouth00 = dq_mouth00.squeeze()
            dqs_mouth00 = dqs_mouth00.squeeze()
            counter = 0
            for ii in inds:
                # at each time step these are vectors of hourly transport in
                # salinity bins (centered at sbins)
                tef_q_mouth0[ii-1,tta,yy] += dq_mouth00[counter]
                tef_qs_mouth0[ii-1,tta,yy] += dqs_mouth00[counter]
                counter += 1
    
        #
        for xx in range(nx):
            s_south0 = s_south[:,xx].squeeze()
            dq_south0 = dq_south[:,xx].squeeze()
            dqs_south0 = dqs_south[:,xx].squeeze()
            
            try:
                s_south00 = s_south0[s_south0.mask==False]
                dq_south00 = dq_south0[dq_south0.mask==False]
                dqs_south00 = dqs_south0[dqs_south0.mask==False]
            except:
                s_south00 = s_south0.flatten()
                dq_south00 = dq_south0.flatten()
                dqs_south00 = dqs_south0.flatten()
            inds = np.digitize(s_south00, sedges, right=True)
            dq_south00 = dq_south00.squeeze()
            dqs_south00 = dqs_south00.squeeze()
            counter = 0
            for ii in inds:
                # at each time step these are vectors of hourly transport in
                # salinity bins (centered at sbins)
                tef_q_south[ii-1,tta,xx] += dq_south00[counter]
                tef_qs_south[ii-1,tta,xx] += dqs_south00[counter]
                counter += 1
            #
            s_plume0 = s_plume[:,xx].squeeze()
            dq_plume0 = dq_plume[:,xx].squeeze()
            dqs_plume0 = dqs_plume[:,xx].squeeze()
            try:
                s_plume00 = s_plume0[s_plume0.mask==False]
                dq_plume00 = dq_plume0[dq_plume0.mask==False]
                dqs_plume00 = dqs_plume0[dqs_plume0.mask==False]
            except:
                s_plume00 = s_plume0.flatten()
                dq_plume00 = dq_plume0.flatten()
                dqs_plume00 = dqs_plume0.flatten()
            inds = np.digitize(s_plume00, sedges, right=True)
            dq_plume00 = dq_plume00.squeeze()
            dqs_plume00 = dqs_plume00.squeeze()
            counter = 0
            for ii in inds:
                # at each time step these are vectors of hourly transport in
                # salinity bins (centered at sbins)
                tef_q_plume[ii-1,tta,xx] += dq_plume00[counter]
                tef_qs_plume[ii-1,tta,xx] += dqs_plume00[counter]
                counter += 1
    
        ds.close()
    
        tta += 1

#%% TEF processing
#sum spatially
tef_q_shelfbox = -np.sum(tef_q_plume,axis=2)+np.sum(tef_q_shelf,axis=2)+np.sum(tef_q_south,axis=2)
tef_qs_shelfbox = -np.sum(tef_qs_plume,axis=2)+np.sum(tef_qs_shelf,axis=2)+np.sum(tef_qs_south,axis=2)
tef_q_mouth = np.sum(tef_q_mouth0,axis=2)
tef_qs_mouth = np.sum(tef_qs_mouth0,axis=2)

# first form tidal averages
tef_q_shelfbox_lp = np.nan * np.ones_like(tef_q_shelfbox)
tef_q_mouth_lp = np.nan * np.ones_like(tef_q_mouth)
tef_qs_shelfbox_lp = np.nan * np.ones_like(tef_qs_shelfbox)
tef_qs_mouth_lp = np.nan * np.ones_like(tef_qs_mouth)
for ii in range(ns):
    tef_q_shelfbox_lp[ii,:] = zfun.filt_godin(tef_q_shelfbox[ii,:])
    tef_q_mouth_lp[ii,:] = zfun.filt_godin(tef_q_mouth[ii,:])
    tef_qs_shelfbox_lp[ii,:] = zfun.filt_godin(tef_qs_shelfbox[ii,:])
    tef_qs_mouth_lp[ii,:] = zfun.filt_godin(tef_qs_mouth[ii,:])


# start by making the low-passed flux arrays sorted
# from high to low salinity
# np.flipud equivalent to m[::-1,:,:], according to documentation
# I prefer that as I know which axis I'm using
rq_shelfbox = np.flipud(tef_q_shelfbox_lp)
rqs_shelfbox = np.flipud(tef_qs_shelfbox_lp)
rq_mouth = np.flipud(tef_q_mouth_lp)
rqs_mouth = np.flipud(tef_qs_mouth_lp)

qcs_mouth = np.cumsum(rq_mouth, axis=0)
qcs_shelfbox = np.cumsum(rq_shelfbox,axis=0)

imax_m = np.argmax(qcs_mouth,axis=0)
imax_sb = np.argmax(qcs_shelfbox,axis=0)

Qin_shelfbox = np.zeros(nt)
QSin_shelfbox = np.zeros(nt)
Qout_shelfbox = np.zeros(nt)
QSout_shelfbox = np.zeros(nt)

Qin_mouth = np.zeros(nt)
QSin_mouth = np.zeros(nt)
Qout_mouth = np.zeros(nt)
QSout_mouth = np.zeros(nt)

for tt in range(nt):
    Qin_shelfbox[tt] = np.sum(rq_shelfbox[:imax_sb[tt], tt])
    Qout_shelfbox[tt] = np.sum(rq_shelfbox[imax_sb[tt]:, tt])
    QSin_shelfbox[tt] = np.sum(rqs_shelfbox[:imax_sb[tt], tt])
    QSout_shelfbox[tt] = np.sum(rqs_shelfbox[imax_sb[tt]:, tt])

    Qin_mouth[tt] = np.sum(rq_mouth[:imax_m[tt], tt])
    Qout_mouth[tt] = np.sum(rq_mouth[imax_m[tt]:, tt])
    QSin_mouth[tt] = np.sum(rqs_mouth[:imax_m[tt], tt])
    QSout_mouth[tt] = np.sum(rqs_mouth[imax_m[tt]:, tt])
    
# form derived quantities
Sin_shelfbox = QSin_shelfbox/Qin_shelfbox
Sout_shelfbox = QSout_shelfbox/Qout_shelfbox
Sin_mouth = QSin_mouth/Qin_mouth
Sout_mouth = QSout_mouth/Qout_mouth


D = dict()
D_list = ['T_arr','V_arr','Salt_arr',
    'Qin_shelfbox','Qout_shelfbox','Sin_shelfbox','Sout_shelfbox','QSin_shelfbox','QSout_shelfbox',
    'Qin_mouth','Qout_mouth','Sin_mouth','Sout_mouth','QSin_mouth','QSout_mouth']

for vn in D_list:
    D[vn] = locals()[vn]

out_dir = Ldir['parent'] + 'etools/TEFdy/'
Lfun.make_dir(out_dir,clean='False')
out_fn = out_dir + 'TEFbalancev3_aestus2_unstrat.p'
pickle.dump(D, open(out_fn, 'wb'))

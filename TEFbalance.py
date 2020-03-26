"""
Modified version of svar to calculate TEF across 20 cross sections
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

Ldir = Lfun.Lstart('aestus1', 'base')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'

Ldirroms = '/pmr4/eab32/LiveOcean_ROMS/'

dir0 = Ldirroms + 'output/' + Ldir['gtagex'] + '/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']

#testing
f_list = f_list[45:120] # control number of days

nt = len(f_list)*24

# initialize intermediate results arrays for TEF quantities
sedges = np.linspace(0, 32, 1001) # original was 35*20 + 1
sbins = sedges[:-1] + np.diff(sedges)/2
ns = len(sbins) # number of salinity bins

# loop over hours
tta = 0 # averages
for f_dir in f_list:
    print(str(tta))
    
    # get arrays for flux claculations from the averages
    a_list = os.listdir(dir0 + f_dir)
    a_list.sort()
    a_list = [x for x in a_list if x[:9]=='ocean_avg']
    for ai in a_list:
        fn = dir0 + f_dir + '/' + ai
        ds = nc.Dataset(fn)
        
        if tta == 0:
            
            Ind = efun.box_inds(-0.6, 0.02, 44.5, 45.5, fn)
            for vn in Ind.keys():
                globals()[vn] = Ind[vn]
        
            G, S, T = zrfun.get_basic_info(fn)
            h = G['h'][latr0:latr1+1, lonr0:lonr1+1]

            ny, nx = h.shape
            
            tef_q_mouth = np.zeros((ns, nt, ny))
            tef_q_plume = np.zeros((ns, nt, nx))
            tef_q_shelf = np.zeros((ns, nt, ny))
            tef_q_south = np.zeros((ns, nt, nx))
            tef_qs_mouth = np.zeros((ns, nt, ny))
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
                tef_q_mouth[ii-1,tta,yy] += dq_mouth00[counter]
                tef_qs_mouth[ii-1,tta,yy] += dqs_mouth00[counter]
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
# first form tidal averages
tef_q_plume_mean = -np.mean(tef_q_plume,axis=1)
tef_qs_plume_mean = -np.mean(tef_qs_plume,axis=1)
tef_q_shelf_mean = np.mean(tef_q_shelf,axis=1)
tef_qs_shelf_mean = np.mean(tef_qs_shelf,axis=1)
tef_q_south_mean = np.mean(tef_q_south,axis=1)
tef_qs_south_mean = np.mean(tef_qs_south,axis=1)
tef_q_mouth_mean = np.mean(tef_q_mouth,axis=1)
tef_qs_mouth_mean = np.mean(tef_qs_mouth,axis=1)

# start by making the low-passed flux arrays sorted
# from high to low salinity
rq_plume = np.flipud(tef_q_plume_mean)
rqs_plume = np.flipud(tef_qs_plume_mean)
rq_shelf = np.flipud(tef_q_shelf_mean)
rqs_shelf = np.flipud(tef_qs_shelf_mean)
rq_mouth = np.flipud(tef_q_mouth_mean)
rqs_mouth = np.flipud(tef_qs_mouth_mean)
rq_south = np.flipud(tef_q_south_mean)
rqs_south = np.flipud(tef_qs_south_mean)
# then form the cumulative sum (the functon Q(s))
qcs_plume = np.cumsum(rq_plume, axis=0)
qcs_shelf = np.cumsum(rq_shelf, axis=0)
qcs_mouth = np.cumsum(rq_mouth, axis=0)
qcs_south = np.cumsum(rq_south, axis=0)
# and find its maximum: this is Qin, and the salinity
# at which it occurs is the "dividing salinity" between
# inflow and outflow that we will use to calculate
# all TEF quantities
imax_plume = np.argmax(qcs_plume, axis=0)
imax_shelf = np.argmax(qcs_shelf, axis=0)
imax_south = np.argmax(qcs_south, axis=0)
imax_mouth = np.argmax(qcs_mouth, axis=0)

Qin_plume = np.zeros(nx)
QSin_plume = np.zeros(nx)
Qout_plume = np.zeros(nx)
QSout_plume = np.zeros(nx)

Qin_shelf = np.zeros(ny)
QSin_shelf = np.zeros(ny)
Qout_shelf = np.zeros(ny)
QSout_shelf = np.zeros(ny)

Qin_south = np.zeros(nx)
QSin_south = np.zeros(nx)
Qout_south = np.zeros(nx)
QSout_south = np.zeros(nx)

Qin_mouth = np.zeros(ny)
QSin_mouth = np.zeros(ny)
Qout_mouth = np.zeros(ny)
QSout_mouth = np.zeros(ny)

for y in range(ny):
    Qin_shelf[y] = np.sum(rq_shelf[:imax_shelf[y], y], axis=0)
    Qout_shelf[y] = np.sum(rq_shelf[imax_shelf[y]:, y], axis=0)
    QSin_shelf[y] = np.sum(rqs_shelf[:imax_shelf[y], y], axis=0)
    QSout_shelf[y] = np.sum(rqs_shelf[imax_shelf[y]:, y], axis=0)
    
    Qin_mouth[y] = np.sum(rq_mouth[:imax_mouth[y], y], axis=0)
    Qout_mouth[y] = np.sum(rq_mouth[imax_mouth[y]:, y], axis=0)
    QSin_mouth[y] = np.sum(rqs_mouth[:imax_mouth[y], y], axis=0)
    QSout_mouth[y] = np.sum(rqs_mouth[imax_mouth[y]:, y], axis=0)

# form derived quantities
Sin_shelf = QSin_shelf/Qin_shelf
Sout_shelf = QSout_shelf/Qout_shelf
Sin_mouth = QSin_mouth/Qin_mouth
Sout_mouth = QSout_mouth/Qout_mouth
    
for x in range(nx):
    Qin_plume[x] = np.sum(rq_plume[:imax_plume[x], x], axis=0)
    Qout_plume[x] = np.sum(rq_plume[imax_plume[x]:, x], axis=0)
    QSin_plume[x] = np.sum(rqs_plume[:imax_plume[x], x], axis=0)
    QSout_plume[x] = np.sum(rqs_plume[imax_plume[x]:, x], axis=0)
    
    Qin_south[x] = np.sum(rq_south[:imax_south[x], x], axis=0)
    Qout_south[x] = np.sum(rq_south[imax_south[x]:, x], axis=0)
    QSin_south[x] = np.sum(rqs_south[:imax_south[x], x], axis=0)
    QSout_south[x] = np.sum(rqs_south[imax_south[x]:, x], axis=0)

# form derived quantities
Sin_plume = QSin_plume/Qin_plume
Sout_plume = QSout_plume/Qout_plume
Sin_south = QSin_south/Qin_south
Sout_south = QSout_south/Qout_south

D = dict()
D_list = ['nx','ny',
    'lon_shelf','lat_shelf','Qin_shelf','Qout_shelf','Sin_shelf','Sout_shelf','QSin_shelf','QSout_shelf',
    'lon_mouth','lat_mouth','Qin_mouth','Qout_mouth','Sin_mouth','Sout_mouth','QSin_mouth','QSout_mouth',
    'lon_south','lat_south','Qin_south','Qout_south','Sin_south','Sout_south','QSin_south','QSout_south',
    'lon_plume','lat_plume','Qin_plume','Qout_plume','Sin_plume','Sout_plume','QSin_plume','QSout_plume']

for vn in D_list:
    D[vn] = locals()[vn]

out_dir = Ldir['parent'] + 'etools/TEFdy/'
Lfun.make_dir(out_dir,clean='False')
out_fn = out_dir + 'TEFbalance.p'
pickle.dump(D, open(out_fn, 'wb'))

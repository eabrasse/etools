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

def TEF_box(lons, lats, t):
    
    Ldir = Lfun.Lstart('aestus1', 'base')
    Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'

    Ldirroms = '/pmr4/eab32/LiveOcean_ROMS/'

    dir0 = Ldirroms + 'output/' + Ldir['gtagex'] + '/'

    f_list = os.listdir(dir0)
    f_list.sort()
    f_list = [x for x in f_list if x[0]=='f']


    if len(t)==1:
        f_list = f_list[t:]
    if len(t)==2:
        f_list = f_list[t[0]:t[1]+1] # control number of days
    else:
        print('t must be of length 1 or 2')
        print('Defaulting to use whole f_list')

    nt = len(f_list)*24

    # initialize intermediate results arrays for TEF quantities
    sedges = np.linspace(0, 35, 1001) # original was 35*20 + 1
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
                
                Ind = efun.box_inds(lons[0], lons[1], lats[0], lats[1], fn)
                for vn in Ind.keys():
                    globals()[vn] = Ind[vn]
            
                G, S, T = zrfun.get_basic_info(fn)
                h = G['h'][latr0:latr1+1, lonr0:lonr1+1]

                ny, nx = h.shape
                
                tef_q0 = np.zeros((ns, nt, ny))
                tef_q1 = np.zeros((ns, nt, ny))
                tef_q2 = np.zeros((ns, nt, ny))
                tef_q3 = np.zeros((ns, nt, ny))
                tef_qs0 = np.zeros((ns, nt, ny))
                tef_qs1 = np.zeros((ns, nt, ny))
                tef_qs2 = np.zeros((ns, nt, ny))
                tef_qs3 = np.zeros((ns, nt, ny))

            # getting fluxes of volume and salt
            dq0 = ds['Huon'][0, : , latu0:latu1+1, lonu0].squeeze()
            dq1 = ds['Huon'][0, : , latu0:latu1+1, lonu1].squeeze()
            dq2 = ds['Hvom'][0, : , latv0, lonv0:lonv1+1].squeeze()
            dq3 = ds['Hvom'][0, : , latv1, lonv0:lonv1+1].squeeze()
            dqs0 = ds['Huon_salt'][0, : , latu0:latu1+1, lonu0].squeeze()
            dqs1 = ds['Huon_salt'][0, : , latu0:latu1+1, lonu1].squeeze()
            dqs2 = ds['Hvom_salt'][0, : , latv0, lonv0:lonv1+1].squeeze()
            dqs3 = ds['Hvom_salt'][0, : , latv1, lonv0:lonv1+1].squeeze()

            # then get the salinity averaged onto the u-grid on both open boundaries
            s0 = (ds['salt'][0, :, latr0:latr1+1, lonr0-1].squeeze() + ds['salt'][0, :, latr0:latr1+1, lonr0].squeeze())/2
            s1 = (ds['salt'][0, :, latr0:latr1+1, lonr1].squeeze() + ds['salt'][0, :, latr0:latr1+1, lonr1+1].squeeze())/2
            s2 = (ds['salt'][0, :, latr0-1, lonr0:lonr1+1].squeeze() + ds['salt'][0, :, latr0, lonr0:lonr1+1].squeeze())/2
            s3 = (ds['salt'][0, :, latr1, lonr0:lonr1+1].squeeze() + ds['salt'][0, :, latr1+1, lonr0:lonr1+1].squeeze())/2

            # TEF variables
            # which are also area integrals at ocean and river ends
            
            for yy in range(ny):
                s00 = s0[:,yy].squeeze()
                dq00 = dq0[:,yy].squeeze()
                dqs00 = dqs0[:,yy].squeeze()
                try:
                    s000 = s00[s00.mask==False] # flattens the array
                    dq000 = dq00[dq00.mask==False]
                    dqs000 = dqs00[dqs00.mask==False]
                except:
                    s000 = s00.flatten()
                    dq000 = dq00.flatten()
                    dqs000 = dqs00.flatten()
                inds = np.digitize(s000, sedges, right=True)
                counter = 0
                for ii in inds:
                    print('ii ');print(ii)
                    print('counter ');print(counter)
                    print('yy ');print(yy)
                    print('tta ');print(tta)
                    tef_q0[ii-1,tta,yy] += dq000[counter]
                    tef_qs0[ii-1,tta,yy] += dqs000[counter]
                    counter += 1
                #
                s11 = s1[:,yy].squeeze()
                dq11 = dq1[:,yy].squeeze()
                dqs11 = dqs1[:,yy].squeeze()
                
                try:
                    s111 = s11[s11.mask==False]
                    dq111 = dq11[dq11.mask==False]
                    dqs111 = dqs11[dqs11.mask==False]

                except:
                    s111 = s11.flatten()
                    dq111 = dq11.flatten()
                    dqs111 = dqs11.flatten()

                inds = np.digitize(s111, sedges, right=True)
                counter = 0
                for ii in inds:
                    # at each time step these are vectors of hourly transport in
                    # salinity bins (centered at sbins)
                    tef_q1[ii-1,tta,yy] += dq111[counter]
                    tef_qs1[ii-1,tta,yy] += dqs111[counter]
                    counter += 1
        
            #
            for xx in range(nx):
                s22 = s2[:,xx].squeeze()
                dq22 = dq2[:,xx].squeeze()
                dqs22 = dqs2[:,xx].squeeze()
                
                try:
                    s222 = s22[dq22.mask==False]
                    dq222 = dq22[dq22.mask==False]
                    dqs222 = dqs22[dqs22.mask==False]
                except:
                    s222 = s22.flatten()
                    dq222 = dq22.flatten()
                    dqs222 = dqs22.flatten()
                inds = np.digitize(s222, sedges, right=True)
                counter = 0
                for ii in inds:
                    # at each time step these are vectors of hourly transport in
                    # salinity bins (centered at sbins)
                    tef_q2[ii-1,tta,xx] += dq222[counter]
                    tef_qs2[ii-1,tta,xx] += dqs222[counter]
                    counter += 1
                #
                s33 = s3[:,xx].squeeze()
                dq33 = dq3[:,xx].squeeze()
                dqs33 = dqs3[:,xx].squeeze()
                try:
                    s333 = s33[dq33.mask==False]
                    dq333 = dq33[dq33.mask==False]
                    dqs333 = dqs33[dqs33.mask==False]
                except:
                    s333 = s33.flatten()
                    dq333 = dq33.flatten()
                    dqs333 = dqs33.flatten()
                inds = np.digitize(s333, sedges, right=True)
                counter = 0
                for ii in inds:
                    # at each time step these are vectors of hourly transport in
                    # salinity bins (centered at sbins)
                    tef_q3[ii-1,tta,xx] += dq333[counter]
                    tef_qs3[ii-1,tta,xx] += dqs333[counter]
                    counter += 1
        
            ds.close()
        
            tta += 1

    # add together the results from the different branches of the shelf
    tef_q0 = np.concatenate((tef_q0, tef_q2, -tef_q3),axis=2)
    tef_qs0 = np.concatenate((tef_qs0, tef_qs2, -tef_qs3),axis=2)
    
    nxy = np.shape(tef_q0)[2]

    #%% TEF processing
    # first form tidal averages
    tef_q0_lp = np.mean(tef_q0,axis=1)
    tef_q1_lp = np.mean(tef_q1,axis=1)
    tef_qs0_lp = np.mean(tef_qs0,axis=1)
    tef_qs1_lp = np.mean(tef_qs1,axis=1)
    
    # start by making the low-passed flux arrays sorted
    # from high to low salinity
    rq0 = np.flipud(tef_q0_lp)
    rqs0 = np.flipud(tef_qs0_lp)
    # then form the cumulative sum (the functon Q(s))
    qcs = np.cumsum(rq0, axis=0)
    # and find its maximum: this is Qin, and the salinity
    # at which it occurs is the "dividing salinity" between
    # inflow and outflow that we will use to calculate
    # all TEF quantities
    imax = np.argmax(qcs, axis=0)
    Qin0 = np.zeros(nxy)
    QSin0 = np.zeros(nxy)
    Qout0 = np.zeros(nxy)
    QSout0 = np.zeros(nxy)
    for xy in range(nxy):
        Qin0[xy] = np.sum(rq0[:imax[xy], xy], axis=0)
        Qout0[xy] = np.sum(rq0[imax[xy]:, xy], axis=0)
        QSin0[xy] = np.sum(rqs0[:imax[xy], xy], axis=0)
        QSout0[xy] = np.sum(rqs0[imax[xy]:, xy], axis=0)
    # then fix masking so that the nan's from the low-pass are retained
    nmask = np.isnan(tef_q0_lp[0,:])
    Qin0[nmask] = np.nan
    QSin0[nmask] = np.nan
    Qout0[nmask] = np.nan
    QSout0[nmask] = np.nan

    # form derived quantities
    Sin0 = QSin0/Qin0
    Sout0 = QSout0/Qout0


    # same steps for the river end
    # OK in this case to do this the original way because
    # it is all freshwater.
    qin = tef_q1_lp.copy()
    qout = tef_q1_lp.copy()
    qsin = tef_qs1_lp.copy()
    qsout = tef_qs1_lp.copy()
    #
    qin[tef_q1_lp>0] = 0 # switch signs compared to open boundary 0
    qout[tef_q1_lp<0] = 0
    qsin[tef_q1_lp>0] = 0
    qsout[tef_q1_lp<0] = 0
    #
    Qin1 = qin.sum(axis=0)
    Qout1 = qout.sum(axis=0)
    QSin1 = qsin.sum(axis=0)
    QSout1 = qsout.sum(axis=0)
    #
    #Sin1 = QSin1/Qin1
    #Sout1 = QSout1/Qout1 # Qout1 = 0
    # Make a dictionary to return variables of interest efficiently
    
    D = dict()
    D_list = ['nxy',
        'Qin0', 'Qout0', 'Sin0', 'Sout0']#
        #'Qin1']
    for vn in D_list:
        D[vn] = locals()[vn]
    
    return D


def TEF_line(lon, t):
    
    Ldir = Lfun.Lstart('aestus1', 'base')
    Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'

    Ldirroms = '/pmr4/eab32/LiveOcean_ROMS/'

    dir0 = Ldirroms + 'output/' + Ldir['gtagex'] + '/'

    f_list = os.listdir(dir0)
    f_list.sort()
    f_list = [x for x in f_list if x[0]=='f']

    if len(t)==1:
        f_list = f_list[t:]
    if len(t)==2:
        f_list = f_list[t[0]:t[1]+1] # control number of days
    else:
        print('t must be of length 1 or 2')
        print('Defaulting to use whole f_list')

    # number of time steps for averages and diagnostics
    nt = len(f_list) * 24

    # initialize result arrays for average files
    # 0 and 1 mean ocean and river ends

    # initialize intermediate results arrays for TEF quantities
    sedges = np.linspace(0, 35, 1001) # original was 35*20 + 1
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
                
                Ind = efun.box_inds(lon, 1.5, 44.9, 45.1, fn)
                for vn in Ind.keys():
                    globals()[vn] = Ind[vn]
                    # print(vn)
            
                G, S, T = zrfun.get_basic_info(fn)
                h = G['h'][latr0:latr1+1, lonr0:lonr1+1]
                ny, nx = h.shape
                
                tef_q0 = np.zeros((ns, nt, ny))
                tef_q1 = np.zeros((ns, nt, ny))
                tef_qs0 = np.zeros((ns, nt, ny))
                tef_qs1 = np.zeros((ns, nt, ny))
                
            # getting fluxes of volume and salt
            dq0 = ds['Huon'][0, : , latu0:latu1+1, lonu0].squeeze()
            dq1 = ds['Huon'][0, : , latu0:latu1+1, lonu1+1].squeeze()
            dqs0 = ds['Huon_salt'][0, : , latu0:latu1+1, lonu0].squeeze()
            dqs1 = ds['Huon_salt'][0, : , latu0:latu1+1, lonu1+1].squeeze()

            # then get the salinity averaged onto the u-grid on both open boundaries
            s0 = (ds['salt'][0, :, latr0:latr1+1, lonr0-1].squeeze() + ds['salt'][0, :, latr0:latr1+1, lonr0].squeeze())/2
            s1 = (ds['salt'][0, :, latr0:latr1+1, lonr1].squeeze() + ds['salt'][0, :, latr0:latr1+1, lonr1+1].squeeze())/2
            
            # TEF variables
            # which are also area integrals at ocean and river ends
            # try:
            for yy in range(ny):
                s00 = s0[:,yy].squeeze()
                dq00 = dq0[:,yy].squeeze()
                dqs00 = dqs0[:,yy].squeeze()
            
                s000 = s00[dq00.mask==False] # flattens the array
                dq000 = dq00[dq00.mask==False]
                dqs000 = dqs00[dqs00.mask==False]

                inds = np.digitize(s000, sedges, right=True)
                counter = 0
                for ii in inds:
                    tef_q0[ii-1,tta,yy] += dq000[counter]
                    tef_qs0[ii-1,tta,yy] += dqs000[counter]
                    counter += 1

                s11 = s1[:,yy].squeeze()
                dq11 = dq1[:,yy].squeeze()
                dqs11 = dqs1[:,yy].squeeze()
                
                s111 = s11[s11.mask==False]
                dq111 = dq11[dq11.mask==False]
                dqs111 = dqs11[dqs11.mask==False]

                inds = np.digitize(s111, sedges, right=True)
                counter = 0
                for ii in inds:
                    # at each time step these are vectors of hourly transport in
                    # salinity bins (centered at sbins)
                    tef_q1[ii-1,tta,yy] += dq111[counter]
                    tef_qs1[ii-1,tta,yy] += dqs111[counter]
                    counter += 1
            ds.close()
        
            tta += 1

    #%% TEF processing
    # first form tidal averages
    tef_q0_lp = np.mean(tef_q0,axis=1)
    tef_q1_lp = np.mean(tef_q1,axis=1)
    tef_qs0_lp = np.mean(tef_qs0,axis=1)
    tef_qs1_lp = np.mean(tef_qs1,axis=1)
    
    # start by making the low-passed flux arrays sorted
    # from high to low salinity
    rq0 = np.flipud(tef_q0_lp)
    rqs0 = np.flipud(tef_qs0_lp)
    # then form the cumulative sum (the functon Q(s))
    qcs = np.cumsum(rq0, axis=0)
    # and find its maximum: this is Qin, and the salinity
    # at which it occurs is the "dividing salinity" between
    # inflow and outflow that we will use to calculate
    # all TEF quantities
    imax = np.argmax(qcs, axis=0)
    Qin0 = np.zeros(ny)
    QSin0 = np.zeros(ny)
    Qout0 = np.zeros(ny)
    QSout0 = np.zeros(ny)
    for yy in range(ny):
        Qin0[yy] = np.sum(rq0[:imax[yy], yy], axis=0)
        Qout0[yy] = np.sum(rq0[imax[yy]:, yy], axis=0)
        QSin0[yy] = np.sum(rqs0[:imax[yy], yy], axis=0)
        QSout0[yy] = np.sum(rqs0[imax[yy]:, yy], axis=0)
    # then fix masking so that the nan's from the low-pass are retained
    nmask = np.isnan(tef_q0_lp[0,:])
    Qin0[nmask] = np.nan
    QSin0[nmask] = np.nan
    Qout0[nmask] = np.nan
    QSout0[nmask] = np.nan
    # form derived quantities
    Sin0 = QSin0/Qin0
    Sout0 = QSout0/Qout0

    # same steps for the river end
    # OK in this case to do this the original way because
    # it is all freshwater.
    qin = tef_q1_lp.copy()
    qout = tef_q1_lp.copy()
    qsin = tef_qs1_lp.copy()
    qsout = tef_qs1_lp.copy()
    #
    qin[tef_q1_lp>0] = 0 # switch signs compared to open boundary 0
    qout[tef_q1_lp<0] = 0
    qsin[tef_q1_lp>0] = 0
    qsout[tef_q1_lp<0] = 0
    #
    Qin1 = qin.sum(axis=0)
    Qout1 = qout.sum(axis=0)
    QSin1 = qsin.sum(axis=0)
    QSout1 = qsout.sum(axis=0)
    #
    #Sin1 = QSin1/Qin1
    
    # Make a dictionary to return variables of interest efficiently
    D = dict()
    D_list = ['ny',
        'Qin0', 'Qout0', 'Sin0', 'Sout0']#,
        #'Qin1']
    for vn in D_list:
        D[vn] = locals()[vn]
    
    return D

"""
Generating a mean state picture of the ocean
"""

import os
import sys
import netCDF4 as nc
import numpy as np

pth = os.path.abspath('/pmr4/eab32/LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
ep = os.path.abspath('/pmr4/eab32/etools')
if ep not in sys.path:
    sys.path.append(ep)
import efun
import zfun
import zrfun
import Lfun
#import pdb
import pickle

import time
start = time.time()
Ldir = Lfun.Lstart('aestus1', 'base')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'

Ldirroms = '/pmr4/eab32/LiveOcean_ROMS/'

dir0 = Ldirroms + 'output/' + Ldir['gtagex'] + '/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']


f_list= f_list[45:91]

# number of time steps for averages and diagnostics
nt = len(f_list) * 24

# loop over hours
tt = 0 # history
tta = 0 # averages
for f_dir in f_list:
    print(str(tt))

    # get volume-integrated quantities from the history files
    fn_list = os.listdir(dir0 + f_dir)
    fn_list.sort()
    d_list = [x for x in fn_list if x[:9]=='ocean_dia']
    h_list = [x for x in fn_list if x[:9]=='ocean_his']
    list_len = len(d_list)
    #
    for i in range(list_len):
        
        fnd = dir0 + f_dir + '/' + d_list[i]
        fnh0 = dir0 + f_dir + '/' + h_list[i]
        fnh1 = dir0 + f_dir + '/' + h_list[i+1]
        dsd = nc.Dataset(fnd)
        dsh0 = nc.Dataset(fnh0)
        dsh1 = nc.Dataset(fnh1)
        
        if tt == 0:
            ot0 = dsh0['ocean_time'][:]
            ot1 = dsh1['ocean_time'][:]
            dt = ot1-ot0
            
            u_prsgrd = dsd['u_prsgrd'][0,:,:,:]
            NZ,NY,NX = np.shape(u_prsgrd)

            u_accel_avg = np.zeros((NZ, NY, NX))
            v_accel_avg = np.zeros((NZ, NY-1, NX+1))
            u_prsgrd_avg = np.zeros((NZ, NY, NX))
            v_prsgrd_avg = np.zeros((NZ, NY-1, NX+1))
            u_vvisc_avg = np.zeros((NZ, NY, NX))
            v_vvisc_avg = np.zeros((NZ, NY-1, NX+1))
            u_cor_avg = np.zeros((NZ, NY, NX))
            v_cor_avg = np.zeros((NZ, NY-1, NX+1))
            u_hadv_avg = np.zeros((NZ, NY, NX))
            v_hadv_avg = np.zeros((NZ, NY-1, NX+1))
            u_vadv_avg = np.zeros((NZ, NY, NX))
            v_vadv_avg = np.zeros((NZ, NY-1, NX+1))
            
            dudt_avg = np.zeros((NZ, NY, NX))
            dvdt_avg = np.zeros((NZ, NY-1, NX+1))
                    
        u_accel = dsd['u_accel'][0, :, :, :].squeeze()
        v_accel = dsd['v_accel'][0, :, :, :].squeeze()
        u_prsgrd = dsd['u_prsgrd'][0, :, :, :].squeeze()
        v_prsgrd = dsd['v_prsgrd'][0, :, :, :].squeeze()
        u_vvisc = dsd['u_vvisc'][0, :, :, :].squeeze()
        v_vvisc = dsd['v_vvisc'][0, :, :, :].squeeze()
        u_cor = dsd['u_cor'][0, :, :, :].squeeze()
        v_cor = dsd['v_cor'][0, :, :, :].squeeze()
        u_hadv = dsd['u_hadv'][0, :, :, :].squeeze()
        v_hadv = dsd['v_hadv'][0, :, :, :].squeeze()
        u_vadv = dsd['u_vadv'][0, :, :, :].squeeze()
        v_vadv = dsd['v_vadv'][0, :, :, :].squeeze()
        
        u0 = dsh0['u'][0,:,:,:].squeeze()
        u1 = dsh1['u'][0,:,:,:].squeeze()
        du = u1-u0
        dudt = du/dt
        v0 = dsh0['v'][0,:,:,:].squeeze()
        v1 = dsh1['v'][0,:,:,:].squeeze()
        dv = v1-v0
        dvdt = dv/dt
        
        u_accel_avg = u_accel_avg + u_accel
        v_accel_avg = v_accel_avg + v_accel
        u_prsgrd_avg = u_prsgrd_avg + u_prsgrd
        v_prsgrd_avg = v_prsgrd_avg + v_prsgrd
        u_vvisc_avg = u_vvisc_avg + u_vvisc
        v_vvisc_avg = v_vvisc_avg + v_vvisc
        u_cor_avg = u_cor_avg + u_cor
        v_cor_avg = v_cor_avg + v_cor
        u_hadv_avg = u_hadv_avg + u_hadv
        v_hadv_avg = v_hadv_avg + v_hadv
        u_vadv_avg = u_vadv_avg + u_vadv
        v_vadv_avg = v_vadv_avg + v_vadv
        
        dudt_avg = dudt_avg + dudt
        dvdt_avg = dvdt_avg + dvdt
    
        dsd.close()
        dsh0.close()
        dsh1.close()
        
        tt += 1

u_accel_avg = u_accel_avg/tt
v_accel_avg = v_accel_avg/tt
u_prsgrd_avg = u_prsgrd_avg/tt
v_prsgrd_avg = v_prsgrd_avg/tt
u_vvisc_avg = u_vvisc_avg/tt
v_vvisc_avg = v_vvisc_avg/tt
u_cor_avg = u_cor_avg/tt
v_cor_avg = v_cor_avg/tt
u_hadv_avg = u_hadv_avg/tt
v_hadv_avg = v_hadv_avg/tt
u_vadv_avg = u_vadv_avg/tt
v_vadv_avg = v_vadv_avg/tt
dudt_avg = dudt_avg/tt
dvdt_avg = dvdt_avg/tt

D = dict()
D['u_accel'] = u_accel_avg
D['v_accel'] = v_accel_avg
D['u_prsgrd'] = u_prsgrd_avg
D['v_prsgrd'] = v_prsgrd_avg
D['u_vvisc'] = u_vvisc_avg
D['v_vvisc'] = v_vvisc_avg
D['u_cor'] = u_cor_avg
D['v_cor'] = v_cor_avg
D['u_hadv'] = u_hadv_avg
D['v_hadv'] = v_hadv_avg
D['u_vadv'] = u_vadv_avg
D['v_vadv'] = v_vadv_avg
D['u_dt'] = dudt_avg
D['v_dt'] = dvdt_avg

out_fn = 'aestus1_base_ae1_average_dia.p'
pickle.dump( D, open(out_fn,'wb'))

end = time.time()
print(end-start)

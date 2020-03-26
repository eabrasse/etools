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
import pickle

Ldir = Lfun.Lstart('aestus2', 'unstrat')
Ldir['gtagex'] = Ldir['gtag'] + '_' + 'ae1'

Ldirroms = '/pmr4/eab32/LiveOcean_ROMS/'

dir0 = Ldirroms + 'output/' + Ldir['gtagex'] + '/'

f_list = os.listdir(dir0)
f_list.sort()
f_list = [x for x in f_list if x[0]=='f']


f_list= f_list[50:120]

# number of time steps for averages and diagnostics
nt = len(f_list) * 24

# loop over hours
tt = 0 # history
tta = 0 # averages
for f_dir in f_list:
    print(str(tt))

    # get volume-integrated quantities from the history files
    a_list = os.listdir(dir0 + f_dir)
    a_list.sort()
    a_list = [x for x in a_list if x[:9]=='ocean_avg']
    #
    for ai in a_list:
        fn = dir0 + f_dir + '/' + ai
        ds = nc.Dataset(fn)
        if tt == 0:
            
            salt = ds['salt'][0,:,:,:]
            NZ,NY,NX = np.shape(salt)

            salt_avg = np.zeros((NZ, NY, NX))
            u_avg = np.zeros((NZ, NY, NX-1))
            v_avg = np.zeros((NZ, NY-1, NX))
            w_avg = np.zeros((NZ+1, NY, NX))
            zeta_avg = np.zeros((NY, NX))
            rho_avg = np.zeros((NZ, NY, NX))
                    
        zeta = ds['zeta'][0, :, :].squeeze()
        salt = ds['salt'][0, :, :,:].squeeze()
        rho = ds['rho'][0,:,:,:].squeeze()
        u = ds['u'][0, :, :, :].squeeze()
        v = ds['v'][0, :, :, :].squeeze()
        w = ds['w'][0, :, :, :].squeeze()
        
        zeta_avg = zeta_avg + zeta
        salt_avg = salt_avg + salt
        rho_avg = rho_avg + rho
        u_avg = u_avg + u
        v_avg = v_avg + v
        w_avg = w_avg + w
    
        ds.close()
        tt += 1

zeta_avg = zeta_avg/tt
salt_avg = salt_avg/tt
rho_avg = rho_avg/tt
u_avg = u_avg/tt
v_avg = v_avg/tt
w_avg = w_avg/tt

D = dict()
D['zeta'] = zeta_avg
D['salt'] = salt_avg
D['u'] = u_avg
D['v'] = v_avg
D['w'] = w_avg
D['rho'] = rho_avg

out_fn = 'aestus2_unstrat_ae1_averages.p'
pickle.dump( D, open(out_fn,'wb'))

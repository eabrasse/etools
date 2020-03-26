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
    d_list = os.listdir(dir0 + f_dir)
    d_list.sort()
    d_list = [x for x in d_list if x[:9]=='ocean_dia']
    #
    for di in d_list:
        fn = dir0 + f_dir + '/' + di
        ds = nc.Dataset(fn)
        if tt == 0:
            
            u_prsgrd = ds['u_prsgrd'][0,:,:,:]
            NZ,NY,NX = np.shape(u_prsgrd)

            u_prsgrd_avg = np.zeros((NZ, NY, NX))
            v_prsgrd_avg = np.zeros((NZ, NY-1, NX+1))
            ubar_prsgrd_avg = np.zeros((NY, NX))
            vbar_prsgrd_avg = np.zeros((NY-1, NX+1))
                    
        u_prsgrd = ds['u_prsgrd'][0, :, :, :].squeeze()
        v_prsgrd = ds['v_prsgrd'][0, :, :, :].squeeze()
        ubar_prsgrd = ds['ubar_prsgrd'][0, :, :].squeeze()
        vbar_prsgrd = ds['vbar_prsgrd'][0, :, :].squeeze()
        
        u_prsgrd_avg = u_prsgrd_avg + u_prsgrd
        v_prsgrd_avg = v_prsgrd_avg + v_prsgrd
        ubar_prsgrd_avg = ubar_prsgrd_avg + ubar_prsgrd
        vbar_prsgrd_avg = vbar_prsgrd_avg + vbar_prsgrd
    
        ds.close()
        tt += 1

u_prsgrd_avg = u_prsgrd_avg/tt
v_prsgrd_avg = v_prsgrd_avg/tt
ubar_prsgrd_avg = ubar_prsgrd_avg/tt
vbar_prsgrd_avg = vbar_prsgrd_avg/tt

D = dict()
D['u_prsgrd'] = u_prsgrd_avg
D['v_prsgrd'] = v_prsgrd_avg
D['ubar_prsgrd'] = ubar_prsgrd_avg
D['vbar_prsgrd'] = vbar_prsgrd_avg

out_fn = 'aestus1_base_ae1_avg_prsgrd.p'
pickle.dump( D, open(out_fn,'wb'))

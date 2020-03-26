"""
Script to call TEFfun to calculate over 20 sections
"""

import os
import sys
import netCDF4 as nc
import pickle
import numpy as np

pth = os.path.abspath('/pmr4/eab32/LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import zfun
import zrfun
import Lfun
import TEFfun_dy as TEFfun#_test as TEFfun

out_dir = '/pmr4/eab32/etools/TEFdy/'
Lfun.make_dir(out_dir)

T0 = 60 # start day
T1 = 65 # end day

nboxes =  2# how many boxes there will be
nsects = nboxes * 2 # no. of estuary sections = no. of boxes, so total no. of sects = twice the no. of boxes

estlons = np.linspace(0.02,0.62,nboxes)
shelflons = np.linspace(-0.59,-0.19,nboxes)
shelftops = np.linspace(45.19,45.59,nboxes)
shelftops = shelftops[::-1]
shelfbots = np.linspace(44.19,44.59)


for i in range(nboxes):
    print('Working on estuary line %.2f' % estlons[i])
    D = TEFfun.TEF_line(estlons[i], [T0, T1])
    D['lon'] = estlons[i]
    out_fn = out_dir + 'TEFline_%i.p' % i
    pickle.dump(D, open(out_fn, 'wb'))
    
    print('Working on shelf box at %.2f' % shelflons[i])
    D = TEFfun.TEF_box([shelflons[i],1.5],[shelfbots[i],shelftops[i]], [T0, T1])
    D['lon'] = shelflons[i]
    D['lats'] = [shelfbots[i], shelftops[i]]
    out_fn = outdir + 'TEFbox_%i.p' % i
    pickle.dump(D,open(out_fn,'wb'))

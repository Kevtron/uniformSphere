#!/usr/bin/env python
import glob
import pickle
import numpy as np

b=glob.glob('p4000_g+3.0_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00/*tau*kappa*.dat')
opa=np.zeros(56)

for i in range(56):
    inf=b[i]
    f=open(inf,'rb')
    kappa=b[i].split('_')[13].split('=')[1].split('.dat')[0]
    array=np.array(pickle.load(f))
    opa[i]=np.reshape(array,(1071,2))[403][0]
    opa[i]*=float(kappa)
f=open('p4000_g+3.0_m0.0_t00_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00/500kappa.dat','wb')
pickle.dump(opa,f)

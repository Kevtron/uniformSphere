#!/usr/bin/env python
import glob
import pickle
import numpy as np

b=glob.glob('sample_data/opacities/*tau*kappa*.dat')
opa=np.zeros(56)

for i in range(56):
    inf=b[i]
    f=open(inf,'rb')
    kappa=b[i].split('_')[13].split('=')[1].split('.dat')[0]
    array=np.array(pickle.load(f))
    opa[i]=np.reshape(array,(1071,2))[403][0]
    opa[i]*=float(kappa)
f=open('sample_data/opacities/500kappa.dat','wb')
pickle.dump(opa,f)

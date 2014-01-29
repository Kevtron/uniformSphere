#!usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

flux = np.genfromtxt('flux.dat')
flux*=2.99e18
mag = 2.5*np.log(flux/8.4e-10)
time=np.arange(0,100.5,0.5)
time*=1.5931e3
time/=86400
print len(time), len(flux)
plt.plot(time,mag)
plt.show()


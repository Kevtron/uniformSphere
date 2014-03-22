#!/usr/bin/env python
#This works because it only takes the rows with eight columns.
import numpy as np
from glob import glob

inf=glob('p*.opa')
for item in inf:
	T=item.split('_')[0].split('p')[1]
	g=item.split('_')[1].split('+')[1]
	print 'p%s-g%s-control.dat'%(T,g)
	inp=np.genfromtxt(item,filling_values=None,skip_header=110,skip_footer=11,invalid_raise=False)
	np.savetxt('p%s-g%s-control.dat'%(T,g),inp,fmt='%.4e')

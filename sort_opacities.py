#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import os
import subprocess
import pickle
work=[[0,0]]*1071
op=glob.glob('sample_data/*.opa')
for f in op:
	wav=np.genfromtxt(f,skip_header=2,skip_footer=20059)
	wavout=np.append(wav,np.genfromtxt(f,skip_header=109,skip_footer=20058))
	wavf=open('sample_data/opacities'+'/'+'wavelengths.dat','wb')
	pickle.dump(wavout,wavf)
	wavf.close()
	inf=open(f,'r').readlines()[110:-10] #Skip the header and the footer
	for i in range(56):
		out=inf[i*358:((i+1)*(358))] #Grab one depth point
		opacity=[]
		control=[]
		for j in range(len(out)):
			if j == 0:
				control.append(out[0].split()[1])
			else: #Opacities
				opacity.append(out[j].split()) #split at spaces
			#	print j,len(opacity)
		for j in range(len(opacity)):
				for k in range(6):
					opacity[j][k]=float(opacity[j][k])
			#		print opacity[j][k]
    		outfile=('sample_data/opacities'+'/'+'%02d'+'_tau='+out[0].split()[1]+'_kappa='+out[0].split()[7]+'.dat') % i
		outf=open(outfile,'wb')
		pickle.dump(opacity,outf)	
		outf.close()

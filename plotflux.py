#!/usr/bin/env python
"""
Lambda = 500 nm
T=5600 K
rho = 4.405e-9 g/cm^3
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import os,sys,glob
import pickle
def IntegrateFlux(tau):
	dmu=0.001
	tmp = np.zeros(2/dmu)
	for x in range(0,1001):
		mu=x/1000.0+0.00000001
		tmp[x]=(1-np.exp(-tau*mu))*mu*dmu
	return sum(tmp)

def makeImage(kappa, filename = None, imagefile = None):
	if filename is None:
		filename = 'image.out'
	f=open(filename)
	lines=f.readlines()
	iformat = int(lines[0])
	nx,ny=lines[1].split()
	nf=int(lines[2])
	sizepix_x,sizepix_y = lines[3].split()
	wavelength=lines[4]
	f.close()
	temp = np.genfromtxt(filename,skip_header = 6)
	image = np.reshape(temp, (int(nx),int(ny)), order="F")
	#image/=1.8763334803e-05
	#image/=1.8763334803e-05
	#image/=1.87808550949e-05#This is just to normalize to the brightest pixel
	flux = np.sum(image)	
	flux = flux*float(sizepix_x)*float(sizepix_y)
	pc=3.0857200e+18#cm
	flux=flux/(pc**2)
	if imagefile is None:
		imagefile = "image.png"
	if not os.path.exists(imagefile):
		plt.clf()
		xi = (np.arange(int(nx))-int(nx)/2+0.5)*float(sizepix_x)
		yi = (np.arange(int(ny))-int(ny)/2+0.5)*float(sizepix_y)
		#make the figure
		print max(image.flatten())#flatten
		fig = plt.figure()
		ax = fig.add_subplot(111,aspect = 'equal')
		cax = ax.imshow(image, cmap=cm.gray, interpolation='nearest', vmin=0, vmax=1, extent=[-float(sizepix_x)*int(nx)/2,float(sizepix_x)*int(nx)/2,-float(sizepix_y)*int(ny)/2,float(sizepix_y)*int(ny)/2])
		cbar = fig.colorbar(cax)
		cbar.set_label('Normalized Intensity')
		tau = kappa*2*1.45042516657e10*6.4e-10
		plt.title(r'Central $\tau$ = %e' % tau)
		plt.xlabel('x (cm)')
		plt.ylabel('y (cm)')	
		fig.savefig(imagefile)
		plt.clf()
	return flux

def planck(T,wav):
	kb=1.3806488e-23
	c=299792458
	h= 6.62606957e-34
	nu=c/wav
	return 2*h*np.power(nu,3)/(np.power(c,2)*(np.exp(h*nu/(kb*T))-1.0))

def numericalCurve(kappa):
	Flux=np.zeros(len(kappa))
	tau=np.zeros(len(kappa))
	for i in range(len(kappa)):
		tau[i]=0.000000001*1.5e10*2.0*kappa[i]# kappa rho d 
		Flux[i]=2*IntegrateFlux(tau[i])
	#fig=plt.figure()
	#plt.semilogx()
	#plt.ylabel(r'$\frac{F_\nu}{\pi B_\nu}$')
	#plt.xlabel(r'$\tau_0$')
	#plt.title('Ex 2.2')
	#plt.axis([0.01,100,0,1.2])
	#plt.axhline(y=1,xmin=0.001,xmax=100,color='k',ls='--')
	#plt.text(0.1,1.0,r'$F_\nu=\pi B_\nu$')
	#plt.semilogx(kappa,Flux)
	return Flux

def plotCurve(kappa):
	#x=np.genfromtxt('flux.dat')
	e=np.genfromtxt('flux64.dat')
	f=np.genfromtxt('flux128.dat')
	g=np.genfromtxt('flux256.dat')
	x,y,z,m=np.genfromtxt("SING001.ascii",usecols=(0,1,2,8),skip_header=12).T
	M=np.sum(m)
        com_x = np.sum(x*m)/M
        com_y = np.sum(y*m)/M
        com_z = np.sum(z*m)/M

        x_dist = np.abs(x-com_x)
        y_dist = np.abs(y-com_y)
        z_dist = np.abs(z-com_z)

        cRsun=6.959e10
        r=np.sqrt(np.power(x_dist,2)+np.power(y_dist,2)+np.power(z_dist,2))
	r*=cRsun
	e/=np.pi*planck(5600,5e-7)*(max(r)/3.08567758e18)**2*(10**7)/(100**2)#Less spherical, different effective area.
	#x/=np.pi*planck(4000,5e-7)*(max(r)/3.08567758e18)**2*(10**7)/(100**2)#Less spherical, different effective area.
	f/=np.pi*planck(5600,5e-7)*(max(r)/3.08567758e18)**2*(10**7)/(100**2) #normalizing by B_vu
	g/=np.pi*planck(5600,5e-7)*(max(r)/3.08567758e18)**2*(10**7)/(100**2)	
	flux=numericalCurve(kappa)
	fig = plt.figure(figsize=(8, 6)) 
	gs = gridspec.GridSpec(3, 3) 
	ax0 = plt.subplot(gs[:,:])
	p1, = ax0.semilogx(kappa,e,'b.')
	p2, = ax0.semilogx(kappa,f,'k.')
	p3, = ax0.semilogx(kappa,g,'r.')
	p4, = ax0.semilogx(kappa,flux, 'k-')
	plt.legend([p1, p2, p3, p4], [r"$64^3$", r"$128^3$", r"$256^3$",r"Eqn 3.7"], loc=0)
	plt.ylabel(r'Normalized Flux (erg cm$^{-2}$ Hz$^{-1}$ s$^{-1}$)')
	ax0.axhline(y=1,xmin=0.01,xmax=10,color='k',ls='-.')
	ax0.text(0.03,1.05,r'$F_\nu=\pi B_\nu$')	
	plt.xlabel(r'$\kappa$ (cm$^2$ g$^{-1}$)')
	#ax1 = plt.subplot(gs[1])
	#res=(e-flux)/flux
	#ax1.semilogx(kappa, res,'b.')
	#resf=(f-flux)/flux
	#print resf
	#ax1.semilogx(kappa, resf,'k.')
	#resg=(g-flux)/flux
	#print resg
	#ax1.semilogx(kappa, resg,'r.')
	#plt.ylabel('Residual (%)')
	#plt.tight_layout()
	plt.savefig('flux-norm.png')
	plt.clf()	
	fig = plt.figure(figsize=(8, 6)) 
	gs = gridspec.GridSpec(3, 3) 
	ax0 = plt.subplot(gs[:,:])
	p1,=ax0.semilogx(kappa,g/max(g),'r.')
 	p2,=ax0.semilogx(kappa,flux, 'k-')
	plt.legend([p1, p2], [r"RADMC3D",r"Eqn 3.7"], loc=0)
	plt.ylabel(r'Normalized Flux (erg cm$^{-2}$ Hz$^{-1}$ s$^{-1}$)')
	plt.xlabel(r'$\kappa$ (cm$^2$ g$^{-1}$)')
	ax0.axhline(y=1,xmin=0.01,xmax=10,color='k',ls='-.')
	ax0.text(0.03,1.05,r'$F_\nu=\pi B_\nu$')	
	plt.savefig('shape.png')
	

if __name__=='__main__':
	a=glob.glob('image/image_*.out')
	b=glob.glob('*tau*.dat')
	d=np.zeros(len(a))
	temp = np.genfromtxt('image/image_055.out',skip_header = 6)
	j=0
	for i in range(len(temp)):
		if (temp[i] > 1e-80) & (temp[i] != 0.0):
	#	print temp[i]
			j+=1
	print j*139949703.21655273*139949703.21655273#Pixel sizes
	opa=pickle.load(open('500kappa.dat','rb'))	

#	for i in range(len(a)):
#		c=a[i].split('.')[0]+'.png'
#		c='png/'+c.split('/')[1]
#		d[i]=float(b[i].split('_')[2].split('.dat')[0].split('=')[1])
#		makeImage(opa[i],a[i],c)
	plotCurve(opa)

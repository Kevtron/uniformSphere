#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import os
import subprocess
from optparse import OptionParser
import pickle

def readDensity(densfile):	
#	n = np.genfromtxt(densfile, skip_header = 18, skip_footer = 15624, invalid_raise = False)
	test = np.genfromtxt(densfile, skip_header = 19)
	nx = 256#128#64#256#128#int(n[0])
	ny = 256#128#64#256#128#int(n[1])
	nz = 256#128#64#256#128#int(n[2])
	temp = np.zeros((nx,ny*nz))
	rhod = np.zeros((nx,ny,nz))
	print nx, ny*nz
	#widgets = ['Read Density: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
         #  ' ', ETA(), ' ', FileTransferSpeed()]
	#pbar = ProgressBar(widgets=widgets, maxval=ny*nz).start()
	for j in range(ny*nz):
		for i in range(nx):
			print i,j
			temp[i,j]=test[j][i]#*4.405e-9
			if temp[i,j] >= 10e-14:
				temp[i,j] = 1.0e-9 
			else:
				temp[i,j] = 0.0
			#pbar.update(j)
	#pbar.finish()
	rhod=np.reshape(temp, (nx,ny,nz), order="F")
	return rhod,nx,ny,nz

def writeDensities(density,nx,ny,nz):
	dust=open('dust_density.inp','w')
	#stars=open('stellarsrc_density.inp','w')
	dust.write('1\n')
	dust.write('%d\n' % (nx*ny*nz))
	dust.write('1\n')
	#stars.write('1\n')
	#stars.write('%d\n' % (nx*ny*nz))
	#stars.write('1\n')
#	widgets = ['Write Density: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
#           ' ', ETA(), ' ', FileTransferSpeed()]
#	pbar = ProgressBar(widgets=widgets, maxval=nz).start()
	for iz in range(nz):
		for iy in range(ny):
			for ix in range(nx):
				dust.write('%3.9f\n' % density[ix,iy,iz])
	#			stars.write('%3.9f\n' % density[ix,iy,iz])
#				pbar.update(iz)
#	pbar.finish()
	dust.close()
#	stars.close()
	return

def makeAMRgrid(nx,ny,nz,sizex,sizey,sizez):
	xi = [-sizex + 2*sizex*(i)/(nx) for i in range(nx+1)] 
	yi = [-sizey + 2*sizey*(i)/(ny) for i in range(ny+1)]
	zi = [-sizey + 2*sizez*(i)/(nz) for i in range(nz+1)]
	amr=open('amr_grid.inp','w')
	amr.write('1\n') 					#iformat
	amr.write('0\n')					#AMR grid style
	amr.write('0\n')					#coordinate system
	amr.write('0\n')					#grid info
	amr.write('1\t1\t1\n')				#include x,y,z
	amr.write('%i\t%i\t%i\n' % (nx,ny,nz))				#size of grid
	[amr.write('%3.9f\n' % xi[i]) for i in range(nx+1)]
	[amr.write('%3.9f\n' % yi[i]) for i in range(ny+1)]
	[amr.write('%3.9f\n' % zi[i]) for i in range(nz+1)]
	amr.close()
	return

def writeControlFiles(nphot):
	radmc3d=open('radmc3d.inp','w')
	dustopac=open('dustopac.inp','w')
	dustopac.write("""
2               Format number of this file
1               Nr of dust species
====================================================================
1               Way in which this dust species is read
0               0=Thermal grain
starstuff       Extension of name of dustkappa_***.inp file
-------------------------------------------------------------------
""")
	dustopac.close()
	radmc3d.write('nphot = %i\n' % nphot)
	radmc3d.write('scattering_mode_max = 0\n')
	radmc3d.write('iranfreqmode = 1\n')
	radmc3d.write('nphot_spec = 100000\n')
	radmc3d.write('incl_dust = 1\n')
	return

def readTemperature(tempfile,MH = 1.67372e-24,KB = 1.38065e-16):
	#skip footer = ny*nz
#	n = np.genfromtxt(tempfile, skip_header = 18, skip_footer = 147456, invalid_raise = False)
	test = np.genfromtxt(tempfile, skip_header = 19)
	nx = 256#128#64#256#128#int(n[0])
	ny = 256#128#64#256#128#int(n[1])
	nz = 256#128#64#256#128#int(n[2])
	array = np.zeros((nx,ny*nz))
	temp = np.zeros((nx,ny,nz))
	for j in range(ny*nz):
		for i in range(nx):
			array[i,j]=test[j][i]#*(2.0*MH)/(3.0*KB)*1.90e15
			array[i,j]=4000
			#	if array[i,j] < 0.25:
	#			array[i,j]=0.0
	
	temp=np.reshape(array, (nx,ny,nz), order="F")
	return temp,nx,ny,nz

def writeTemperature(temperature,nx,ny,nz):
	temp=open('dust_temperature.dat','w')
	temp.write('1\n')
	temp.write('%d\n' % (nx*ny*nz))
	temp.write('1\n')
#	widgets = ['Write Temperature: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
#          ' ', ETA(), ' ', FileTransferSpeed()]
#	tbar = ProgressBar(widgets=widgets, maxval=nz).start()
	for iz in range(nz):
		for iy in range(ny):
			for ix in range(nx):
				temp.write('%3.9f\n' % temperature[ix,iy,iz])
#				tbar.update(iz)
#	tbar.finish()
	temp.close()
	return
	
def makeWavelength():	
	lambda1=0.01
	lambda2=7.0
	lambda3=25.0
	lambda4=1e4
	n12=1000
	n23=100
	n34=30
	lam12 = [lambda1*(lambda2/lambda1)**(i/(1.0*n12)) for i in range(n12)]
	lam23 = [lambda2*(lambda3/lambda2)**(i/(1.0*n23)) for i in range(n23)]
	lam34 = [lambda3*(lambda4/lambda3)**(i/(1.0*n34)) for i in range(n34)]
	[lam12.append(lam23[i]) for i in range(n23)]
	[lam12.append(lam34[i]) for i in range(n34)]
	#
	# Print wavelength.pinp
	#
	f=open('wavelength_micron.inp','w')
	f.write('%d\n' % len(lam12))
#	widgets = ['Write Wavelength: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
#          ' ', ETA(), ' ', FileTransferSpeed()]
#	wbar = ProgressBar(widgets=widgets, maxval=len(lam12)).start()
	for i in range(len(lam12)):
		f.write('%f\n' % lam12[i])
#		wbar.update(i)
#	wbar.finish()
	f.close()
	return lam12,len(lam12)

def stellarTemplate(lam,nlam,nx,ny,nz,tstar,RS,MS):
	template=open('stellarsrc_templates.inp','w')
	template.write('2\n')
	template.write('1\n') 
	template.write('%d\n' % nlam)
	[template.write('%f\n' % lam[i]) for i in range(nlam)]
	template.write('-%d\n' % tstar)
	template.write('%f\n' % RS)
	template.write('%f\n' % MS)
	template.close()
	return
	
def dustOpacity(filename):
	
	opacity = np.genfromtxt(filename, skip_header = 13715, skip_footer = 6096) #opacities
	opacity = flatten(opacity)
	cont = np.genfromtxt(filename, skip_header = 13714, skip_footer = 6453) #control line
	a = np.genfromtxt(filename, skip_header = 2, skip_footer = 20059) #wavelengths
	b = np.genfromtxt(filename, skip_header = 109, skip_footer = 20058) #final wavelength, on its own line
	wav = np.append(a,b)
	ofile = open('dustkappa_starstuff.inp','w')	
	totalopacity = cont[7]
	print "optical depth = %f " % cont[1]
	print "opacity at 500 nm = %f" % totalopacity
	ofile.write('2\n')
	ofile.write('%d\n' % len(wav))
	[ofile.write('%3.9f %3.9f %3.9f\n' % (wav[i]/10000, opacity[2*i]*totalopacity, opacity[2*i+1]*totalopacity)) for i in range(len(wav))]
	ofile.close()
	return

def dustPreprocessed(filename,wavfile):
	inf = open(filename,'rb')
	wavin = open(wavfile,'rb')
	opacity = []
	wav = []
	opacity = pickle.load(inf)	
	wav = pickle.load(wavin)
	opacity = sum(opacity,[])
	ofile = open('dustkappa_starstuff.inp','w')
	tau=float(filename.split('_')[1].split('=')[1])
	kappa=float(filename.split('_')[2].split('=')[1].split('.dat')[0])
	ofile.write('2\n')
	ofile.write('%d\n' % len(wav))
#	print len(opacity),len(wav)
	[ofile.write('%3.9f %3.9f %3.9f\n' % (wav[i]/10000, opacity[2*i]*kappa, opacity[2*i+1]*kappa)) for i in range(len(wav))]
	ofile.close()

def doModel(densfile,tempfile,opacityfile,wavfile):
	AU = 1.49598e13 # Astronomical Unit 	[cm]
	pc = 3.08572e18 # Parsec		[cm]
	MS = 1.98892e33 # Solar Mass		[g]
	TS = 5.78e3     # Solar Temp		[K]
	LS = 3.8525e33	# Solar Luminosity	[erg/s]
	RS = 6.96e10	# Solar radius		[cm]
	MH = 1.67372e-24# Hydrogren Atom 	[g]
	KB = 1.38065e-16# Boltzmann constant 	[erg/k]
	#
	# Monte Carlo parameters
	#
	nphot = 1000000
	rho,nx,ny,nz=readDensity(densfile)
	print 'Density write'
	writeDensities(rho,nx,ny,nz)
	print 'Densities Done'
	makeAMRgrid(nx,ny,nz,2e10,2e10,2e10)
	writeControlFiles(nphot)
	print 'Control and Grid Done'
	temp,nx,ny,nz=readTemperature(tempfile,MH,KB)
	writeTemperature(temp,nx,ny,nz)
	print 'Temp Done'
	wavelength,nlam=makeWavelength()
	print 'Starting Dust'
	dustPreprocessed(opacityfile,wavfile)
	print 'Done'
	return

def doSet(densfiles,ufiles,opacity,wave):
	if not os.path.exists("./image/"):
		os.mkdir("./image/")	
#	widgets = ['Running entire set: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
#          ' ', ETA(), ' ', FileTransferSpeed()]
#	sbar = ProgressBar(widgets=widgets, maxval=len(densfiles)).start()
	for i in range(len(opacity)):
		doModel(densfile[0],ufile[0],opacity[i],wave[0])
		subprocess.call(['radmc3d', 'image', 'lambda', '0.5', 'npix', '500', 'incl', '90'])
		os.rename('image.out','./image/image_%03d.out' % i)
#		sbar.update(i)
#	sbar.finish()	
	return

def makeImage(filename = None,imagefile = None):
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
		fig = plt.figure()
		ax = fig.add_subplot(111,aspect = 'equal')
		cax = ax.imshow(image,cmap=plt.cm.bone, interpolation='nearest',extent=[-float(sizepix_x)*int(nx)/2,float(sizepix_x)*int(nx)/2,-float(sizepix_y)*int(ny)/2,float(sizepix_y)*int(ny)/2])
		ax.contour(image,cmap=plt.cm.bone)
		cbar = fig.colorbar(cax)
		plt.title(r'fliux at 1 pc = %s erg/cm^2/s/Hz' % flux)
		fig.savefig(imagefile)
	return flux

def lightCurve(substring = None):
	files = glob.glob("./image/image_*.out")
	flux = np.array([])
#	widgets = ['Making lightcurve: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
 #          ' ', ETA(), ' ', FileTransferSpeed()]
#	lbar = ProgressBar(widgets=widgets, maxval=len(files)).start()
	if not os.path.exists("./png/"):
		os.mkdir("./png/")
	for i in range(len(files)):
		a=makeImage(files[i], "./png/image_%03d.png" % i)
		flux = np.append(flux,a)
		print a
#		lbar.update(i)
#	ebar.finish()
	fig = plt.figure()
	#plt.ylim([0,0.5e-35]) 
	ax=fig.add_subplot(111,aspect = 'equal')
	#flux*=2.99e18
	#mag = 2.5*np.log(flux/8.4e-10)
	#time=np.arange(0,100.5,0.5)
	#time*=1.5931e3
	#time/=86400
	#ax.plot(time, mag)
	ax.plot(flux)
	plt.savefig("flux.png")
	f=open('flux.dat','w')
	for i in range(len(flux)):
		f.write("%e\n" % (flux[i]))#,time[i]))
	f.close()

if __name__=="__main__":
	densfile = glob.glob('SING*density*')
	ufile = glob.glob('SING*u*')
	opacity = glob.glob('opacities/*tau*')
	wav = glob.glob('wavelengths.dat')
	#densfile[0],ufile[0],opacity[0],wav[0]
	print densfile[0]
	doModel(densfile[0],ufile[0],opacity[0],wav[0])
	#makeImage()
	doSet(densfile,ufile,opacity,wav)
	lightCurve()

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
	test = np.genfromtxt(densfile, skip_header = 19)
	nx = 256#hard coded, could probably stand to read from file.
	ny = 256
	nz = 256
	temp = np.zeros((nx,ny*nz))
	rhod = np.zeros((nx,ny,nz))
	for j in range(ny*nz): #Set the density equal to 10e-9 inside the sphere, 0.0 outside the sphere.
		for i in range(nx):
			print i,j
			temp[i,j]=test[j][i]
			if temp[i,j] >= 10e-14:
				temp[i,j] = 1.0e-9 
			else:
				temp[i,j] = 0.0
	rhod=np.reshape(temp, (nx,ny,nz), order="F")
	return rhod,nx,ny,nz

def writeDensities(density,nx,ny,nz):
	dust=open('dust_density.inp','w')
	dust.write('1\n')
	dust.write('%d\n' % (nx*ny*nz))
	dust.write('1\n')
	for iz in range(nz):
		for iy in range(ny):
			for ix in range(nx):
				dust.write('%3.9f\n' % density[ix,iy,iz]) #Write the radmc3d dust_density.inp file
	dust.close()
	return

def makeAMRgrid(nx,ny,nz,sizex,sizey,sizez): #Make an AMR grid for radmc3d
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

def writeControlFiles(nphot): #Write the required control file for radmc3d
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
	test = np.genfromtxt(tempfile, skip_header = 19)
	nx = 256
	ny = 256
	nz = 256
	array = np.zeros((nx,ny*nz))
	temp = np.zeros((nx,ny,nz))
	for j in range(ny*nz):
		for i in range(nx):
			array[i,j]=4000 #set the temperature everywhere to 4000K. A bit lazy, but because the density is strictly zero outside the sphere, it doesn't matter for this test.
	temp=np.reshape(array, (nx,ny,nz), order="F")
	return temp,nx,ny,nz

def writeTemperature(temperature,nx,ny,nz):
	temp=open('dust_temperature.dat','w')
	temp.write('1\n')
	temp.write('%d\n' % (nx*ny*nz))
	temp.write('1\n')
	for iz in range(nz):
		for iy in range(ny):
			for ix in range(nx):
				temp.write('%3.9f\n' % temperature[ix,iy,iz])#write the "dust_temperature.dat" file for radmc3d.
	temp.close()
	return
	
def makeWavelength():	#Make the wavelength file. wavelength range is mostly arbitrary.
	lambda1=0.01
	lambda2=1e4
	n12=1000
	lam12 = [lambda1*(lambda2/lambda1)**(i/(1.0*n12)) for i in range(n12)]
	f=open('wavelength_micron.inp','w')
	f.write('%d\n' % len(lam12))
	for i in range(len(lam12)):
		f.write('%f\n' % lam12[i])
	f.close()
	return lam12,len(lam12)

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
	[ofile.write('%3.9f %3.9f %3.9f\n' % (wav[i]/10000, opacity[2*i]*kappa, opacity[2*i+1]*kappa)) for i in range(len(wav))]
	ofile.close()

def doModel(densfile,tempfile,opacityfile,wavfile):
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
	for i in range(len(opacity)):
		doModel(densfile[0],ufile[0],opacity[i],wave[0])
		subprocess.call(['radmc3d', 'image', 'lambda', '0.5', 'npix', '500', 'incl', '90'])
		os.rename('image.out','./image/image_%03d.out' % i)
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
	if not os.path.exists("./png/"):
		os.mkdir("./png/")
	for i in range(len(files)):
		a=makeImage(files[i], "./png/image_%03d.png" % i)
		flux = np.append(flux,a)
		print a
	fig = plt.figure()
	ax=fig.add_subplot(111,aspect = 'equal')
	ax.plot(flux)
	plt.savefig("flux.png")
	f=open('flux.dat','w')
	for i in range(len(flux)):
		f.write("%e\n" % (flux[i]))#,time[i]))
	f.close()

if __name__=="__main__":
	#Automatically runs through the entire set of opacities, calculating images and fluxes for use in comparing to theoretical results.
	densfile = glob.glob('SING*density*')
	ufile = glob.glob('SING*u*')
	opacity = glob.glob('opacities/*tau*')
	wav = glob.glob('wavelengths.dat')
	doModel(densfile[0],ufile[0],opacity[0],wav[0])
	doSet(densfile,ufile,opacity,wav)
	lightCurve()

#Uniform Sphere Project#
==========

##Introduction##

I wrote this code to set up a series of initial conditions as a test case for my Master's Thesis project. The code sets up a uniform density sphere in a cubic grid and then applies a wavelength dependent global opacity from a [stellar atmosphere model](http://marcs.astro.uu.se/index.php "MARCS"). The sphere is then imaged using the [RADMC3D](http://www.ita.uni-heidelberg.de/~dullemond/software/radmc-3d/ "RADMC3D") code. The atmosphere model and uniform density are chosen such that the sphere becomes optically thick and eventually attains uniform brightness as the opacity increases. This case is designed to test that the code works as expected when given stellar atmosphere opacities instead of the usual dust opacity.

#sort_opacities.py#

This code reads a [.opa file](http://marcs.astro.uu.se/documents.php?doc=opac "MARCS Opacity files") and separates out each individual depth point, writing the wavelength dependent opacity to file. This allows for easy input to the RADMC3D code when increasing the opacity.

#500kappa.py#

Because the MARCS opacities are normalized by the opacity at 500 nm, and because this opacity is recorded in the control line before each block of ~1000 pairs, it was easier to just grab each control opacity and write it to file, to account for the normalization.

#uniformSphere.py#

This code takes density and internal energy fields of an SPH simulation which have been interpolated to a regular cartesian grid using the [SPLASH](http://users.monash.edu.au/~dprice/splash/ "SPLASH") tool and rewrites them in the format expected by RADMC3D as dust_density.inp and dust_temperature.inp. It also reads one file generated by sort_opacities.py and rewrites it as "dustkappa_starstuff". It also produces other generic input files required by RADMC3D. After this it calls RADMC3D to produce a synthetic image of the sphere and then loops over all existing opacity files. 

#plotflux.py#

This code reads all the image*.out files produced by radmc3d and plots a curve of brightness as a function of global opacity. 

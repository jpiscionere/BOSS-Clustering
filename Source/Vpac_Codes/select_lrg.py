#!/usr/bin/python
#
#  extract the table information from a fits file and print out as a space-delimited ASCII table
#
# Martin Paegert and Joshua Pepper 
# Last modified: 2011-05-13
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

import pyfits
import numpy as np 
## Read in the fits file
hdulist = pyfits.open('object_sdss_imaging.fits')


## Uncomment the following line to view the info about the table
print hdulist.info()

## Read the tabular portion of the fits file into the variable 'table'.  This assumes that the table of interest is located in extension 1 
table = hdulist[1].data
#mask= 17.5 < table.field('CMODELMAG')[3] <19.9  
#newtbdata=table[mask]
#hdu=pyfits.BinTableHDU(newtbdata)



ra=table.field('RA')
dec=table.field('DEC')
petro_r=table.field('PETROR50')
petro_r_r=petro_r[:,2]


fiberflux=table.field('FIBERFLUX')

fiberflux_g=fiberflux[:,1]
fiberflux_r=fiberflux[:,2]
fiberflux_i=fiberflux[:,3]


modelflux=table.field('MODELFLUX')
petroflux=table.field('PETROFLUX')
psfflux=table.field('PSFFLUX')


psf_r=psfflux[:,2]
modelflux_r=modelflux[:,2]
petroflux_r=petroflux[:,2]

surface_brightness=petroflux_r + 2.5*np.log10(2*np.pi*petro_r_r**2)
c_perp=(fiberflux_r-fiberflux_i) - (fiberflux_g-fiberflux_r)/4.0 - 0.18
c_perp=np.abs(c_perp)

C_constant=0.7

c_par=C_constant*(fiberflux_g-fiberflux_r) + (1.0-C_constant)*4.0*((fiberflux_r-fiberflux_i)-0.18)








petro_cut_par= 13.1 + c_par/0.3
petro_cut=19.2
c_perp_cut=0.2
mu_cut=24.2
psf_model_cut=0.3


array=np.column_stack((ra,dec,surface_brightness,petroflux_r,modelflux_r,psf_r,c_perp,petro_cut_par))
np.savetxt('lrg_selection.txt',array,delimiter='\t',newline='\n')

ids=np.where((petroflux_r < petro_cut_par) & (petroflux_r < petro_cut) & (c_perp < c_perp_cut) & (surface_brightness < mu_cut) & ( psf_r - modelflux_r < psf_model_cut))
 

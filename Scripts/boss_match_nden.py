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
import astropy
from astropy import cosmology
from astropy.cosmology import get_current
from astropy.cosmology import LambdaCDM

def Calculate_Magnitude( flux_r, extinction, redshift):
	#Calculate apparent magnitude

	m_r=22.5-2.5*np.log10(flux_r)-extinction_r

	#Calculate Mag Correction
	###THERE IS A HUGE CAVEAT THAT THIS IS THE OLD WAY OF CALCULATING LRG MAGS AND WILL NOT WORK FOR Z > 0.6 GALAXIES
	interp_me=np.genfromtxt("/home/piscioja/Clustering/Boss/Source/Vpac_Codes/correction.txt",names=('z_','delta_g','g_r'))

	delta_g=np.interp(redshift,interp_me[('z_')],interp_me[('delta_g')])
	g_r=np.interp(redshift,interp_me[('z_')],interp_me[('g_r')])
	z_cal=0.2

	#Distance Modulus

	DM=cosmo.distmod(redshift)

	#Calculate Absolute Magnitude

	M_g = m_r - DM - delta_g + g_r - z_cal

	return M_g	 


## Read in the fits file
hdulist = pyfits.open('/hd0/Research/Clustering/Boss/dr11/dr11v2/cmass-dr11v2-N-Anderson.dat.fits')


## Uncomment the following line to view the info about the table
#print hdulist.info()

## Read the tabular portion of the fits file into the variable 'table'.  This assumes that the table of interest is located in extension 1 
table = hdulist[1].data

cosmo = LambdaCDM(H0=100, Om0=0.274, Ode0=0.726)

ra=table.field('RA')
dec=table.field('DEC')
redshift=table.field('Z')
fibcol=table.field('WEIGHT_CP')
poly=table.field('IPOLY')
sector=table.field('ISECT')
extinction=table.field('EXTINCTION')
extinction_r=extinction[:,3]
flux=table.field('MODELFLUX')
flux_r=flux[:,3]

M_g=Calculate_Magnitude(flux_r,extinction_r,redshift)

array=np.column_stack((ra,dec,redshift,fibcol,poly,M_g))
array=array[ (array[:,2] > 0.43) & (array[:,2] < 0.55)  ]

nden_target=0.00013
sky_coverage=2.036048 #in rad
percent_sky=sky_coverage / ( 4 * np.pi )


big_distance=cosmo.comoving_distance(np.max(array[:,2]))
small_distance=cosmo.comoving_distance(np.min(array[:,2]))



max_volume=4./3. * np.pi * cosmo.comoving_distance(np.max(array[:,2]))**3
min_volume=4./3. * np.pi * cosmo.comoving_distance(np.min(array[:,2]))**3

volume= percent_sky * (max_volume-min_volume)



number_of_galaxies=nden_target*volume

print number_of_galaxies

#final_array=array[array[:,5].argsort()[::-1][:number_of_galaxies]]
#Magnitudes are Negative Goddammit
final_array=array[array[:,5].argsort()[:][:number_of_galaxies]]


np.savetxt('nden_matched_array.txt',final_array,fmt="%lf %lf %lf %d %d %lf", delimiter='\t',newline='\n') 

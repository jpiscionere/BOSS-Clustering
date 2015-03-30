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
from astropy.cosmology import FlatLambdaCDM
cosmo=FlatLambdaCDM(H0=100,Om0=0.266)

## Read in the fits file
hdulist = pyfits.open('/hd0/Research/Clustering/Boss/dr11/dr11v2/cmass-dr11v2-N-Anderson.dat.fits')


## Uncomment the following line to view the info about the table
print hdulist.info()

## Read the tabular portion of the fits file into the variable 'table'.  This assumes that the table of interest is located in extension 1 
table = hdulist[1].data


ra=table.field('RA')
dec=table.field('DEC')
redshift=table.field('Z')
polygon=table.field('IPOLY')
extinction=table.field('EXTINCTION')

extinction_g=extinction[:,1]
extinction_r=extinction[:,2]
extinction_i=extinction[:,3]

fiberflux=table.field('FIBER2FLUX')
fiberflux_i= 22.5-2.5*np.log10(fiberflux[:,3]) - extinction_i


modelflux=table.field('MODELFLUX')
modelflux_g= 22.5-2.5*np.log10( modelflux[:,1]) - extinction_g
modelflux_r= 22.5-2.5*np.log10( modelflux[:,2]) - extinction_r
modelflux_i= 22.5-2.5*np.log10( modelflux[:,3]) - extinction_i


d_perp=(modelflux_r-modelflux_i) - (modelflux_g-modelflux_r)/8.



i_cmod_cut= 19.86 + 1.6*(d_perp - 0.8)
d_perp_cut=0.55


weight_cp=table.field('WEIGHT_CP')
icollided=table.field('ICOLLIDED')

distance_modulus=cosmo.distmod(redshift)
Mag_i=modelflux_i - distance_modulus - (-0.5)


array=np.column_stack((ra,dec,redshift,weight_cp,polygon,modelflux_g,modelflux_r,modelflux_i,fiberflux_i,distance_modulus,Mag_i))
dimensions=str("ra dec redshift weight_cp polygon modelflux_g modelflux_r modelflux_i fiberflux_i distance_modulus Mag_i")


np.savetxt('/hd0/Research/Clustering/Boss/dr11/dr11v2/dr11v2_all.out',array,delimiter='\t',newline='\n',header=str(dimensions),comments=' ')

ids=np.where(( modelflux_i > 17.5  ) & 
	( modelflux_i < 19.9) & 
	( modelflux_r -  modelflux_i < 2 ) & 
	( d_perp > d_perp_cut) & 
	( fiberflux_i < 21.5 ) &
	( modelflux_i < i_cmod_cut))

array_filter=array[ids]

bin=(0.43,0.5,0.55,0.6,0.7)

np.savetxt('/hd0/Research/Clustering/Boss/dr11/dr11v2/dr11v2_imaging.txt',array_filter,delimiter='\t',newline='\n')



path=('/hd0/Research/Clustering/Boss/dr11/dr11v2/')


for x in range(0,4):
	array_list=array_filter[np.where((array_filter[:,2] > bin[x]) & (array_filter[:,2] < bin[x + 1]))]	
	np.savetxt(str(path) + "bin" + str(x+1) + "_selection.txt",array_list,delimiter='\t',newline='\n')


array_output=array_filter[:,range(0,5)]

np.savetxt('/hd0/Research/Clustering/Boss/dr11/dr11v2/Di.dr11v2.out.selection',array_output,delimiter='\t',newline='\n')


for x in range(0,4):
        array_list=array_output[np.where((array_filter[:,2] > bin[x]) & (array_filter[:,2] < bin[x + 1]))]
        np.savetxt(str(path) + "bin" + str(x+1) + "_Ds_dr11v2_selection.out",array_list,delimiter='\t',newline='\n')

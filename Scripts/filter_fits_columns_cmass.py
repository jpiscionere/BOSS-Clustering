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



array=np.column_stack((ra,dec,redshift,polygon,modelflux_i,modelflux_r))

ids=np.where(( modelflux_i > 17.5  ) & 
	( modelflux_i < 19.9) & 
	( modelflux_r -  modelflux_i < 2 ) & 
	( d_perp > d_perp_cut) & 
	( fiberflux_i < 21.5 ) &
	( modelflux_i < i_cmod_cut))

array_filter=array[ids]

bin=(0.43,0.5,0.55,0.6,0.7)

np.savetxt('/hd0/Research/Clustering/Boss/dr11/dr11v2/dr11v2_imaging.txt',array_filter,delimiter='\t',newline='\n')

array=np.column_stack((ra_,dec_,redshift_,weight_cp_,polygon_))

path=('/hd0/Research/Clustering/Boss/dr11/dr11v2/')


for x in range(0,3):
	array_list=array[np.where((array[:,2] > bin[x]) & (array[:,2] < bin[x + 1]))]	
	np.savetxt(str(path) + "bin" + str(x) + "selection.text",array_list,delimiter='\t',newline='\n')
 		


 

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

modelflux=table.field('MODELFLUX')




modelflux_g= 22.5-2.5*np.log10( modelflux[:,1]) - extinction_g
modelflux_r= 22.5-2.5*np.log10( modelflux[:,2]) - extinction_r
modelflux_i= 22.5-2.5*np.log10( modelflux[:,3]) - extinction_i

weight_cp=table.field('WEIGHT_CP')
icollided=table.field('ICOLLIDED') + 1


array=np.column_stack((ra,dec,redshift,icollided,polygon,modelflux_g,extinction_g,modelflux_r,extinction_r,modelflux_i,extinction_i))
np.savetxt('/hd0/Research/Clustering/Boss/dr11/dr11v2/dr11v2_mag_extinction.txt',array,fmt='%lf',delimiter='\t',newline='\n') 

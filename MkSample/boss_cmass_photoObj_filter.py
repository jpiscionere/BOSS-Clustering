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
hdulist = pyfits.open('cmass-dr11v2-N-Anderson.dat.fits')


## Uncomment the following line to view the info about the table
print hdulist.info()

## Read the tabular portion of the fits file into the variable 'table'.  This assumes that the table of interest is located in extension 1 
table = hdulist[1].data
#mask= 17.5 < table.field('CMODELMAG')[3] <19.9  
#newtbdata=table[mask]
#hdu=pyfits.BinTableHDU(newtbdata)



ra=table.field('RA')
dec=table.field('DEC')
redshift=table.field('Z')
fibcol=table.field('WEIGHT_CP')
poly=table.field('IPOLY')
sector=table.field('ISECT')




array=np.column_stack((ra,dec,redshift,fibcol,poly,sector))

np.savetxt('cmass-dr11v2-N-Anderson.dat.ascii',array,delimiter='\t',newline='\n')







             

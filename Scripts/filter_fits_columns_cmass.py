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


ra=table.field('RA')
dec=table.field('DEC')
redshift=table.field('Z')

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

ids=np.where(( modelflux_i > 17.5  ) & 
	( modelflux_i < 19.9) & 
	( modelflux_r -  modelflux_i < 2 ) & 
	( d_perp > d_perp_cut) & 
	( fiberflux_i < 21.5 ) &
	( modelflux_i < i_cmod_cut))



ra_=ra[ids]
dec_=dec[ids]
redshift_=modelflux_i[ids]
weight_cp_=modelflux_r[ids]
modelflux_i_=modelflux_i[ids]
modelflux_r_=modelflux_r[ids]

array=np.column_stack((ra_,dec_,modelflux_i_,modelflux_r_))
np.savetxt('dr11v2_imaging.txt',array,delimiter='\t',newline='\n')




 

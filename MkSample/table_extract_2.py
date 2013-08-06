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

## Read in the fits file
hdulist = pyfits.open('cmass-dr11v2-N-Anderson.dat.fits')

## Uncomment the following line to view the info about the table
# print hdulist.info()

## Read the tabular portion of the fits file into the variable 'table'.  This assumes that the table of interest is located in extension 1 
table = hdulist[1].data

of = open('cmass-dr11v2-N-Anderson.dat','w')
for line in table:
  for col in line:
    ## To use an alternate delimiter, replace the space inside the quotes with the desired delimiting character/string
    of.write(str(col) + ' ')
  of.write("\n")

of.close()



             

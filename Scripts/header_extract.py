#!/usr/bin/python

import pyfits
import numpy as np

#this assumes that the header you're interested in is not the primary header
hdulist=pyfits.open("cmass-dr11v2-N-Anderson.dat.fits")

#header=pyfits.getheader("cmass-dr11v2-N-Anderson.dat.fits")
prihdr=hdulist[1].header

print prihdr[3]

of = open('cmass-dr11v2-N-Anderson.head','w')
for col in prihdr:
       of.write("{}\n".format(prihdr[col]))

of.close()

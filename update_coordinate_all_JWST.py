
from astropy.io import fits
from astropy import wcs

from astropy.io import fits

import glob
files = glob.glob("*fits*")

i=0
while i<len(files):
    
    with fits.open(files[i], memmap=False) as hdulist:
        w = wcs.WCS(hdulist[1].header)
        hdr = hdulist[1].header
    
    fits.setval(files[i], 'CRVAL1', ext=1, value=hdulist[1].header['CRVAL1']-(-0.61/3600))    
    fits.setval(files[i], 'CRVAL2', ext=1, value=hdulist[1].header['CRVAL2']-(0.12/3600))    
    i=i+1
    


    

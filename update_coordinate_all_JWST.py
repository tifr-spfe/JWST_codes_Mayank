import time
import numpy as np
import scipy


from astropy.io import fits
from astropy import wcs
import astropy.units as u
from astropy.stats import sigma_clip
from astropy.utils.data import download_file
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import regions 
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

from astropy.io import fits

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Ubuntu'
plt.rcParams['font.monospace'] = 'Ubuntu Mono'
plt.rcParams['font.size'] = 35
plt.rcParams['axes.labelsize'] = 45
plt.rcParams['xtick.labelsize'] = 45
plt.rcParams['ytick.labelsize'] =45
plt.rcParams['legend.fontsize'] = 45
plt.rcParams['figure.titlesize'] = 24
plt.rcParams['axes.labelweight']='bold'



plt.rcParams['axes.linewidth'] = 5
plt.rcParams['xtick.major.size'] = 15
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['ytick.major.size'] = 15
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['ytick.major.width'] = 3
plt.rcParams['xtick.minor.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.5

plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True


plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = False
plt.rcParams['xtick.major.pad']='10'
plt.rcParams['ytick.major.pad']='10'

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
    


    

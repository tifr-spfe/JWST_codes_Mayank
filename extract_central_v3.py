import time
import numpy as np
import scipy
from regions import CircleSkyRegion
from astropy.io import fits
from astropy import wcs
import astropy.units as u
from astropy.stats import sigma_clip
from astropy.utils.data import download_file
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from spectral_cube import SpectralCube
import regions 
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
plt.rcParams['font.size'] = 25
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] =25
plt.rcParams['legend.fontsize'] = 25
plt.rcParams['figure.titlesize'] = 12
plt.rcParams['axes.labelweight']='bold'
plt.rcParams['axes.linewidth'] = 3
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['xtick.major.width'] = 5
plt.rcParams['ytick.major.width'] = 3
plt.rcParams['xtick.minor.width'] = 5
plt.rcParams['ytick.minor.width'] = 3
import warnings
warnings.filterwarnings("ignore")
from astropy.coordinates import SkyCoord as sc
wave =[]
flux=[]
fluxerr=[]
bg=[]

sp=[3.6,4.5,5.8,8,24]
sf=[0.59,2.5,4.5,3.1,90]
rr=1.2*u.arcsec
rr1=0.7*u.arcsec


def open_file(filename,cx,cy,rr):
    with fits.open(filename, memmap=False) as hdulist:
        sci = hdulist['SCI'].data
        err = hdulist['ERR'].data
        w = wcs.WCS(hdulist[1].header)
        hdr = hdulist[1].header
    hdulist[2].header=hdulist[1].header
    nw = hdulist[1].shape[0]
    cc1=sc(ra=cx,dec=cy,unit='deg',frame='icrs')
    wave = (np.arange(nw) + hdr["CRPIX3"]) * hdr["CDELT3"] + hdr["CRVAL3"]
    midwave = wave[int(nw / 2)]
    cdelt1 = hdr["CDELT1"]
    cdelt2 = hdr["CDELT2"]
    #h1= (rr1*1.22*(midwave*1e-6/6.5)*(180*3600/np.pi)*u.arcsec)
    px_area = cdelt1 * cdelt2 * 3.0461741978670859934e-4 *1e6
    h1=rr
    print(midwave,h1)
    w,f=extract_flux(hdulist[1],cc1,h1,px_area)
    w,ferr=extract_flux2(hdulist[2],cc1,h1,px_area)
    return w,f,ferr
    #print(hdr)
    
 

    
    
def extract_flux(hdulist,cc1,h1,px_area):
    cube = SpectralCube.read(hdulist)
    region =regions.CircleSkyRegion(center=cc1, radius=h1)    
    subcube = cube.subcube_from_regions([region]) 
    spectrum = subcube.sum(axis=(1, 2)) 
    c_f = np.array(spectrum) 
    c_w = np.array(spectrum.spectral_axis)
    
    return c_w,c_f*px_area



def extract_flux2(hdulist,cc1,h1,px_area):
    cube = SpectralCube.read(hdulist)
    region =regions.CircleSkyRegion(center=cc1, radius=h1)
    subcube = cube.subcube_from_regions([region]) 
    subcube=(subcube**2)
    spectrum = (subcube.sum(axis=(1, 2)))**0.5
    c_f = np.array(spectrum) 
    c_w = np.array(spectrum.spectral_axis)
    
    return c_w,c_f*px_area
    

files = ["jw01802-o015_t012_nirspec_g395m-f290lp_crop1_s3d.fits",
    "Level3_ch1-short_s3d.fits",
    "Level3_ch1-medium_s3d.fits",
    "Level3_ch1-long_s3d.fits",
    "Level3_ch2-short_s3d.fits",
    "Level3_ch2-medium_s3d.fits",
    "Level3_ch2-long_s3d.fits",
    "Level3_ch3-short_s3d.fits",
    "Level3_ch3-medium_s3d.fits",
    "Level3_ch3-long_s3d.fits",
    "Level3_ch4-short_s3d.fits",
    "Level3_ch4-medium_s3d.fits",
    "Level3_ch4-long_s3d.fits",
]



files2 = ["jw01802-o015_t012_nirspec_g395m-f290lp_crop1_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch1-short_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch1-medium_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch1-long_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch2-short_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch2-medium_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch2-long_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch3-short_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch3-medium_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch3-long_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch4-short_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch4-medium_s3d.fits",
    "IRAS16253_Cubes_NoBKG_sub/Level3_ch4-long_s3d.fits",
]


def find_neighboring_differences(arr):
    differences = []
    for i in range(len(arr) - 1):
        diff = arr[i + 1] - arr[i]
        differences.append(diff)
    return differences

#RA=247.0901709 #JWST
#Dec=-24.6065499 #JWST


RA = 166.6932418   
Dec =-77.3757983

RA_bg=166.6908945
Dec_bg=-77.3756768


i=1
while(i<len(files)):
    plt.figure(figsize=(14,10))
    w,f,ferr=open_file(files[i],RA,Dec,rr)
    w,fb,ferrb=open_file(files[i],RA,Dec,rr1)
    f=np.array(f)
    fb=np.array(fb)
    
    #w1,f1,ferr1=open_file(files2[i],RA,Dec)
    wave.append(w)
    flux.append(f)
    fluxerr.append(ferr)

       


    #plt.plot(w1,f1*1000,c='r',label='No BG sub')
    plt.plot(w,f*1000,c='b',label='Target')
    plt.plot(w,fb*1000,c='k',label='BG')
    #plt.plot(w,fb*1000,label='BG ',zorder=100)
    #plt.plot(w,ferr*1000,label='err',zorder=100)
    #w,f,fb=open_file(files2[i],cx[i],cy[i],cxa[i],cy[i],cxb[i],cyb[i],3.33,4.95,flg[i])
    #plt.plot(w,f*1000,label='stage 3 BG sub')
    plt.xlabel('Wavelength ($\mu$m)')
    plt.ylabel('Flux (mJy)')
    plt.legend()
    #plt.scatter(sp,sf,s=100,c='k')
    #plt.yscale('log')
    
    plt.savefig(str(files[i])+'_err.png')
    
    """
    plt.figure(figsize=[12,12])
    wn=find_neighboring_differences(w)
    plt.plot(w[:-1],np.round(wn,5))
    plt.xlabel('Wavelength ($\mu$m)')
    plt.ylabel('$\Delta \lambda$ ($\mu$m)')
    plt.savefig(str(files[i])+'_lamb.png')
    
    plt.figure(figsize=[12,12])
    vel =(wn/w[:-1])*3e5
    plt.plot(w[:-1],vel)
    plt.xlabel('Wavelength ($\mu$m)')
    plt.ylabel('$\Delta V$ ($\mu$m)')
    plt.savefig(str(files[i])+'_vel.png')
    """
    i=i+1

#%%
xxx=[]

def scale(waves,spec1ds,flxerr):
    ii=0
    osubs_left = np.where(waves[ii]>np.min(waves[ii+1]))
    osubs_right = np.where(waves[ii+1]<np.max(waves[ii]))
    if(len(osubs_left)==0 or  len(osubs_right)==0):
        scale=1
    else:
        scale = np.median(spec1ds[ii+1][osubs_right])/np.median(spec1ds[ii][osubs_left])
    spec1ds[ii] *= scale
    flxerr[ii] *= scale
    print(scale)
    xxx.append(scale)
    ii=1
    
    while (ii<len(waves)):
        osubs_left = np.where(waves[ii-1]>np.min(waves[ii]))
        #print(len(osubs_left[0]),np.mean(waves[ii]))
        osubs_right = np.where(waves[ii]<np.max(waves[ii-1]))
        if(len(osubs_left[0])==0 or  len(osubs_right[0])==0):
            scale=scale
        else:
            scale = np.median(spec1ds[ii-1][osubs_left])/np.median(spec1ds[ii][osubs_right])
        spec1ds[ii] *= scale
        flxerr[ii] *= scale
        #spec1ds[ii][spec1ds[ii]<1e-3] =1e-3
        ii+=1
        print(scale)
        xxx.append(scale)
        
    return spec1ds,flxerr 



#%%

import pandas as pd
spec1ds,flxerr =  scale(wave,flux,fluxerr) #flux,fluxerr

from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker
waves_flat = np.concatenate(wave)
spec1d_flat = np.concatenate(spec1ds)
flxerr_flat=np.concatenate(flxerr)
ssubs = np.argsort(waves_flat)
wave_all = waves_flat[ssubs]
flux_all = spec1d_flat[ssubs]
flxerr_all=flxerr_flat[ssubs]

flux_all=flux_all[wave_all>3.1]

flxerr_all=flxerr_all[wave_all>3.1]

wave_all=wave_all[wave_all>3.1]

fig, ax =plt.subplots(figsize=(16,12))







df=pd.DataFrame()

wave_all=np.round(wave_all,decimals=4) 
df['Wavelength'] = wave_all
df['Flux_mJy'] = (flux_all)*1000
df['Flux_err'] = (flxerr_all)*1000


df.to_csv('Knot.csv')



#plt.plot(wave,bg)



plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Flux (mJy)')
plt.plot(wave_all,(flux_all)*1000)
#plt.errorbar(wave_all,(flux_all)*1000,yerr=(flxerr_all)*1000,ecolor='r')

#plt.scatter(sp,sf,s=100,c='k',zorder=110)
plt.xscale('log')
plt.yscale('log')
#plt.legend()

for axis in [ax.xaxis]:
    axis.set_major_formatter(ScalarFormatter())
    axis.set_minor_formatter(ScalarFormatter())

ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))

df1=pd.read_csv('H2_line_list.csv')

y=(1000)

x=df1['wave']
y=np.ones(len(x))*300

x2=[5.3403,17.936,24.519,25.988,4.889,4.115,4.434839]
y2=np.ones(len(x2))*1000

x3=[6.636,6.9852,12.813,25.249]
y3=np.ones(len(x3))*700

plt.scatter(4.05228, 60,marker='|', s=200,c='m',linewidths=3,label='$HI$')

plt.scatter(x, y,marker='|', s=200,c='k',linewidths=3,label='$H_2$')
plt.scatter(x2, y2,marker='|', s=200,c='g',linewidths=3,label='$[FeII]$')
plt.scatter(x3, y3,marker='|', s=200,c='r',linewidths=3,label='FS lines')

plt.text(5.8,11,'H$_2$O bend',size=20)

plt.text(6.7,14,'CH$_3$OH',size=20)

plt.text(8.5,8.5,'Silicates',size=20)

plt.text(15.1,37,'CO$_2$',size=20)

plt.xlim(4.8,29)
plt.ylim(5,2000)

plt.legend()
plt.savefig('IRAS_16253_spec.png')


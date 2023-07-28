#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:19:52 2023

@author: mnarang
"""

import astropy.units as u
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import pvextractor



import numpy as np
import scipy


from astropy.io import fits
from astropy import wcs
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from spectral_cube import SpectralCube

from astropy.io import fits

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


filename = "Level3_ch1-short_s3d.fits"
H1=5.3401687 +(0.00017812912619416128)
ac=0
with fits.open(filename, memmap=False) as hdulist:
    sci = hdulist['SCI'].data
    err = hdulist['ERR'].data
    w = wcs.WCS(hdulist[1].header)
    hdr = hdulist[1].header
    hdulist[2].header=hdulist[1].header
    #print(hdr)

res=2700
    
cube = SpectralCube.read(hdulist[1])
cube_err = SpectralCube.read(hdulist[2])
    
new_mask = (cube !=0* u.MJy / u.sr) 

cube = cube.with_mask(new_mask)  
cube_err=cube_err.with_mask(new_mask)
    


w1=(1-(2/res))*H1
w2=(1+(2/res))*H1
subcube = cube.spectral_slab(w1*u.micron, w2*u.micron)
subcube_err = cube_err.spectral_slab(w1*u.micron, w2*u.micron)

w1=(1-(4/res))*H1
w2=(1-(3/res))*H1
subcube2 = cube.spectral_slab(w1*u.micron, w2*u.micron)
med1 = subcube2.median(axis=0)  

w1=(1+(3/res))*H1
w2=(1+(4/res))*H1
subcube2 = cube.spectral_slab(w1*u.micron, w2*u.micron)
med2 = subcube2.median(axis=0)  

med=(med2+med1)/2


cube=subcube-med

cube0=cube.moment(order=0)  


w1=(1-(2.5/res))*H1
w2=(1+(2.5/res))*H1
subcube = cube.spectral_slab(w1*u.micron, w2*u.micron)
print(w1,w2)
w1=(1-(3.5/res))*H1
w2=(1-(2.5/res))*H1
subcube2 = cube.spectral_slab(w1*u.micron, w2*u.micron)
med1 = subcube2.median(axis=0)  
print(w1,w2)
w1=(1+(2.5/res))*H1
w2=(1+(3.5/res))*H1
print(w1,w2)
subcube2 = cube.spectral_slab(w1*u.micron, w2*u.micron)
med2 = subcube2.median(axis=0)  

med=(med2+med1)/2

sub_f=subcube-med
moment_0 = sub_f.moment(order=0) 

# In[2]:


vel_cube = cube.with_spectral_unit(unit = 'km/s', velocity_convention='radio', rest_value=H1*u.micron)
vel_cube_err = subcube_err.with_spectral_unit(unit = 'km/s', velocity_convention='radio', rest_value=H1*u.micron)



dd = cube
from astropy.wcs import WCS

ww = WCS(dd.header).celestial

from astropy.coordinates import SkyCoord

from pvextractor import extract_pv_slice, Path



c2a=[166.6925932,-77.3762467]
c1a=[166.6937514,-77.3753579]
skypath2 = Path(SkyCoord([c1a[0],c2a[0]]*u.deg, [c1a[1],c2a[1]]*u.deg, frame='fk5'), width=1.4*u.arcsec)




RA = 166.6932418   
Dec =-77.3757983

from astropy.visualization import ImageNormalize
from astropy.visualization import LogStretch
from matplotlib.colors import LogNorm
coord = SkyCoord(ra=RA,dec=Dec, unit=(u.deg, u.deg)) # replace 'RA' and 'Dec' with the actual coordinates

norm = ImageNormalize(vmin=0.3,vmax=10,stretch=LogStretch())
x, y = ww.world_to_pixel(coord)

plt.figure(figsize=[12,10])
ax = plt.subplot(111, projection=ww)
ax.imshow(moment_0.data,interpolation='bilinear',cmap='Blues',norm =norm )

#plt.scatter(x,y,c='k',s=510)
#plt.scatter(x,y,c='w',s=200)
#ax.add_artist(path5.as_artist(0.5, wcs=ww), )

skypath2.show_on_axis(ax, spacing=1, color='r',linestyle=':',linewidth=0.3,alpha=0.2)
plt.show()

plt.xlabel('RA')
plt.ylabel('Dec')
plt.savefig('pv_Jet_fe_534_region.png',bbox_inches="tight")
# In[5]:

    
from pvextractor import extract_pv_slice, extract_pv_slice_err, Path
def gaussian_plus_line(x, amp, cen, wid, b):
    return amp * np.exp(-(x-cen)**2 / (2*wid**2)) + b



pvdiagram = extract_pv_slice(cube=vel_cube, path=skypath2,spacing=1)
pvdiagram_err= extract_pv_slice_err(cube=vel_cube_err, path=skypath2,spacing=1)

ww = WCS(pvdiagram.header)




cc = SkyCoord(RA*u.deg, Dec*u.deg, frame='icrs')
#cc = cc[0]
c1 = SkyCoord(c1a[0]*u.deg, c1a[1]*u.deg, frame='icrs')
c2 = SkyCoord(c2a[0]*u.deg, c2a[1]*u.deg, frame='icrs')
sep = c1.separation(cc)
#%%
pvdiagram = extract_pv_slice(cube=vel_cube, path=skypath2)
ww = WCS(pvdiagram.header)

new_wcs = ww
new_wcs.wcs.crval = np.array([-sep.value*189,  ww.wcs.crval[1]])
new_wcs.wcs.cdelt = np.array([ ww.wcs.cdelt[0]*189,  ww.wcs.cdelt[1]])



data=pvdiagram.data
data_err= pvdiagram_err.data

data[data<0]=1



#%%



plt.figure(figsize=[12,12])
#fits.setval("iras16253_rob-updatedWCS.fits",'CRVAL1', value=hdr_cr1, ext=1)
ax = plt.subplot(111, projection=new_wcs)
ax.imshow(data_err,norm=LogNorm())

#plt.scatter(0,4,s=100,c='k')

levels = np.array([0.1,0.2,0.3,0.5,0.7,0.9,0.95])*np.max(data)
#plt.contour(data,levels=levels,colors='r')
ax0 = ax.coords[0]
ax0.set_format_unit(u.arcsec)
ax1 = ax.coords[1]
ax1.set_ticks(np.arange(-150.0, 150.0, 30)*u.km/u.s)
ax1.set_format_unit(1000*u.km/u.s)

ax.set_aspect(1.5)

ax.set_ylabel("Velocity [km/s]")
ax.set_xlabel("Offset [au]")
plt.savefig('PV_err5.png')





#%%

#from lmfit import Model

import lmfit

x=np.linspace(0,np.shape(data)[0]-1,np.shape(data)[0])
plt.figure(figsize=[12,12])
xx=[]
ya=[]
yae=[]
ya_err=[]
yae_err=[]
x2=np.linspace(0,np.shape(data)[0]-1,100)
ii=ac
while(ii<np.shape(data)[1]):
#while(ii<26):
    model = lmfit.Model(gaussian_plus_line)
    def objective(params, x, data, errors):
        model_values = gaussian_plus_line(x, **params.valuesdict())
        residuals = (data - model_values) / errors
        return residuals
    
    params = model.make_params(amp=np.max(data[:,ii])*1.0, cen=np.argmax(data[:,ii])*1.0, b=0)
    params.add('wid', value=1, min=0.5,max=4)
    #result = model.fit(data[:,ii], params, x=x)
    
    minimizer = lmfit.Minimizer(objective, params, fcn_args=(x, data[:,ii], data_err[:,ii]))
    result = minimizer.minimize()


    
    plt.errorbar(x,data[:,ii],yerr=data_err[:,ii])
    xx.append(ii)
    ya.append(result.params['cen'].value)
    yae.append(result.params['wid'].value)
    ya_err.append(result.params['cen'].stderr)
    yae_err.append(result.params['wid'].stderr)
    plt.plot(x2, gaussian_plus_line(x2, **result.params), 'r-', label='Fitted Curve')
    plt.text(x2[0], 1.05*np.max(data[:,ii]), 'slices = ' + str(ii),size=15)
    plt.text(x2[0], 0.9*np.max(data[:,ii]), 'FWHM = ' + str(np.around(result.params['wid'].value*2.355,2)) + '$\pm$'
             + str(np.around(result.params['wid'].stderr*2.355,2))+'um',size=15)
    plt.text(x2[0], 0.8*np.max(data[:,ii]), 'Cen Vel = ' + str(np.around(result.params['cen'].value,2)) + 
             '$\pm$'+ str(np.around(result.params['cen'].stderr,2)) +'um',size=15)
    #plt.plot(x, result.best_fit, label='Fit')
    #ym = gauss(x2, popt[0], popt[1], popt[2])
    #plt.plot(x2, ym, c='r', label='Best fit')
    plt.savefig('test/PV1_'+str(ii)+'.png')
    plt.clf()
    ii=ii+1

#%%
plt.figure(figsize=[12,12])

ax = plt.subplot(111, projection=new_wcs)
ax.imshow(pvdiagram.data)
x=np.linspace(0,np.shape(data)[0]-1,np.shape(data)[0])
#plt.scatter(0,4,s=100,c='k')

levels = np.array([0.1,0.2,0.3,0.5,0.7,0.9,0.95])*np.max(data)
#plt.contour(data,levels=levels,colors='r')
ax0 = ax.coords[0]
ax0.set_format_unit(u.arcsec)
ax1 = ax.coords[1]
ax1.set_ticks(np.arange(-450.0, 450.0, 50)*u.km/u.s)
ax1.set_format_unit(1000*u.km/u.s)


ax.set_ylabel("Velocity [km/s]")
ax.set_xlabel("Offset [au]")

ax.set_aspect(2)

yae=np.array(yae)*2.355

xx=np.array(xx)
ya=np.array(ya)
ya_err=np.array(ya_err)
yae=np.array(yae)
yae_err=np.array(yae_err)

plt.errorbar(xx,ya,yerr=ya_err,color='k',ms=5,fmt='--o')

plt.tight_layout()
plt.show()
plt.savefig('PV_5.png')

#%%

plt.figure()

nx,ny=new_wcs.pixel_to_world_values(xx,ya)

plt.plot(nx*3600,ny/1000,'--o')
plt.ylabel("Velocity [km/s]")
plt.xlabel("Offset [arcsec]")

#%%

plt.figure(figsize=[12,12])




ii=ac

a,x1=new_wcs.pixel_to_world_values(0,x)
x1=x1/1000.0
x2=np.linspace(0,np.shape(data)[0]-1,100)
a,x2=new_wcs.pixel_to_world_values(0,x2)

x2=x2/1000.0

ya2=[]
ya2_err=[]
yae2=[]
yae2_err=[]


while(ii<np.shape(data)[1]):
    plt.plot(x1,data[:,ii])
    model = lmfit.Model(gaussian_plus_line)
    def objective(params, x, data, errors):
        model_values = gaussian_plus_line(x, **params.valuesdict())
        residuals = (data - model_values) / errors
        return residuals
    
    params = model.make_params(amp=np.max(data[:,ii])*1.0, cen=np.argmax(data[:,ii])*1.0, b=0)
    params.add('wid', value=44, min=20,max=80)
    params.add('amp', value=np.max(data[:,ii])*1.0, min=1)
    #result = model.fit(data[:,ii], params, x=x)
    
    minimizer = lmfit.Minimizer(objective, params, fcn_args=(x1, data[:,ii], data_err[:,ii]))
    result = minimizer.minimize()
    ya2.append(result.params['cen'].value)
    ya2_err.append(result.params['cen'].stderr)
    yae2.append(2.355*result.params['wid'].value)
    yae2_err.append(2.355*result.params['wid'].stderr)
    plt.errorbar(x1, data[:,ii], data_err[:,ii], fmt='bo', label='Data with Errors')
    plt.plot(x2, gaussian_plus_line(x2, **result.params), 'r-', label='Fitted Curve')
    plt.text(x2[0], 1.05*np.max(data[:,ii]), 'slices = ' + str(ii),size=15)
    plt.text(x2[0], 0.9*np.max(data[:,ii]), 'FWHM = ' + str(np.around(result.params['wid'].value*2.355,2)) + '$\pm$'
             + str(np.around(result.params['wid'].stderr*2.355,2))+'km/s',size=15)
    plt.text(x2[0], 0.8*np.max(data[:,ii]), 'Cen Vel = ' + str(np.around(result.params['cen'].value,2)) + 
             '$\pm$'+ str(np.around(result.params['cen'].stderr,2)) +'km/s',size=15)
    plt.xlabel('Velocity (km/s)')
    plt.ylabel('Flux')


    plt.clf()
    ii=ii+1

nx=np.array(nx)
ya2=np.array(ya2)
ya2_err=np.array(ya2_err)
yae2=np.array(yae2)
yae2_err=np.array(yae2_err)

#%%
plt.figure(figsize=[12,10])
plt.errorbar(nx*3600,yae2,yerr=yae2_err)
plt.xlabel("Offset [au]")
plt.ylabel("FWHM (km/s)")

#%%
plt.figure(figsize=[12,10])
#plt.scatter(nx[(yae2/yae2_err)>3]*3600,ya2)
plt.text(-270,40,'(b)')
plt.errorbar(nx*3600,ya2,yerr=ya2_err)
plt.ylabel("Velocity [km/s]")
plt.xlabel("Offset [au]")
plt.savefig('Velocity_5.png',bbox_inches="tight")

import pandas as pd

df2=pd.DataFrame()

df2['Sep'] = nx*3600

df2['Vel'] = ya2
df2['Vel_err'] = ya2_err

df2['FWHM'] = yae2
df2['FWHM_err'] = yae2_err

df2.to_csv('PV_Fe_534.csv')
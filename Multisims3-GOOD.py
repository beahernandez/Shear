
# coding: utf-8

# In[1]:

import sys
import os
import math
import numpy
import galsim
import numpy as np
from astropy.io import ascii


# In[2]:

sky_level = 0.8   #euclid
pixel_scale = 0.05  #1               # arcsec / pixel
noise_variance = (sky_level*pixel_scale**2)**2            # ADU^2

gal_flux_min =1.e4 # 3. #1.e4     # Range for galaxy flux
gal_flux_max =1.e6# 3.e1 #1.e6  
gal_hlr_min = 0.6       # arcsec
gal_hlr_max = 1.3       # arcsec
gal_e_min = 0.          # Range for ellipticity
gal_e_max = 0.8

psf_fwhm = 0.1 #0.65         # arcsec
psf_beta = 5 #Moffat profile

nx_tiles = 20                   #
ny_tiles = 20                   #
stamp_xsize = 100                #
stamp_ysize = 100                #

gal_ellip_max = 0.8 
gal_ellip_rms = 0.3 

disk_n = 1.5
disk_r0 = 0.85
#2.3

random_seed = 6424512
#gal_resolution = 0.98 
gal_re = 0.1 #0.85

gal_g1=0.05
gal_g2=0.07

shift_radius_sq = 1.0 
print 'done'

#SB Error: fourierDraw() requires an FFT that is too large, 6144
#If you can handle the large FFT, you may update gsparams.maximum_fft_size.
# In[ ]:

gal_g1=-0.4
gal_g2=-0.45
#cat_file_name ='galsim_default_input.asc'
#cat = galsim.Catalog(cat_file_name)
shears1=[]
shears2=[]
numbers=[]
ellipticities=[]
size=[]
flux=[]
mag=[]
sersic_indices=[]
print 'hey'
indices=[]

gara=galsim.GaussianDeviate(random_seed,mean=2,sigma=1)
gasn=galsim.GaussianDeviate(random_seed,mean=0.05,sigma=1.) #ORIGINAL:gasn=galsim.GaussianDeviate(random_seed,mean=20,sigma=0.5) #0.5 #BRIGHT:gasn=galsim.GaussianDeviate(random_seed,mean=0.05,sigma=1.) #

#gara=galsim.GaussianDeviate(mean=3,sigma=3.)
#gasn=galsim.GaussianDeviate(mean=1.,sigma=2.)

for n in range(50):
    tmp_indices=[]
    images = []

    gal_image = galsim.ImageF(stamp_xsize * nx_tiles-1 , stamp_ysize * ny_tiles-1,
                              scale=pixel_scale)
    gal_image_rot = galsim.ImageF(stamp_xsize * nx_tiles-1 , stamp_ysize * ny_tiles-1,
                              scale=pixel_scale)

    psf_image = galsim.ImageF(stamp_xsize * nx_tiles-1 , stamp_ysize * ny_tiles-1,
                              scale=pixel_scale)
    k=1
    big_fft_params = galsim.GSParams(maximum_fft_size=20000)

    gal_g1+=0.016
    gal_g2+=0.018
    print gal_g1
    shears1.append(gal_g1)
    shears2.append(gal_g2)
    numbers.append('%03d'%n)
    for iy in range(ny_tiles):
        for ix in range(nx_tiles):
            
            ud = galsim.UniformDeviate(random_seed)
            beta = ud() * 2. * math.pi * galsim.radians
            gd = galsim.GaussianDeviate(ud, sigma=gal_ellip_rms)
            val = gd()
            ellip = math.fabs(val)
            if ellip>gal_ellip_max:
                ellip=gal_ellip_max
            #print ellip
            #ellipticities.append(ellip)
            b = galsim.BoundsI(ix*stamp_xsize +1, (ix+1)*stamp_xsize -1, 
                               iy*stamp_ysize +1, (iy+1)*stamp_ysize-1)
            sub_gal_image = gal_image[b]
            sub_gal_image_rot = gal_image_rot[b]
            sub_psf_image = psf_image[b]
        
            rsq = 2 * shift_radius_sq #Random shift
            while (rsq > shift_radius_sq):
                dx = (2*ud()-1) * shift_radius_sq
                dy = (2*ud()-1) * shift_radius_sq
                rsq = dx**2 + dy**2 
        
            psf = galsim.Gaussian(fwhm = psf_fwhm, gsparams=big_fft_params)
            #psf = galsim.Moffat(beta=psf_beta, flux=1., fwhm=psf_fwhm)
        
            psf_re = psf.getHalfLightRadius()
            #gal_re = psf_re * gal_resolution
            
            reu = np.abs(gara())# 0.2+ud()*()
            #print reu*disk_r0
        
            gal_flux=np.abs(gasn())#gal_flux_min+ud()*(gal_flux_max-gal_flux_min)
            
            gal_mag=-2.5*np.log10(gal_flux/(0.8*0.05**2))+31.5
            print gal_flux,gal_mag
            
            sersic_idx=0.5+ud()*(4.-0.5)
            re=gal_re*reu
            bulge = galsim.Sersic(n=3, half_light_radius=re)
            disk = galsim.Exponential(half_light_radius=re*2)
            gal = bulge + disk

            #gal = galsim.Sersic(sersic_idx, flux=gal_flux, half_light_radius=re)
            
            tmp_indices.append(sersic_idx)
            
            #gal = galsim.DeVaucouleurs(flux=gal_flux, half_light_radius=re)
            gal = gal.shear(e=ellip, beta=beta)
            gal_rot=gal
            gal_rot=gal_rot.rotate(90*galsim.degrees)

        
            gal = gal.shear(g1=gal_g1, g2=gal_g2) #Add shear
          
            gal_rot = gal_rot.shear(g1=gal_g1, g2=gal_g2) #Add shear
    
            final = galsim.Convolve([psf, gal]) #Final galaxy and psf
            final = final.shift(dx,dy)
        
            final_rot = galsim.Convolve([psf, gal_rot]) #Final galaxy and psf
            final_rot = final_rot.shift(dx,dy)
        
            psf_flux=gal_flux_min*10**(ud()*(np.log10(gal_flux_max/gal_flux_min)))
            #print psf_flux
            star = psf.withFlux(psf_flux)
            star = star.shift(dx,dy)
        
            sub_gal_image += sky_level * pixel_scale**2
        
            final.drawImage(sub_gal_image, scale=pixel_scale)  #Draws image
        
            final_rot.drawImage(sub_gal_image_rot, scale=pixel_scale)  #Draws image
        
            star.drawImage(sub_psf_image, scale=pixel_scale)
        
            ellipticities.append(ellip)         
            size.append(re)
            sersic_indices.append(sersic_idx)
            flux.append(gal_flux)
            mag.append(gal_mag)
        
   
            k+=1
    
    t_image=[ellipticities,size,sersic_indices,flux,mag] #, ellipticities]
    ascii.write(t_image, 'HST_disk_bulge_2_sims%03d_galaxies_00_input.txt'%n)
    noise = galsim.GaussianNoise(ud, sigma=math.sqrt(noise_variance))
    #noise = galsim.CCDNoise(ud, gain=gain, read_noise=read_noise)
    gal_image.addNoise(noise)
    gal_image_rot.addNoise(noise)
    psf_image.addNoise(noise)
    if n==1:
        psf_file_name='HST_disk_bulge_2_sims_starfield.fits'
        psf_image.write(psf_file_name)
        
    gal_file_name='HST_disk_bulge_2_sims%03d_galaxies_00.fits'%n
    gal_file_name_rot='HST_disk_bulge_2_sims%03d_galaxies_90.fits'%n

    
    gal_image_rot.write(gal_file_name_rot)
    gal_image.write(gal_file_name)   
    
    indices.append(tmp_indices)

    
ascii.write(indices)


# In[17]:

tables=[numbers,shears1,shears2] #, ellipticities]
ascii.write(tables, 'input_shear_disk_bulge_2.txt')
ascii.write([ellipticities,size,flux], 'input_parameters_disk_bulge_2.txt')



# In[ ]:




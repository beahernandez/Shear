from astropy.io import fits
from scipy.ndimage.filters import gaussian_filter
import os

from skimage import io, color
import matplotlib.pyplot as plt
import numpy as np
from skimage import exposure
import pylab

def convolve2d_sub(image, kernel):
    # This function which takes an image and a kernel 
    # and returns the convolution of them
    # Args:
    #   image: a numpy array of size [image_height, image_width].
    #   kernel: a numpy array of size [kernel_height, kernel_width].
    # Returns:
    #   a numpy array of size [image_height, image_width] (convolution output).
    
    kernel = np.flipud(np.fliplr(kernel))    # Flip the kernel
    output = np.zeros_like(image)            # convolution output
    # Add zero padding to the input image
    image_padded = np.zeros((image.shape[0] + 8, image.shape[1] + 8))   
    image_padded[4:-4, 4:-4] = image
    for x_tmp in range(image.shape[1]):     # Loop over every pixel of the image
        for y_tmp in range(image.shape[0]):
            
            x = 3*(x_tmp/3)
            y = 3*(y_tmp/3)
            print x_tmp, x,y_tmp,y
            # element-wise multiplication of the kernel and the image
            output[y_tmp,x_tmp]=(kernel*image_padded[y:y+9,x:x+9]).sum()        
    return output


psf=fits.getdata('/export/data1/bhernandez/PhD/PSF_HST/Tiny_Tim/result00_subsampled3.fits')
kernel=[[0.0161,0.0161,0.0161, 0.0714,0.0714,0.0714,0.0161,0.0161, 0.0161],[0.0161,0.0161,0.0161, 0.0714,0.0714,0.0714,0.0161,0.0161, 0.0161],[0.0161,0.0161,0.0161, 0.0714,0.0714,0.0714,0.0161,0.0161, 0.0161],[0.0714,0.0714,0.0714,0.6500,0.6500,0.6500,0.0714,0.0714,0.0714],[0.0714,0.0714,0.0714,0.6500,0.6500,0.6500,0.0714,0.0714,0.0714],[0.0714,0.0714,0.0714,0.6500,0.6500,0.6500,0.0714,0.0714,0.0714],[0.0161,0.0161,0.0161, 0.0714,0.0714,0.0714,0.0161,0.0161, 0.0161],[0.0161,0.0161,0.0161, 0.0714,0.0714,0.0714,0.0161,0.0161, 0.0161],[0.0161,0.0161,0.0161, 0.0714,0.0714,0.0714,0.0161,0.0161, 0.0161]]
smoothed=convolve2d_sub(psf,kernel)#gaussian_filter(psf, 2)

hdu = fits.PrimaryHDU(smoothed)
hdul = fits.HDUList([hdu])

os.system('rm /export/data1/bhernandez/PhD/PSF_HST/Tiny_Tim/result00_subsampled3_cdgstride.fits')
hdul.writeto('/export/data1/bhernandez/PhD/PSF_HST/Tiny_Tim/result00_subsampled3_cdgstride.fits')
from astropy.io import fits
from astropy.io import ascii
from matplotlib import pylab as plt
import numpy as np
import glob
import os
from numpy import random

#psf='I'

import argparse

#psf='I'

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path', default = '',
                        type = str, help = 'Path to simulations')

args = parser.parse_args()
path=args.path
fi=glob.glob(path+'*_galaxies_90_filter_calib.cat')
fi.sort()

if not os.path.exists(path+'/final_cats_matched/'):
    os.mkdir(path+'/final_cats_matched/')
    
n1=-1

for f in fi:
    n1+=1
    num='%03d'%n1
    #num=f[81:84] #f[61:64]#num=f[67:70] #num=f[59:62] #f[63:66] #
    print num
    name=glob.glob(path+'/*'+num+'_galaxies_00_filter_calib.cat')
    hdulist = fits.open(name[0])
    cat_0=hdulist[3].data
    print len(cat_0)
    #cat_0=cat_0[cat_0['snratio']>5.] #cat_0['snratio'] cat_0['FLUX_ISO']/cat_0['FLUXERR_ISO']>5.] #
    #hdulist.close()
    #hdulist = fits.open('/vol/euclid1/euclid1_3/bhernandez/galsims/ellip_Poisson/sim_001_24.5_a45/sim_001_24.5_a45_filter_calib.cat')
    #cat_45=hdulist[3].data
    #hdulist.close()
    name=glob.glob(path+'/*'+num+'_galaxies_90_filter_calib.cat')
    hdulist = fits.open(name[0])
    cat_90=hdulist[3].data
    #cat_90=cat_90[cat_90['snratio']>5.]# cat_90['FLUX_ISO']/cat_90['FLUXERR_ISO']>5.] #
    hdulist.close()
    #hdulist = fits.open('/vol/euclid1/euclid1_3/bhernandez/galsims/ellip_Poisson/sim_001_24.5_a135/sim_001_24.5_a135_filter_calib.cat')
    #cat_135=hdulist[3].data
    #hdulist.close()
    #hdulist_er = fits.open('/vol/euclid1/euclid1_3/bhernandez/galsims/lens9/input9_sims'+num+'_galaxies_00/input9_sims'+num+'_galaxies_00_filter_calib.cat')
    #hdulist_er_rot = fits.open('/vol/euclid1/euclid1_3/bhernandez/galsims/lens9/input9_sims'+num+'_galaxies_90/input9_sims'+num+'_galaxies_90_filter_calib.cat')
    #errors=hdulist_er[3].data
    #errors_rot=hdulist_er_rot[3].data
    
    l=[]
    e1iso_0=[]
    e2iso_0=[]
    e1iso_45=[]
    e2iso_45=[]
    e1iso_90=[]
    e2iso_90=[]
    e1iso_135=[]
    e2iso_135=[]
    NrIn_ar=[]
    deltae1=[]
    deltae2=[]
    n=0
    
    for obj in range(len(cat_0)):
	
	x=cat_0['x'][obj]
	y=cat_0['y'][obj]
	
	#NrIn=cat_0['SeqNr'][obj]
	#idx45 = np.where((x-10<cat_45['x']) & (cat_45['x']<x+10) & (y-10<cat_45['y'])&(cat_45['y']<y+10))
	idx90 = np.where((x-10<cat_90['x']) & (cat_90['x']<x+10) & (y-10<cat_90['y'])&(cat_90['y']<y+10))
	#idx135 = np.where((x-10<cat_135['x']) & (cat_135['x']<x+10) & (y-10<cat_135['y'])&(cat_135['y']<y+10))
	
	#idx=cat_rot['NrIn'].index(NrIn)
	#print np.size(idx45), np.size(idx90), np.size(idx135), n
	#if np.size(idx45)>0 and np.size(idx90)>0 and np.size(idx135)>0:
	print np.size(idx90), n
	if np.size(idx90)>0:
            
            
            m_corr=-0.078*(cat_0['snratio'][obj]/2.)**(-0.38)
	    e1iso_0.append(cat_0['e1iso'][obj]/(1+m_corr))
	    e2iso_0.append(cat_0['e2iso'][obj]/(1+m_corr))
	
	    #e1iso_45.append(cat_45['e1iso'][idx45])
	    #e2iso_45.append(cat_45['e2iso'][idx45])
	    
	    m_corr=-0.078*(cat_90['snratio'][idx90]/2.)**(-0.38)
	    
	    e1iso_90.append(cat_90['e1iso'][idx90]/(1+m_corr))
	    e2iso_90.append(cat_90['e2iso'][idx90]/(1+m_corr))
	    
	    #e1iso_135.append(cat_135['e1iso'][idx135])
	    #e2iso_135.append(cat_135['e2iso'][idx135])
	    
	    n+=1
	    #deltae1.append((errors['deltae1'][obj]+errors_rot['deltae1'][idx])/2.)
	    #deltae2.append((errors['deltae2'][obj]+errors_rot['deltae1'][idx])/2.)
	
	    #NrIn_ar.append(NrIn)


    #col1 = fits.Column(name='NrIn', format='I4', array=NrIn_ar)
    col2 = fits.Column(name='e1iso_0', format='E', array=e1iso_0)
    #col3 = fits.Column(name='e1iso_45', format='E', array=e1iso_45)
    col4 = fits.Column(name='e1iso_90', format='E', array=e1iso_90)
    #col5 = fits.Column(name='e1iso_135', format='E', array=e1iso_135)
    
    col6 = fits.Column(name='e2iso_0', format='E', array=e2iso_0)
    #col7 = fits.Column(name='e2iso_45', format='E', array=e2iso_45)
    col8 = fits.Column(name='e2iso_90', format='E', array=e2iso_90)
    #col9 = fits.Column(name='e2iso_135', format='E', array=e2iso_135)
    
    ###First approximation to errors###
    #col6 = fits.Column(name='deltae1', format='E', array=deltae1)
    #col7 = fits.Column(name='deltae2', format='E', array=deltae2)
    #cols = fits.ColDefs([col2,col3,col4,col5,col6,col7,col8,col9])
    cols = fits.ColDefs([col2,col4,col6,col8])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    hdu=fits.PrimaryHDU()
	
    hdulist = fits.HDUList([hdu,tbhdu])
    #os.remove('/vol/euclid1/euclid1_3/bhernandez/galsims/ellip_Poisson/final_cats_matched/test'+num+'.fits')
    hdulist.writeto(path+'/final_cats_matched/test_nocut_corr_'+num+'.fits')
    
#files=glob.glob()    
#data=fits.getdata('/vol/euclid1/euclid1_3/bhernandez/galsims/ellip_Poisson/final_cats_matched/test'+num+'.fits') #'/vol/euclid1/euclid1_3/bhernandez/Nsims/final_cats_matched/test.fits')
#e1=[]
#e2=[]
#for itera in range(100):
##Results2=np.polynomial.polynomial.polyfit(e1_in, e1_iso-e1_in,deg=1,w=1./e1_err)
    #boot=random.randint(0, high=len(data[0]), size=len(data[0]))
    #for b in boot:
	#e1.append((data['e1iso_0'][b]+data['e1iso_90'][b])/2.)
	#e2.append((data['e2iso_0'][b]+data['e2iso_90'][b])/2.)
	##e1.append((data['e1iso_0'][b]+data['e1iso_45'][b]+data['e1iso_90'][b]+data['e1iso_135'][b])/4.)
	##e2.append((data['e2iso_0'][b]+data['e2iso_45'][b]+data['e2iso_90'][b]+data['e2iso_135'][b])/4.)

#print np.mean(e1), np.sqrt(np.var(e1))/len(data[0]), np.mean(e2), np.sqrt(np.var(e2))/len(data[0])
    
    

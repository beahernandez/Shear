import numpy as np
from astropy.io import fits
from astropy.io import ascii
import glob
from matplotlib import pylab as plt
from scipy import stats

import argparse
from matplotlib.axes import Axes as ax

import matplotlib
####Configuring matplotlib##############################
matplotlib.rcParams['xtick.major.size']= 6
matplotlib.rcParams['ytick.major.size']= 6
matplotlib.rcParams['ytick.minor.size']= 4
matplotlib.rcParams['xtick.minor.size']= 4
matplotlib.rcParams['xtick.major.width']= 1.
matplotlib.rcParams['ytick.major.width']= 1.
matplotlib.rcParams['ytick.minor.width']= 1.
matplotlib.rcParams['xtick.minor.width']= 1.
matplotlib.rcParams['axes.linewidth']= 2.0
matplotlib.rcParams['lines.linewidth']=2
matplotlib.rcParams['lines.markersize']=8
matplotlib.rcParams['lines.markeredgewidth']=1
matplotlib.rcParams['xtick.labelsize']= 15
matplotlib.rcParams['ytick.labelsize']= 15
###########################################################


matplotlib.rcParams['patch.linewidth']=2.

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path', default = '',
                        type = str, help = 'Path to simulations')
parser.add_argument('-r', '--ref', default = True,
                        type = bool, help = 'Include reference plot')

args = parser.parse_args()
path=args.path
ref=args.ref
fi=glob.glob(path+'*filter_calib.cat') #00_
fi.sort()

sn=[]
r=[]
mag=[]

###Zp=26.7 #For F814W
###t_exp=100. #In seconds
###gain=14. #HST like

#From fitting:
A=2.52722419076 
Zp=27.1007978301

#Creates the data distribution
n1=-1
for f in fi:
    n1+=1
    num='%03d'%n1
    hdulist = fits.open(f)
    data=hdulist[3].data
    #print len(data)
    data=data[data['FLUX_ISO']>(Zp-2.5*(np.log10(24.2)/A))]
    hdulist.close()
    for n in range(len(data['snratio'])):
        sn.append(data['FLUX_ISO'][n]/data['FLUXERR_ISO'][n]) #sn.append(data['snratio'][n])#data['FLUX_ISO'][n]/data['FLUXERR_ISO'][n])#
        r.append(data['rh'][n])
        mag.append(Zp-2.5*(np.log10(data['FLUX_ISO'][n]/A))) #gain/t_exp)))
sn=np.array(sn)
r=np.array(r)
mag=np.array(mag)
    

    
#Uses the real catalogue as comparison
tmp=ascii.read('/vol/euclid1/euclid1_3/bhernandez/galsims/HAWKI/setup/good_catalogue.cat')
hawki_cat= tmp[(tmp['col6']>21.0) & (tmp['col6']<24.2) & (tmp['col7']>10.) & (np.sqrt(tmp['col3'] **2 +tmp['col4']**2)<1.0)]

#Plots sn
plt.hist(sn[sn<80], bins=50, normed=True,histtype='stepfilled',color='0.5', label='Simulations')
if ref==True:
    plt.hist(hawki_cat['col7'][hawki_cat['col7']<50], bins=50, normed=True,histtype='step',color='r',label='reference')
plt.xlabel(r'$S/N_{ksb}$',size=20)
plt.ylabel(r'Normalized histogram',size=20)
plt.legend(numpoints=1)
plt.savefig(path+'sn_distribution_ksb_mod.eps')
plt.show()

#Plots radius
plt.hist(r, bins=50,normed=True,histtype='stepfilled',color='0.5', label='Simulations')
if ref==True:
    plt.hist(hawki_cat['col8'], bins=50, normed=True,histtype='step',color='r',label='reference')
plt.xlabel(r'Radius [pix]',size=20)
plt.ylabel(r'Normalized histogram',size=20)
plt.xlim(0,14)
plt.legend(numpoints=1)
plt.savefig(path+'fluxradius_distribution_ksb_mod.eps')
plt.show()

#Plots magnitudes
plt.hist(mag, bins=50,normed=True,histtype='stepfilled',color='0.5', label='Simulations')
if ref==True:
    plt.hist(hawki_cat['col6'], bins=50, normed=True,histtype='step',color='r',label='reference')
plt.xlabel(r'MAG_AUTO',size=20)
plt.ylabel(r'Normalized histogram',size=20)
#plt.xlim(0,14)
plt.legend(numpoints=1)
plt.savefig(path+'mag_distribution_ksb_mod.eps')
plt.show()

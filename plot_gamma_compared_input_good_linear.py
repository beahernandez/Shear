from astropy.io import ascii
from astropy.io import fits
from matplotlib import pylab as plt
from numpy import random
from scipy import stats
from scipy import optimize
import emcee
import numpy as np
import glob
import argparse
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
matplotlib.rcParams['lines.linewidth']=1.5
matplotlib.rcParams['lines.markersize']=8
matplotlib.rcParams['lines.markeredgewidth']=1
matplotlib.rcParams['xtick.labelsize']= 15
matplotlib.rcParams['ytick.labelsize']= 15
############################################################


##################################################################
## This code reads the matched catalogs and finds the average,  ##
## comparing it to the input to obtain the bias.                ##
## Bootstrap is used to compute the errors.                     ##
## A fitting procedure is used to compute bias.                 ##
##################################################################


#Reads the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path', default = '',
                        type = str, help = 'Path to simulations')

args = parser.parse_args()
path=args.path

#Finds the catalogs and sorts them by name
fi=glob.glob(path+'*_galaxies_90_filter_calib.cat')
fi.sort()

#Reads the input
input_cats=ascii.read(path+'input_shear.txt')

e1_in=input_cats['col2']#[7]
e2_in=input_cats['col3']#[8]


#Defines empty arrays needed for the loops
q_1=[]
q_2=[]
slopes_1=[]
intercepts_1=[]
slopes_2=[]
intercepts_2=[]
slopes_1_lin=[]
intercepts_1_lin=[]
slopes_2_lin=[]
intercepts_2_lin=[]

e1_iso=[]
e1_err=[]
e2_iso=[]
e2_err=[]

e1_iso_subs=[]
e1_err_subs=[]
e2_iso_subs=[]
e2_err_subs=[]

n1=-1 #Defines the counter for the loop. Starts in -1 so we can put the counter at the beginning.

#Loops over the files
for f in fi:
    print f
    n1+=1 #Counter
    num='%03d'%n1 #Formatting of the counter so we can use it to open files
    
    #Reads the matched catalogs
    data=fits.getdata(path+'final_cats_matched/test_nocut_corr_'+num+'.fits') #'/vol/euclid1/euclid1_3/bhernandez/galsims/lens9/final_cats_matched/input9_sims'+num+'.fits')
    
    e1_fin=[]
    e2_fin=[]
    
    #The mean is not computed using bootstrapping but it is just the mean, so we keep the original values here (matched)
    e1_or=(data['e1iso_0']+data['e1iso_90'])/2.
    e2_or=(data['e2iso_0']+data['e2iso_90'])/2.
    
    ###Bootstrapping####
    for itera in range(100):  
	boot=random.randint(0, high=len(data['e1iso_0']), size=len(data['e1iso_0'])) #Random indices for the bootstrapping

        #Matching of the bootstrap sample
        e1=(data['e1iso_0'][boot]+data['e1iso_90'][boot])/2.
        e2=(data['e2iso_0'][boot]+data['e2iso_90'][boot])/2.
        
        #Mean values of this sample:
	e1_fin.append(np.mean(e1))
	e2_fin.append(np.mean(e2))
    
    #Finds the catalogs for the SN# NOT USED HERE
    #name=glob.glob(path+'*'+num+'_galaxies_00_filter_calib.cat')
    #hdulist = fits.open(name[0]) #/vol/euclid1/euclid1_3/bhernandez/galsims/lens9/input9_sims'+num+'_galaxies_00/input9_sims'+num+'_galaxies_00_filter_calib.cat')
    #cat=hdulist[3].data
    #hdulist.close()
    #SN=cat['snratio'][b]
    #m_corr=-0.078*(SN/2.)**(-0.38)
    
    #Final results for this file
    e1_iso.append(np.mean(e1_or))
    e2_iso.append(np.mean(e2_or))
    e1_err.append(np.sqrt(np.var(e1_fin)))#/np.sqrt(len(data['e1iso_0'])))
    e2_err.append(np.sqrt(np.var(e2_fin)))#/np.sqrt(len(data['e2iso_0'])))


######



##################SN correction#### NOT USED HERE
#################sn_list=[]

#################n1=0

#################for f in fi:
    #################n1+=1
    #################num='%03d'%n1
    ##################num=f[81:84] #f[61:64]#num=f[67:70] #num=f[59:62] #f[63:66] #
    ####################Plotting comparing to input shear###
    #################name=glob.glob(path+'*'+num+'_galaxies_00_filter_calib.cat')
    #################hdulist = fits.open(name[0]) #/vol/euclid1/euclid1_3/bhernandez/galsims/lens9/input9_sims'+num+'_galaxies_00/input9_sims'+num+'_galaxies_00_filter_calib.cat')
    #################cat=hdulist[3].data
    #################hdulist.close()
    #################sn_list.append(2*cat['snratio'])
    
#################sn_list=np.concatenate(sn_list).ravel()
#################sn_list.sort()

#################SN=np.mean(sn_list)


#################m_corr=-0.078*(SN/2.)**(-0.38)

##################



#Initiation of the plotting for the first axis
fig=plt.figure()
plt.errorbar(e1_in,e1_iso-e1_in,yerr=e1_err,fmt='ko', mfc='none')
plt.xlabel(r"$g_1^{\rm{inp}}$",size=25)
plt.ylabel(r"$g_1-g_1^{\rm{inp}}$",size=25) 

#Defines the arrays
e1_in_lin=e1_in #[(e1_in>-0.2) & (e1_in<0.2)]
e1_iso=np.array(e1_iso)
e1_iso_lin=e1_iso #[(e1_in>-0.2) & (e1_in<0.2)]

#Bootstrap of the fit
for itera in range(100):

    boot=random.randint(0, high=len(e1_in_lin), size=len(e1_in_lin)) #Bootstrap samples

    #Finds the values of the sample
    x_lin=e1_in_lin[boot]
    y_lin=e1_iso_lin[boot]-e1_in_lin[boot]
    
    #Performs a fit to these values
    Results_lin=stats.linregress(x_lin,y_lin)
    
    #Saves fit so we can compute the average and obtain errors
    slopes_1_lin.append(Results_lin[0])
    intercepts_1_lin.append(Results_lin[1])
    
    
#Final mean and error of the samples of the fit
m_lin=np.mean(slopes_1_lin)
m_err_lin=np.sqrt(np.var(slopes_1_lin))
n_lin=np.mean(intercepts_1_lin)
n_err_lin=np.sqrt(np.var(intercepts_1_lin))


#Continue plotting   
plt.tight_layout()
plt.plot(e1_in,m_lin*e1_in+n_lin,'b--')
plt.savefig(path+'final_cats_matched/input_comparison_gamma_1_lin_nocut_corr.eps')
#plt.show()
plt.close(fig)

#Print the fit parameters
print 'e1 lin',m_lin,m_err_lin,n_lin,n_err_lin

#Saves values of the fit to write to file later
e1_measured=[m_lin,m_err_lin,n_lin,n_err_lin]


#Initiation of the plotting for the second axis
fig=plt.figure()
plt.errorbar(e2_in,e2_iso-e2_in,yerr=e2_err,fmt='ko', mfc='none')
plt.xlabel(r"$g_2^{\rm{inp}}$",size=25)
plt.ylabel(r"$g_2-g_2^{\rm{inp}}$",size=25) 

#Defines the arrays
e2_in_lin=e2_in#[(e2_in>-0.2) & (e2_in<0.2)]
e2_iso=np.array(e2_iso)
e2_iso_lin=e2_iso#[(e2_in>-0.2) & (e2_in<0.2)]

#Bootstrap of the fit
for itera in range(100):

    boot=random.randint(0, high=len(e2_in_lin), size=len(e2_in_lin)) #Bootstrap samples
    
    #Finds the values of the sample
    x_lin=e2_in_lin[boot]
    y_lin=e2_iso_lin[boot]-e2_in_lin[boot]
    
    #Performs a fit to these values
    Results_lin=stats.linregress(x_lin,y_lin)
    
    #Saves fit so we can compute the average and obtain errors
    slopes_2_lin.append(Results_lin[0])
    intercepts_2_lin.append(Results_lin[1])
    

#Final mean and error of the samples of the fit
m_lin=np.mean(slopes_2_lin)
m_err_lin=np.sqrt(np.var(slopes_2_lin))
n_lin=np.mean(intercepts_2_lin)
n_err_lin=np.sqrt(np.var(intercepts_2_lin))

#Continue plotting
plt.tight_layout()
plt.plot(e2_in,m_lin*e2_in+n_lin,'b--')
plt.savefig(path+'final_cats_matched/input_comparison_gamma_2_lin_nocut_corr.eps')
#plt.show()
plt.close(fig)

#Print the fit parameters
print 'e2 lin',m_lin,m_err_lin,n_lin,n_err_lin

#Saves values of the fit to write to file later
e2_measured=[m_lin,m_err_lin,n_lin,n_err_lin]

#Prints to file to save the fit parameters
with open(path+'values_lin_nocut_corr.txt', 'w') as f:#_ksbsn5
	print >> f, e1_measured, e2_measured

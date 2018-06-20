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




#####################################################
#This code reads the matched catalogs and finds the average, comparing it to the input to obtain the bias.
#This is done for the different intrinsic ellipticities. The catalogs should be in the same path (given as an input) under three different folders named 020, 025, 030, 035. Then they are plot together so we can compare them.
#####################################################

#Reads the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path', default = '',
                        type = str, help = 'Path to simulations')

args = parser.parse_args()
path=args.path



#Defines the empty arrays used in the bootstraping
M1=[]
M1_err=[]
N1=[]
N1_err=[]
M2=[]
M2_err=[]
N2=[]
N2_err=[]

E1_iso=[]
E1_err=[]
E2_iso=[]
E2_err=[]



#Loop for the different ellipticities.
for ellip in ['020','025','030','035']:
    
    
    #Finds the catalogs and sorts them by name
    fi=glob.glob(path+'KSB_cats_'+ellip+'/*_galaxies_90_filter_calib.cat')
    fi.sort()
    
    #Reads the input
    input_cats=ascii.read(path+'KSB_cats_'+ellip+'/input_shear.txt')

    e1_in=input_cats['col2']#[7]
    e2_in=input_cats['col3']#[8]
    
    #Defines more empty arrays needed for the loops
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
        data=fits.getdata(path+'KSB_cats_'+ellip+'/final_cats_matched/test_nocut_corr_'+num+'.fits') #'/vol/euclid1/euclid1_3/bhernandez/galsims/lens9/final_cats_matched/input9_sims'+num+'.fits')
    
        
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
            #print f, itera
            
        #Final mean and error of all samples
        e1_iso.append(np.mean(e1_or))
        e2_iso.append(np.mean(e2_or))
        e1_err.append(np.sqrt(np.var(e1_fin)))#/np.sqrt(len(data['e1iso_0'])))
        e2_err.append(np.sqrt(np.var(e2_fin)))#/np.sqrt(len(data['e2iso_0'])))


    ######



    ##SN correction####NOT USED HERE
    #sn_list=[]

    #n1=-1

    #for f in fi:
        #n1+=1
        #num='%03d'%n1
        ##num=f[81:84] #f[61:64]#num=f[67:70] #num=f[59:62] #f[63:66] #
        ####Plotting comparing to input shear###
        #name=glob.glob(path+'KSB_cats'+ellip+'/*'+num+'_galaxies_00_filter_calib.cat')
        #hdulist = fits.open(name[0]) #/vol/euclid1/euclid1_3/bhernandez/galsims/lens9/input9_sims'+num+'_galaxies_00/input9_sims'+num+'_galaxies_00_filter_calib.cat')
        #cat=hdulist[3].data
        #hdulist.close()
        #sn_list.append(2*cat['snratio'])
        
    #sn_list=np.concatenate(sn_list).ravel()
    #sn_list.sort()

    #SN=np.mean(sn_list)


    #m_corr=-0.078*(SN/2.)**(-0.38)

    ###################



    
    ###########Plotting and bias determination########
    
    #For first axis
    e1_in_lin=e1_in #[(e1_in>-0.2) & (e1_in<0.2)]
    e1_iso=np.array(e1_iso)
    e1_iso_lin=e1_iso #[(e1_in>-0.2) & (e1_in<0.2)]
    
    #Bootstrap of the fit
    for itera in range(100):
    
        boot=random.randint(0, high=len(e1_in_lin), size=len(e1_in_lin)) #Bootstrap sample

        #Finds the values of the sample
        x_lin=e1_in_lin[boot]
        y_lin=e1_iso_lin[boot]-e1_in_lin[boot]
        
        #Performs a fit to these values
        Results_lin=stats.linregress(x_lin,y_lin)
        
        #Saves fit so we can compute the average and obtain errors
        slopes_1_lin.append(Results_lin[0])
        intercepts_1_lin.append(Results_lin[1])

        
    #Final mean and error of the samples of the fit
    m1_lin=np.mean(slopes_1_lin)
    m1_err_lin=np.sqrt(np.var(slopes_1_lin))
    n1_lin=np.mean(intercepts_1_lin)
    n1_err_lin=np.sqrt(np.var(intercepts_1_lin))

    #For second axis
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
    m2_lin=np.mean(slopes_2_lin)
    m2_err_lin=np.sqrt(np.var(slopes_2_lin))
    n2_lin=np.mean(intercepts_2_lin)
    n2_err_lin=np.sqrt(np.var(intercepts_2_lin))
    
    #Arrays to save the values of the individual fits for each of the different ellipticities.
    M1.append(m1_lin)
    M1_err.append(m1_err_lin)
    N1.append(n1_lin)
    N1_err.append(n1_err_lin)
    
    M2.append(m2_lin)
    M2_err.append(m2_err_lin)
    N2.append(n2_lin)
    N2_err.append(n2_err_lin)
    
    #Lists with fit parameters for writting output
    e1_measured=[m1_lin,m1_err_lin,n1_lin,n1_err_lin]
    e2_measured=[m2_lin,m2_err_lin,n2_lin,n2_err_lin]

    #Output for each of the ellipticities individually
    with open(path+'KSB_cats_'+ellip+'/values_lin_nocut_comp_corr.txt', 'w') as f:
            print >> f, e1_measured, e2_measured
     
    #Arrays with the data to plot
    E1_iso.append(e1_iso)
    E1_err.append(e1_err)
    E2_iso.append(e2_iso)
    E2_err.append(e2_err)
    
    



#Plotting for first axis
fig=plt.figure()
plt.errorbar(e1_in,E1_iso[0]-e1_in,yerr=E1_err[0],color='#bfbfbf',fmt='.', mfc='none',label='0.20',elinewidth=1)#000000
plt.errorbar(e1_in+0.003,E1_iso[1]-e1_in,yerr=E1_err[1],color='#9999ff',fmt='x', mfc='none',label='0.25',elinewidth=1)#0000ff
plt.errorbar(e1_in+0.006,E1_iso[2]-e1_in,yerr=E1_err[2],color='#ff9999',fmt='+', mfc='none',label='0.30',elinewidth=1)#ff9999#ff0000
plt.errorbar(e1_in+0.009,E1_iso[3]-e1_in,yerr=E1_err[3],color='#80ffaa',fmt='1', mfc='none',label='0.35',elinewidth=1)#009933
plt.xlabel(r"$g_1^{\rm{inp}}$",size=20)
plt.ylabel(r"$g_1-g_1^{\rm{inp}}$",size=20) 

plt.plot(e1_in,M1[0]*e1_in+N1[0],'k--')
plt.plot(e1_in,M1[1]*e1_in+N1[1],'b--')
plt.plot(e1_in,M1[2]*e1_in+N1[2],'r--')
plt.plot(e1_in,M1[3]*e1_in+N1[3],'g--')

plt.tight_layout()

plt.legend(numpoints=1)
plt.savefig(path+'/ell_1_lin_nocut_corr.eps')
#plt.show()
plt.close(fig)


#Plotting for second axis
fig=plt.figure()
plt.errorbar(e2_in,E2_iso[0]-e2_in,yerr=E2_err[0],color='#bfbfbf',fmt='.', mfc='none',label='0.20',elinewidth=1)#000000
plt.errorbar(e2_in+0.003,E2_iso[1]-e2_in,yerr=E2_err[1],color='#9999ff',fmt='x', mfc='none',label='0.25',elinewidth=1)#0000ff
plt.errorbar(e2_in+0.006,E2_iso[2]-e2_in,yerr=E2_err[2],color='#ff9999',fmt='+', mfc='none',label='0.30',elinewidth=1)#ff9999#ff0000
plt.errorbar(e2_in+0.009,E2_iso[3]-e2_in,yerr=E2_err[3],color='#80ffaa',fmt='1', mfc='none',label='0.35',elinewidth=1)#009933
plt.xlabel(r"$g_2^{\rm{inp}}$",size=20)
plt.ylabel(r"$g_2-g_2^{\rm{inp}}$",size=20) 

plt.plot(e2_in,M2[0]*e2_in+N2[0],'k--')
plt.plot(e2_in,M2[1]*e2_in+N2[1],'b--')
plt.plot(e2_in,M2[2]*e2_in+N2[2],'r--')
plt.plot(e2_in,M2[3]*e2_in+N2[3],'g--')

plt.tight_layout()
plt.legend(numpoints=1)
plt.savefig(path+'/ell_2_lin_nocut_corr.eps')
#plt.show()
plt.close(fig)


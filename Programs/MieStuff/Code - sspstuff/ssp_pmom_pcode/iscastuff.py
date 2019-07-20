# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 20:32:14 2016

@author: Steven
"""

import pandas as pd
import numpy as np


def loadiscp(iscafilename, p11filename):
    # Read the isca file
    isca = pd.read_csv(iscafilename,delim_whitespace=True,header=None)
    a,b = np.shape(isca)
    NANG = 498 # number of scattering angles
    Nreff = 189 # number of sizes
    
    if (a==9261):
        print ('This is a 16-99 um file')    
        NIORS = 49 # number of wavelengths
    elif (a==74844):
        print ('This is a 0-15 um file')    
        NIORS = 396
    elif (a==84105):
        print ('This is a 0-99 um file')    
        NIORS = 445 # 49+396
    else:
        print ('Bad input size')
        return
    
    #print (NIORS*Nreff)


    # Extract un-multiplexed single-scattering parameters in the isca file
    wlen_long = isca[0]
    maxdim_long = isca[1]
    volume_long = isca[2]
    parea_long = isca[3]
    qext_long = isca[4]
    w0_long = isca[5]
    asym_long = isca[6]
    
    # Calculate the effective radius: V/A(proj) = 4/3 pi r^3 / pi r^2 so reff = 3 V/(A(proj)x4) = 3 V/(A(surf))
    # In other words, what we're calling here "reff" is actually the equivalent V/A radius,
    reff_long = 3*volume_long/(parea_long*4)
    
    # Mesh and list: Wavemumber
    wlen_mesh = np.reshape(wlen_long,(NIORS,Nreff)).T; #print (wlen_mesh[0:3,0:3])
    wnum_mesh = 1e4/wlen_mesh; #print (wnum_mesh[0,:])
    wnumlist = wnum_mesh[0,:]; #print (wnumlist)
    wnum_list = np.zeros((NIORS,1)); wnum_list[:,0] = wnumlist; #print ("wnum", "\n", wnum_list, wnum_list.shape)
    #print (np.shape(wnum_mesh))
    
    # Mesh and list: reff
    reff_mesh = np.reshape(reff_long,(NIORS,Nreff)).T; #print reff_mesh[0:3,0:3]
    refflist = reff_mesh[:,0]
    reff_list = np.zeros((Nreff,1)); reff_list[:,0] = refflist; #print ("reff", "\n", reff_list, reff_list.shape)
    #print (np.shape(reff_mesh))
    #print np.max(reff_mesh)
    
    # Mesh: Single-scattering albedo
    w0_mesh = np.reshape(w0_long,(NIORS,Nreff)).T; #print (w0_mesh[0:3,0:3])
    #print (np.shape(w0_mesh))
    
    
    # Mesh: Extinction efficiency
    qext_mesh = np.reshape(qext_long,(NIORS,Nreff)).T; #print (qext_mesh[0:3,0:3])
    #print (np.shape(qext_mesh))
    
    # Mesh: Asymmetry parameter
    asym_mesh = np.reshape(asym_long,(NIORS,Nreff)).T; #print (asym_mesh[0:3,0:3])
    #print (np.shape(asym_mesh))
    
    # List: Maxdim
    maxdim_mesh = np.reshape(maxdim_long,(NIORS,Nreff)).T; #print (maxdim_mesh[0:3,0:3])
    maxdimlist = maxdim_mesh[:,0]
    maxdim_list = np.zeros((Nreff,1)); maxdim_list[:,0] = maxdimlist;
    #print (np.shape(maxdim_mesh))
        
    # List: Volume
    volume_mesh = np.reshape(volume_long,(NIORS,Nreff)).T; #print (volume_mesh[0:3,0:3])
    volumelist = volume_mesh[:,0]
    volume_list = np.zeros((Nreff,1)); volume_list[:,0] = volumelist;
    #print (np.shape(volume_mesh))
    
    # List: Projected area
    parea_mesh = np.reshape(parea_long,(NIORS,Nreff)).T; #print (parea_mesh[0:3,0:3])
    parealist = parea_mesh[:,0]
    parea_list = np.zeros((Nreff,1)); parea_list[:,0] = parealist;
    #print (np.shape(parea_mesh))
        
    # Reading the phase function
    p11_frame = pd.read_csv(p11filename,delim_whitespace=True,header=None)
    #print (np.shape(p11_frame))
    
    # Extracting the phase angles (first row)
    ppa = p11_frame.ix[0,:]; #print(ppa[0:5])
    #print (np.shape(ppa))
    
    # Extracting the phase functions (after the first row)
    p11_long = np.array(p11_frame.ix[1:,:].T); #print (np.shape(p11_long))
    p11 = np.reshape(p11_long,(NANG,NIORS,Nreff)); # Number of angles, wavelengths, and particle sizes
    #print (np.shape(p11))  
    
    
    return \
    ppa, p11, \
    wnum_mesh, reff_mesh, w0_mesh, qext_mesh, asym_mesh, \
    wnum_list, reff_list, maxdim_list, volume_list, parea_list
    
    
    
    
    
    
    
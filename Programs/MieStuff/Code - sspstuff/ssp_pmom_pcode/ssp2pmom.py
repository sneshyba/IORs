
# coding: utf-8

# In[1]:

# Get resources
import numpy as np
get_ipython().magic(u'load_ext oct2py.ipython')


# In[3]:

# Get filenames
# filename = 'ssp_getpsd_T273_S331.txt'
# outfilename = 'test.nc'
filename =    "../IORs/iors_water_263 (Hybrid)/ssp_getpsd_T263_S331.txt"
outfilename = "../IORs/iors_water_263 (Hybrid)/ssp_getpsd_T263_S331_pmom.nc"

# This runs getiops like a script; input is
# %%octave -i filename -o f,Pnrm_psd_straight,iops_psd_straight,NANG,NIORS,Nreff,Nfrompsd,ppa
get_ipython().magic(u'run getiops.ipynb')


# In[4]:

Pnrm_psd = Pnrm_psd_straight.reshape(NANG,NIORS,Nreff)
print Pnrm_psd.shape

# print Pnrm_psd_straight[0:10,0] # first r-value, first wavelength
# print Pnrm_psd_straight[0:10,1] # second r-value, first wavelength
# print Pnrm_psd_straight[0:10,40] # first r-value, second wavelength
# print Pnrm_psd_straight[0:10,41] # second r-value, second wavelength

# print Pnrm_psd[0:10,0,0] # So this means 1st index is angle, second is wavelength, third is r
# print Pnrm_psd[0:10,0,1]
# print Pnrm_psd[0:10,1,0]
# print Pnrm_psd[0:10,1,1]


# In[5]:

iops_psd = iops_psd_straight.reshape(13,NIORS,Nreff)
print iops_psd.shape

# print iops_psd_straight[:,0] # first r-value, first wavelength
# print iops_psd_straight[:,1] # second r-value, first wavelength
# print iops_psd_straight[:,40] # first r-value, second wavelength
# print iops_psd_straight[:,41] # second r-value, second wavelength

# print iops_psd[:,0,0]
# print iops_psd[:,0,1]
# print iops_psd[:,1,0]
# print iops_psd[:,1,1]


# In[6]:

refflist = iops_psd[2,0,:]; print refflist
wnumlist = iops_psd[1,:,0]; #print wnumlist


# In[12]:

import pmomstuff3 as pm


# In[8]:

pmomarray, Npmomarray = pm.ssp2pmom(Pnrm_psd,refflist,wnumlist,ppa)


# In[9]:

wnum_mesh = np.squeeze(iops_psd.T[:,:,1]); #print "wnum", "\n", wnum_mesh[0,:], "\n", wnum_mesh[1,:], "\n"
reff_mesh = np.squeeze(iops_psd.T[:,:,2]); #print "reff", "\n", reff_mesh[:,0], "\n", reff_mesh[:,1], "\n"
w0_mesh   = np.squeeze(iops_psd.T[:,:,6]); #print "w0",   "\n", w0_mesh[:,0],   "\n", w0_mesh[:,1],   "\n"
qext_mesh = np.squeeze(iops_psd.T[:,:,8]); #print "qext", "\n", qext_mesh[:,0], "\n", qext_mesh[:,1], "\n"
reff_list = np.zeros((Nreff,1)); reff_list[:,0] = refflist; #print "reff", "\n", reff_list, reff_list.shape
wnum_list = np.zeros((NIORS,1)); wnum_list[:,0] = wnumlist; #print "wnum", "\n", wnum_list, wnum_list.shape


# In[10]:

# Save the moments as a netcdf file
pm.pmomsave(outfilename,Npmomarray,pmomarray,wnum_mesh,reff_mesh,w0_mesh,qext_mesh,wnum_list,reff_list)


# In[15]:

# # This shows how to load it
# Npmomarray_x, pmomarray_x, wnum_mesh_x, reff_mesh_x, w0_mesh_x, qext_mesh_x, wnum_list_x, reff_list_x = \
# pm.pmomload(outfilename)


# In[16]:

# # This shows how to graph it
# import matplotlib.pyplot as plt
# %matplotlib inline

# i_reff = 39
# i_wnum = 50

# N_fromfortran = Npmomarray[i_reff,i_wnum]; print "from fortran", N_fromfortran
# n_fromfortran = [i for i in range(N_fromfortran)] 
# p_fromfortran = np.squeeze(pmomarray[i_reff,i_wnum,n_fromfortran])

# plt.figure()
# plt.semilogy(n_fromfortran,p_fromfortran)
# plt.show()


# In[ ]:




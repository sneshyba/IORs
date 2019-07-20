
# coding: utf-8

# In[38]:

# Get resources
import numpy as np
import pmomstuff3 as pm
import scipy.io
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')


# In[39]:

# Get filenames
matNpmomfilename =    "../IORs/iors_water_273 (Hybrid)/ssp_getpsd_T273_S331_Npmom.mat"
matpmomfilename =    "../IORs/iors_water_273 (Hybrid)/ssp_getpsd_T273_S331_pmom.mat"
ncfilename = "../IORs/iors_water_273 (Hybrid)/ssp_getpsd_T273_S331_pmom.nc"


# In[40]:

# Get the netcdf data
Npmomarray, pmomarray, wnum_mesh, reff_mesh, w0_mesh, qext_mesh, wnum_list, reff_list = pm.pmomload(ncfilename)


# In[41]:

# Get the matlab-generated moments
matNpmom = scipy.io.loadmat(matNpmomfilename)
matpmom = scipy.io.loadmat(matpmomfilename)


# In[42]:

print matpmom['pmomarray'].shape
print matNpmom['Npmomarray'].shape
print pmomarray.shape
print Npmomarray.shape
print Npmomarray[30,:]


# In[47]:

i_reff = 5
i_wnum = 190

N_frommatlab = matNpmom['Npmomarray'][i_wnum,i_reff]; print "from matlab", N_frommatlab
n_frommatlab = [i for i in range(N_frommatlab)]
p_frommatlab = np.squeeze(matpmom['pmomarray'][n_frommatlab,i_wnum,i_reff])

N_fromfortran = Npmomarray[i_reff,i_wnum]; print "from fortran", N_fromfortran
n_fromfortran = [i for i in range(N_fromfortran)] 
p_fromfortran = np.squeeze(pmomarray[i_reff,i_wnum,n_fromfortran])

plt.figure()
plt.semilogy(n_frommatlab,p_frommatlab,'o',n_fromfortran,p_fromfortran)
plt.legend(['matlab', 'fortran'])
plt.show()


# In[ ]:




# pmomstuff4.py

import numpy

# convert ssp to pmom
def ssp2pmom(Pnrm_psd,refflist,reffmax, wnumlist,ppa):

  # S. Neshyba, Jan. 2012
  # sspfile, pmomfile are names of files
  # Nreff, if supplied, specifies an override for the maximum number of radii
  # Niors, if supplied, specifies an override for the maximum number of wavenumbers
  # Calls external functions getmompf, readsspfileblind2, ssp2pmom (as a backup)

  # Imports
  import moms

  
  # Make local copy / limit the ranges if desired
  Nang, Niors, Nreff = Pnrm_psd.shape
  #Nreff = len(refflist) # Overrides the result from Pnrm_psd (for debugging)
  print ("Nreff=", Nreff)
  print ("Niors=", Niors)  
  print ("Nang=", Nang)

  # Call the ssp2pmom code for comparison
  #pmom_matlab = mlab.ssp2pmom(ssp);

  # Parameters for getmompf
  thetadeg = 0.25
  maxtableg = 1000 # This parameter does not actually seem to be used by getmompf, but we'll use it here for the pmom array

  # Set aside some memory to hold results
  pmomarray = numpy.zeros((Nreff,Niors,maxtableg))
  Npmomarray = numpy.zeros((Nreff,Niors),dtype=int)
  

  # Loop over the radii
  for ireff in range(Nreff):
  #for ireff in range(2):

   # Extract this radius
   reff = refflist[ireff]

   # Reporting
   print ("reff = ", reff)
   
   if (reff < reffmax):

       # Loop over the wavenumbers
       for iiors in range(Niors):
    
        # Extract this wavenumber
        nu = wnumlist[iiors]
        wlen = 1e4/nu
           
        # calculate the size parameter for getmompf
        sizeparam = 2*numpy.pi*reff/wlen
    
        # Call getmompf to normalize, and again to get the legendre moments
        phasefun = numpy.squeeze(Pnrm_psd[:,iiors,ireff])
        legen, nleg, newphase = moms.getmompf(phasefun, ppa, Nang, maxtableg, sizeparam, thetadeg, 0, 1)
        #if (nleg > maxtableg): 
        #    print ("nleg = ", nleg)
        #[legen,nleg,newnewphase]=moms.getmompf(newphase, ssp.ppa, ssp.Nang, maxtableg, sizeparam, thetadeg, 0, 0)
        
        # Save these moments
        ncollected = min(nleg,maxtableg)
        pmomarray[ireff,iiors,0:ncollected] = legen[0:ncollected]
        Npmomarray[ireff,iiors] = ncollected

  # return it
  return pmomarray, Npmomarray

# just testing
def hithere():
	print ('hi there everybody')
	return


# Save the moments
#def pmomsave(datafile,Npmomarray,pmomarray,wnum_mesh,reff_mesh,w0_mesh,qext_mesh,wnum_list,reff_list):
#	from scipy.io.netcdf import NetCDFFile as DS
#	import numpy
#	(Nreff,Niors,Nmoms) = pmomarray.shape
#	(Nreff_list, dummy) = reff_list.shape
#	(Niors_list, dummy) = reff_list.shape
#	nc = DS(datafile,'w')
#	nc.createDimension('reff',Nreff)
#	nc.createDimension('iors',Niors)
#	nc.createDimension('moms',Nmoms)
#	nc.createDimension('one', 1)
#	data = nc.createVariable('Npmomarray',numpy.dtype('int32').char,  ('reff','iors')); 		data[:]=Npmomarray
#	data = nc.createVariable('pmomarray', numpy.dtype('float64').char,('reff','iors','moms')); data[:]=pmomarray
#	data = nc.createVariable('wnum_mesh', numpy.dtype('float64').char,('reff','iors')); 		data[:]= wnum_mesh[0:Nreff,0:Niors]
#	data = nc.createVariable('reff_mesh', numpy.dtype('float64').char,('reff','iors')); 		data[:]= reff_mesh[0:Nreff,0:Niors]
#	data = nc.createVariable('w0_mesh',   numpy.dtype('float64').char,('reff','iors')); 		data[:]=   w0_mesh[0:Nreff,0:Niors]
#	data = nc.createVariable('qext_mesh', numpy.dtype('float64').char,('reff','iors')); 		data[:]= qext_mesh[0:Nreff,0:Niors]
#	data = nc.createVariable('reff_list', numpy.dtype('float64').char,('reff','one')); 	 	data[:]= reff_list[0:Nreff,:]
#	data = nc.createVariable('wnum_list', numpy.dtype('float64').char,('iors','one')); 		data[:]= wnum_list[0:Niors,:]
#	nc.close()
#	return

def pmomsave2(\
                datafile, \
                Npmomarray, pmomarray, \
                wnum_mesh, reff_mesh, w0_mesh, qext_mesh, asym_mesh, \
                wnum_list, reff_list, maxdim_list, volume_list, parea_list):
            
	from scipy.io.netcdf import NetCDFFile as DS
	import numpy
	(Nreff,Niors,Nmoms) = pmomarray.shape
	nc = DS(datafile,'w')
	nc.createDimension('reff',Nreff)
	nc.createDimension('iors',Niors)
	nc.createDimension('moms',Nmoms)
	nc.createDimension('one', 1)
	data = nc.createVariable('Npmomarray',  numpy.dtype('int32').char,  ('reff','iors')); 		data[:]= Npmomarray
	data = nc.createVariable('pmomarray',   numpy.dtype('float64').char,('reff','iors','moms'));  data[:]= pmomarray
	data = nc.createVariable('wnum_mesh',   numpy.dtype('float64').char,('reff','iors')); 		data[:]= wnum_mesh[0:Nreff,0:Niors]
	data = nc.createVariable('reff_mesh',   numpy.dtype('float64').char,('reff','iors')); 		data[:]= reff_mesh[0:Nreff,0:Niors]
	data = nc.createVariable('w0_mesh',     numpy.dtype('float64').char,('reff','iors')); 		data[:]=   w0_mesh[0:Nreff,0:Niors]
	data = nc.createVariable('qext_mesh',   numpy.dtype('float64').char,('reff','iors')); 		data[:]= qext_mesh[0:Nreff,0:Niors]
	data = nc.createVariable('asym_mesh',   numpy.dtype('float64').char,('reff','iors')); 		data[:]= asym_mesh[0:Nreff,0:Niors]
	data = nc.createVariable('wnum_list',   numpy.dtype('float64').char,('iors','one')); 		data[:]= wnum_list[0:Niors,:]
	data = nc.createVariable('reff_list',   numpy.dtype('float64').char,('reff','one')); 	 	data[:]= reff_list[0:Nreff,:]
	data = nc.createVariable('maxdim_list', numpy.dtype('float64').char,('reff','one')); 	 	data[:]= maxdim_list[0:Nreff,:]
	data = nc.createVariable('volume_list', numpy.dtype('float64').char,('reff','one')); 	 	data[:]= volume_list[0:Nreff,:]
	data = nc.createVariable('parea_list',  numpy.dtype('float64').char,('reff','one')); 	 	data[:]= parea_list[0:Nreff,:]
	nc.close()
	return


# Load the moments
def pmomload2(datafile):
	#(Npmomarray,pmomarray,wnum_mesh,reff_mesh,w0_mesh,qext_mesh,wnum_list,reff_list) = pmomload(pmomfile)
	from scipy.io.netcdf import NetCDFFile as DS
	nc = DS(datafile,'r')
	Npmomarray = nc.variables['Npmomarray'][:]; Npmomarray = Npmomarray.astype('int32')
	pmomarray  = nc.variables['pmomarray'][:]; pmomarray = pmomarray.astype('float64')
	wnum_mesh  = nc.variables['wnum_mesh'][:]; wnum_mesh = wnum_mesh.astype('float64')
	reff_mesh  = nc.variables['reff_mesh'][:]; reff_mesh = reff_mesh.astype('float64')
	w0_mesh    = nc.variables['w0_mesh'][:];   w0_mesh   = w0_mesh.astype('float64')
	qext_mesh  = nc.variables['qext_mesh'][:]; qext_mesh = qext_mesh.astype('float64')
	asym_mesh  = nc.variables['asym_mesh'][:]; asym_mesh = asym_mesh.astype('float64')
	wnum_list  = nc.variables['wnum_list'][:]; wnum_list = wnum_list.astype('float64')
	reff_list  = nc.variables['reff_list'][:]; reff_list = reff_list.astype('float64')
	maxdim_list  = nc.variables['maxdim_list'][:]; maxdim_list = maxdim_list.astype('float64')
	volume_list  = nc.variables['volume_list'][:]; volume_list = volume_list.astype('float64')
	parea_list  = nc.variables['parea_list'][:]; parea_list = parea_list.astype('float64')
	nc.close()
	return  Npmomarray, pmomarray, \
             wnum_mesh, reff_mesh, w0_mesh, qext_mesh, asym_mesh, \
             wnum_list, reff_list, maxdim_list, volume_list, parea_list
 
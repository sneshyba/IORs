# pmomstuff.py

# convert ssp to pmom
def ssp2pmom(sspfile,pmomfile,Nreff=-1,Niors=-1):

  # S. Neshyba, Jan. 2012
  # sspfile, pmomfile are names of files
  # Nreff, if supplied, specifies an override for the maximum number of radii
  # Niors, if supplied, specifies an override for the maximum number of wavenumbers
  # Calls external functions getmompf, readsspfileblind2, ssp2pmom (as a backup)

  # Imports
  from mlabwrap import mlab
  import numpy
  import disort
  import moms

  # Get ssp data
  ssp = mlab.readsspfileblind2(sspfile)
  print ssp

  # Make local copy / limit the ranges if desired
  if Nreff<0:
	Nreff = int(round(ssp.Nreff))
  if Niors<0:
	Niors = int(round(ssp.Niors))
  Nang = int(round(ssp.Nang))
  print Nreff, Niors, Nang

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

   # Extract this radius
   reff = ssp.reff[ireff]

   # Loop over the wavenumbers
   for iiors in range(Niors):

    # Extract this wavenumber
    nu = ssp.wnum[iiors]
    wlen = 1e4/nu

    # calculate the size parameter for getmompf
    sizeparam = 2*numpy.pi*reff/wlen

    # Call getmompf to normalize, and again to get the legendre moments
    [legen,nleg,newphase]=moms.getmompf(ssp.Pnrm[:,ireff,iiors], ssp.ppa, ssp.Nang, maxtableg, sizeparam, thetadeg, 0, 1)
    [legen,nleg,newnewphase]=moms.getmompf(newphase, ssp.ppa, ssp.Nang, maxtableg, sizeparam, thetadeg, 0, 0)
    
    mlab.plotpmom(legen[0:nleg],reff,nu,Nang)
    #mlab.plot(ssp.mu,newphase,ssp.mu,ssp.Pnrm[:,ireff,iiors],'o')
    print reff, nu, nleg

    # Save these moments
    ncollected = min(nleg,maxtableg)
    pmomarray[ireff,iiors,0:ncollected] = legen[0:ncollected]
    Npmomarray[ireff,iiors] = ncollected

  # Save to file
  pmomsave(pmomfile,Npmomarray,pmomarray,ssp.wnum_mesh,ssp.reff_mesh,ssp.w0_mesh,ssp.qext_mesh,ssp.wnum,ssp.reff)


  return

# just testing
def hithere():
	print 'hi there'
	return


# Save the moments
def pmomsave(datafile,Npmomarray,pmomarray,wnum_mesh,reff_mesh,w0_mesh,qext_mesh,wnum_list,reff_list):
	from scipy.io.netcdf import NetCDFFile as DS
	import numpy
	(Nreff,Niors,Nmoms) = pmomarray.shape
	(Nreff_list, dummy) = reff_list.shape
	(Niors_list, dummy) = reff_list.shape
	nc = DS(datafile,'w')
	nc.createDimension('reff',Nreff)
	nc.createDimension('iors',Niors)
	nc.createDimension('moms',Nmoms)
	nc.createDimension('one', 1)
	data = nc.createVariable('Npmomarray',numpy.dtype('int32').char,  ('reff','iors')); 		data[:]=Npmomarray
	data = nc.createVariable('pmomarray', numpy.dtype('float64').char,('reff','iors','moms')); 	data[:]=pmomarray
	data = nc.createVariable('wnum_mesh', numpy.dtype('float64').char,('reff','iors')); 		data[:]= wnum_mesh[0:Nreff,0:Niors]
	data = nc.createVariable('reff_mesh', numpy.dtype('float64').char,('reff','iors')); 		data[:]= reff_mesh[0:Nreff,0:Niors]
	data = nc.createVariable('w0_mesh',   numpy.dtype('float64').char,('reff','iors')); 		data[:]=   w0_mesh[0:Nreff,0:Niors]
	data = nc.createVariable('qext_mesh', numpy.dtype('float64').char,('reff','iors')); 		data[:]= qext_mesh[0:Nreff,0:Niors]
	data = nc.createVariable('reff_list', numpy.dtype('float64').char,('reff','one')); 	 	data[:]= reff_list[0:Nreff,:]
	data = nc.createVariable('wnum_list', numpy.dtype('float64').char,('iors','one')); 		data[:]= wnum_list[0:Niors,:]
	nc.close()
	return

# Load the moments
def pmomload(datafile):
	#(Npmomarray,pmomarray,wnum_mesh,reff_mesh,w0_mesh,qext_mesh,wnum_list,reff_list) = pmomload(pmomfile)
	from scipy.io.netcdf import NetCDFFile as DS
	nc = DS(datafile,'r')
	Npmomarray = nc.variables['Npmomarray'][:]; Npmomarray = Npmomarray.astype('int32')
	pmomarray  = nc.variables['pmomarray'][:]; pmomarray = pmomarray.astype('float64')
	wnum_mesh  = nc.variables['wnum_mesh'][:]; wnum_mesh = wnum_mesh.astype('float64')
	reff_mesh  = nc.variables['reff_mesh'][:]; reff_mesh = reff_mesh.astype('float64')
	w0_mesh    = nc.variables['w0_mesh'][:]; w0_mesh = w0_mesh.astype('float64')
	qext_mesh  = nc.variables['qext_mesh'][:]; qext_mesh = qext_mesh.astype('float64')
	reff_list  = nc.variables['reff_list'][:]; reff_list = reff_list.astype('float64')
	wnum_list  = nc.variables['wnum_list'][:]; wnum_list = wnum_list.astype('float64')
	nc.close()
	return Npmomarray, pmomarray,wnum_mesh,reff_mesh,w0_mesh,qext_mesh,wnum_list,reff_list	

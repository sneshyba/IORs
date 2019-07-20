def unpack(sspdata):
    Pmom_mesh = sspdata.Pmom
    reff_mesh = sspdata.reff_mesh
    wnum_mesh = sspdata.wnum_mesh
    qext_mesh = sspdata.qext_mesh
    w0_mesh = sspdata.w0_mesh
    NMOM = sspdata.lmax
    return(Pmom_mesh,reff_mesh,wnum_mesh,qext_mesh,w0_mesh,NMOM)

def combine(DTAUC,PMOM,SSALB,cldyLYR,cldyTAU_vis,cldyqext,cldyPmom,cldyw0):
    import copy
    prevTAU = DTAUC[cldyLYR]
    prevSSALB = SSALB[cldyLYR]
    prevPmom = PMOM[:,cldyLYR]
    cldyTAU = cldyTAU_vis*cldyqext/2
    newTAU = cldyTAU + prevTAU 
    newSSALB = (prevSSALB*prevTAU + cldyw0*cldyTAU)/newTAU
    newPmom = (prevSSALB*prevTAU*prevPmom + cldyw0*cldyTAU*cldyPmom)/(newTAU*newSSALB)
    SSALB_out = copy.deepcopy(SSALB); SSALB_out[cldyLYR] = newSSALB
    DTAUC_out = copy.deepcopy(DTAUC); DTAUC_out[cldyLYR] = newTAU
    PMOM_out  = copy.deepcopy(PMOM);   PMOM_out[:,cldyLYR] = newPmom
    return(DTAUC_out,PMOM_out,SSALB_out)

def combine2(DTAUC_gas, SSALB_gas, NPmom_gas, Pmom_gas, DTAUC_wat, SSALB_wat, NPmom_wat, Pmom_wat):
	import numpy
	NLYR, Pmom_gas_dim = Pmom_gas.shape
	NLYR, Pmom_wat_dim = Pmom_wat.shape
	DTAUC_new = DTAUC_gas + DTAUC_wat
	SSALB_new = numpy.zeros((NLYR,1))
	NPmom_new = numpy.zeros((NLYR,1)).astype('int32')
	for iLYR in range(NLYR):
		SSALB_new[iLYR] = (SSALB_gas[iLYR]*DTAUC_gas[iLYR] + SSALB_wat[iLYR]*DTAUC_wat[iLYR])/DTAUC_new[iLYR]
	if (Pmom_wat_dim > Pmom_gas_dim): 
		Pmom_new = numpy.zeros((NLYR,Pmom_wat_dim))
		Pmom_new_dim = Pmom_wat_dim
	else:
		Pmom_new = numpy.zeros((NLYR,Pmom_gas_dim))
		Pmom_new_dim = Pmom_gas_dim
	for iLYR in range(NLYR):
		NPmom_new[iLYR] = 0
		denom = DTAUC_new[iLYR]*SSALB_new[iLYR]
		#print 'denom: ', denom.shape, denom
		if (denom > 0.):
			for iPmom in range(Pmom_new_dim):
				#print 'Got a non-zero denom at ', iLYR
				if (iPmom < NPmom_wat and iPmom < NPmom_gas):
					Pmom_new[iLYR,iPmom] = (SSALB_gas[iLYR]*DTAUC_gas[iLYR]*Pmom_gas[iLYR,iPmom] + SSALB_wat[iLYR]*DTAUC_wat[iLYR]*Pmom_wat[iLYR,iPmom])/denom
					NPmom_new[iLYR] = NPmom_new[iLYR] + 1
					#print 'Incrementing NPmom_new to ', NPmom_new[iLYR], ' at ', iLYR
				elif (iPmom < NPmom_wat):
					Pmom_new[iLYR,iPmom] = (SSALB_wat[iLYR]*DTAUC_wat[iLYR]*Pmom_wat[iLYR,iPmom])/denom
					NPmom_new[iLYR] = NPmom_new[iLYR] + 1
					#print 'Incrementing NPmom_new to ', NPmom_new[iLYR], ' at ', iLYR
				elif (iPmom < NPmom_gas):
					Pmom_new[iLYR,iPmom] = (SSALB_gas[iLYR]*DTAUC_gas[iLYR]*Pmom_gas[iLYR,iPmom])/denom
					NPmom_new[iLYR] = NPmom_new[iLYR] + 1
					#print 'Incrementing NPmom_new to ', NPmom_new[iLYR], ' at ', iLYR
				else:
					break
					#print 'ran out of moments for layer ', iLYR
	NPmom_new = max(NPmom_new)
	return(DTAUC_new, SSALB_new, NPmom_new, Pmom_new)

def sspinterp3(reff_mesh,wnum_mesh,qext_mesh,w0_mesh,NPmom_mesh,Pmom_mesh,reff_i,wnum_i):
    import numpy
    from scipy.interpolate import interp2d
    #print "Here we are in sspinterp3"
    
    # Set up empty output arrays
    M = reff_i.size;
    qext_i = numpy.zeros((M,)) 
    w0_i = numpy.zeros((M,)) 
    NPmom_i_fp = numpy.zeros((M,))
    Pmom_i = 0
    
    #Extract single columns from the 2d reff and wnum arrays
    reff_vec = reff_mesh[:,0]; #print(reff_vec); #reff_mesh(:,1)'; #size(reff_vec), reff_vec(1:5)
    wnum_vec = wnum_mesh[0];   #print(wnum_vec); #wnum_mesh(1,:);  #size(wnum_vec), wnum_vec(1:5)
    #print(reff_mesh)
    #print(wnum_mesh)
    #print(reff_i, wnum_i)

    # Interpolate qext 
    #print(qext_mesh)
    f = interp2d(reff_vec,wnum_vec,qext_mesh.T)
    for i in range(0,M):
        qext_i[i] = f(reff_i[i],wnum_i[i])
    #print (qext_i)
 
    # Interpolate w0
    #print(w0_mesh)
    f = interp2d(reff_vec,wnum_vec,w0_mesh.T)
    for i in range(0,M):
        w0_i[i] = f(reff_i[i],wnum_i[i]); 
    #print (w0_i)

    # Even get an interpolated number of moments!
    #print(NPmom_mesh)
    f = interp2d(reff_vec,wnum_vec,NPmom_mesh.T);
    for i in range(0,M):
        NPmom_i_fp[i] = f(reff_i[i],wnum_i[i])
    NPmom_i = numpy.rint(NPmom_i_fp).astype(int)
    #print(NPmom_i)
    NPmom_max = max(NPmom_i); #print NPmom_max

    # Loop over all the moments to do the same
    Pmom_i = numpy.zeros((M,NPmom_max));
    for j in range(0,NPmom_max):
        #print (Pmom_mesh[:,:,j])
        f = interp2d(reff_vec,wnum_vec,Pmom_mesh[:,:,j].T)
        for i in range(0,M):
            Pmom_i[i,j] = f(reff_i[i],wnum_i[i])
        #print(Pmom_i[:,j])
        
    # Done
    return(qext_i,w0_i,NPmom_i,Pmom_i)

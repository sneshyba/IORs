# get them all
everything: moms.so

# This is the legendre moments fortran library
moms.so: getmompf.f getmom.f
	f2py -c -m moms getmompf.f getmom.f  --fcompiler=gfortran


README file for acrg_tdmcmc

Files in this directory:

acrg_hbtdmcmc_uncorr.f90  - Fortran script for running uncorrelated transdimensional, hierarchical reversible jump mcmc
acrg_hbtdmcmc_evencorr.f90  - Fortran script for running temporally correlated transdimensional, hierarchical reversible jump mcmc. Assumes every data point is evenly spaced in time.
acrg_hbtdmcmc_corr.f90  - Fortran script for running temporally correlated transdimensional, hierarchical reversible jump mcmc. Suitable for any separation of data points.

tdmcmc_inputs.py - Python template script to create right inputs for template_tdmcmc.f90

run_tdmcmc.py - Python script that reads in all specified inputs and calls fortran subroutine, writes ouput file and returns dataset of outputs 
                   - Other than major changes you shouldn't have to touch this on a day-to-day basis

tdmcmc_post_process.py - Python script containing routines for plotting tdmcmmc output and generating country totals
post_process_template.py - Template file to show how you might call functions contained in tdmcmc_post_process.py


***** Templates and inputs should be copied to a local area for you to edit *****************
****** They exist for reference only, at some point you have to actually do some work yourself!!!!! *************************

In order for the run_tdmcmc.py script to work you need to have compiled both Fortran subroutines using f2py.
By default the .so files must be called:
tdmcmc_uncorr.so, tdmcmc_evencorr.so, tdmcmc_corr.so 

But you don't have to stick with these defaults

***************** Compiling uncorrelated version with f2py ****************************************

The file tdmcmc_uncorr.so is created with:

f2py -c -m tdmcmc_uncorr acrg_hbtdmcmc_uncorr.f90

This compiles with gfortran, for a single processor. 

For parallel tempering, to compile for multiple processors using OMP do:

f2py -c -m tdmcmc_uncorr_pt --f90flags='-fopenmp' -lgomp acrg_hbtdmcmc_uncorr.f90

To compile using intel compiler use: 

f2py -c -m tdmcmc_uncorr --fcompiler=intelem acrg_hbtdmcmc_uncorr.f90

If you get a segmentation fault for an OMP run, then you will need to set the memory allocation using:

export KMP_STACKSIZE=512m

Sometimes this doesn't seem to work either, in which case also try

export OMP_STACKSIZE=512m
export GOMP_STACKSIZE=512m

This appears to be foolproof, but not sure which of these commands is key. Stick them in your .bashrc so you don't have to type them every time.



***************** Compiling evencorr hierarchical version with f2py ****************************************

Exactly the same as before, just change the file names. So:

f2py -c -m tdmcmc_evencorr acrg_hbtdmcmc_evencorr.f90

***************** Compiling correlated hierarchical version with f2py ****************************************

f2py -L/usr/lib64 -llapack -c -m tdmcmc_corr_s acrg_hbtdmcmc_corr.f90

f2py -L/usr/lib64 -llapack -c -m tdmcmc_corr_pt --f90flags='-fopenmp' -lgomp acrg_hbtdmcmc_corr.f90

***************** Compiling correlated hierarchical version with ifort ****************************************

f2py -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -c -m tdmcmc_corr_pt --fcompiler=intelem --f90flags='-fast -openmp' -liomp5 acrg_hbtdmcmc_corr.f90

To get this to work on air you'll need to copy the following into a terminal, or put in your bashrc:

export LD_PRELOAD=/opt/intel/mkl/lib/intel64/libmkl_core.so:/opt/intel/mkl/lib/intel64/libmkl_sequential.so

For some reason it works on snowy without this. Not sure about non-Bristol servers...

Also make sure you've got the following line in your .bashrc to set up the correct intle environments:

source /opt/intel/bin/compilervars.sh intel64

This should be faster than gfortran, my tests have shown run time improvements of ~2x to 4x faster.

******************************************************************************
For Edinburgh:

f2py -L/geos/usrgeos/intel/Compiler/xe2016/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -c -m test_corr --fcompiler=intelem --f90flags='-openmp' -liomp5 acrg_hbtdmcmc_corr.f90

# 02/2022 openmp now should be -qopenmp
f2py -L/geos/usrgeos/intel/Compiler/xe2016/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -c -m test_corr --fcompiler=intelem --f90flags='-qopenmp' -liomp5 acrg_hbtdmcmc_corr.f90

********************************* Running from python *************************************************

import module_name

Where module_name is the name of module_name.so

When calling the trasdimensional mcmc function follow what's in the template file, but if using a different module name make sure you change that.

************** If using openMP then you have to run from ipython not Spyder *********************************




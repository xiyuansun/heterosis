Will Landau
will.landau@gmail.com
will-landau.com
Iowa State University
July 2013


========= LICENSE =======================

This is free software. You can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.


========= ABOUT =========================

This readme explains the software that fits the model given in 
doc/writeup/writeup.Rnw (which may be compiled with knitR to produce a 
more easily readable pdf, doc/writeup/writeup.pdf). This package, at 
minimum, should include the following directories and files.

data/
  test/
    smallData.txt
    smallGroup.txt

doc/
  README.txt
  writeup/
    joint.csv
    writeup.bib
    writeup.Rnw

include/
  cpu/
    Chain.h
    Config.h
    constants.h
    functions.h

  gpu/
    Chain.h
    Config.h
    constants.h
    functions.h

Makefile

R/
  gelman-factors.r

sh/
  gelman-factors.sh

src/
  cpu/
	allocChain.c
	alp.c
	c.c
	config.c
	d.c
	del.c
	dic.c
	eps.c
	eta.c
	freeChain.c
	getopts.c
	interimResults.c
	logLik.c
	main.c
	mcmc.c
	mu.c
	newChain.c
	phi.c
	phiAlpDel.c
	phiAlpDelJoint.c
	piAlp.c
	piDel.c
	printArrays.c
	printChain.c
	printConfig.c
	printHeaders.c
	rbeta.c
	readData.c
	readGrp.c
	resetChain.c
	rgamma.c
	rnormal.c
	runChain.c
	runiform.c
	sigAlp.c
	sigC.c
	sigDel.c
	sigPhi.c
	summarizeChain.c
	tau.c
	theAlp.c
	theDel.c
	thePhi.c
	updateDICprep.c

  gpu/
	allocChain.cu
	alp.cu
	c.cu
	chainDeviceToHost.cu
	config.cu
	d.cu
	del.cu
	dic.cu
	eps.cu
	eta.cu
	freeChain.cu
	getopts.cu
	interimResults.cu
	main.cu
	mcmc.cu
	newChain.cu
	phi.cu
	phiAlpDel.cu
	phiAlpDelJoint.cu
	piAlp.cu
	piDel.cu
	printArrays.cu
	printChain.cu
	printConfig.cu
	printHeaders.cu
	rbeta.cu
	readData.cu
	readGrp.cu
	resetChain.cu
	rgamma.cu
	rnormal.cu
	runChain.cu
	runiform.cu
	sigAlp.cu
	sigC.cu
	sigDel.cu
	sigPhi.cu
	summarizeChain.cu
	tau.cu
	theAlp.cu
	theDel.cu
	thePhi.cu
	updateDICprep.cu
    

========= REQUIREMENTS ==================

This package contains two versions: a standard C version
and a CUDA C, GPU-accelerated version that takes advantage 
of opportunities to parallelize the MCMC across genes. In
addition, a GNU Bash script, sh/gelman-factors.sh, and an 
R script R/gelman-factors.c, are included for computing 
Gelman potential scale reduction factors on output.

General requirements:
- A Linux/Unix or Mac operating system. Windows systems 
  must run this program from within a Linux/Unix-style 
  command line interface such as Cygwin (cygwin.com) or 
  MinGW (mingw.org/).
- GNU Bash version 4.1.2 or higher.
- The GNU C library version 2.12 or higher.

Regular C version requirements:
- The GNU Compiler Collection (gcc) version 4.2.1 or higher.

CUDA C version requirements:
- A Compute Unified Device Architecture (CUDA) graphics 
  device with compute capability 2.0 or higher.
- CUDA Production Release version 4.2 or above.

gelman-factors.sh and gelman-factors.r requriements:
- The R language version 2.15.3.
- The coda package in R, available for download from within R at 
  http://cran.us.r-project.org/.


========= SETUP =========================

Open a Linux/Unix-based command line interface program (Terminal
on Mac and Linux). Be sure that your current working directory
is the root directory of this package. Compile the regular C 
version with

$ make

or equivalently,

$ make cpu

Either command creates the following files and directories:

bin/
  mcmc
  
obj/
  cpu/
    allocChain.o
    alp.o
    ...

The binary, mcmc, is the main program. See USAGE below for 
usage details.

If you have the required CUDA tools listed in REQUIREMENTS, 
you can compile the gpu version with

$ make gpu

This command creates the following files and directories:

bin/
  gpu-mcmc
  
obj/
  cpu/
    allocChain.o
    alp.o
    ...

The binary, gpu-mcmc, is the GPU-accelerated version of the 
main program. See USAGE below for usage details.

Before using gelman-factors.sh and gelman-factors.r to compute 
Gelman factors on the output of the main program, you must install 
R version 2.15.3 or higher and the coda package. Download and 
install R from http://cran.us.r-project.org/. Then, to download 
the coda package, open R with

$ R

Begin the download and installation of the coda package with

> install.packages("coda")

and follow the directions. If you like (though not required), 
you can load the coda package into your current R session with

> library(coda)

and explore its contents with

> help(package = coda)

To read about the Gelman diagnostics function that gelman-factors.r 
uses, type

> ?gelman.diag


========= USAGE =========================

Open a Linux/Unix-based command line interface program (Terminal 
on Mac and Linux). Be sure that your current working directory
is the root directory of this package. After compilation, 
change the current working directory to the bin/ directory.

$ cd bin/

and run the regular C version with

$ ./mcmc [OPTION1] [OPTION2] ...

The usage of the GPU-accelerated version is analogous.

$ cd bin/
$ ./gpu-mcmc [OPTION1] [OPTION2] ...

For either version, type any of the options in OPTIONS
in place of [OPTION1] [OPTION2] ... For example,

$ ./mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt --hyper

or

$ ./gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt --hyper

To compute Gelman factors on the output of either mcmc or gpu-mcmc, 
change the current working directory to the sh/ directory and run

$ ./gelman-factors.sh [OUTPUT_DIR]

where [OUTPUT_DIR] is the root directory of the output of the main program.
Note: in order to compute the Gelman factors on the output, you must use
the --hyper or --parms options when you run mcmc or gpu-mcmc. See OPTIONS
for details.


========= OPTIONS =======================
   

    
    
    
    
    
    
    
    
    
    





 

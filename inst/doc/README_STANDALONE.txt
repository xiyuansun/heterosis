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

  This guide explains the C/CUDA C standalone program 
  behind the R package, heterosis. This program fits the model given in 
  inst/doc/method/method.Rnw (which may be compiled with knitR to produce 
  a more easily readable pdf, inst/doc/method/method.pdf). At minimum,
  you should have the following directories and files.

  data/
    test/
      smallData.txt
      smallGroup.txt

  inst/
    doc/
      README.txt
      method/
        joint.csv
        method.bib
        method.Rnw

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
  addition, an R script R/gelman-factors.r, is included for 
  computing Gelman potential scale reduction factors on output.

  General requirements:
  - A command line interface from which to launch the main 
    program. On Mac and Linux, use Terminal. Windows users 
    may consider Cygwin (cygwin.com) or MinGW (mingw.org/).
  - GNU Bash version 4.1.2 or higher.
  - The GNU C library version 2.12 or higher.

  Regular C version requirements:
  - The GNU Compiler Collection (gcc) version 4.2.1 or higher.

  CUDA C version requirements:
  - A Compute Unified Device Architecture (CUDA) graphics 
	device with compute capability 2.0 or higher.
  - CUDA Production Release version 5.0 or above.

  gelman-factors.r requirements:
  - The R language version 2.15.0.
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
  change the current working directory to the sh/ directory open R 
  (version 2.15.0 or later)

  $ R
  ...
  > 

  load the R script, gelman-factors.r

  > source("R/gelman-factors.r")  

  and then run

  > gelmanFactors([OUTPUT_DIR], [PARMS])

  where [OUTPUT_DIR] is the root directory of the output of the main program.
  [OUTPUT_DIR] defaults to out/ within the current working directory. Set [PARMS]
  to True to compute Gelman factors on all the parameters and to False to compute
  Gelman factors on the hyperparameters only. [PARMS] defaults to False because
  computing Gelman factors on all the parameters is extremely time-consuming.
  Note: in order to compute the Gelman factors on the output, you must use
  the --hyper or --parms options when you run mcmc or gpu-mcmc. Otherwise,
  there will be no output for gelmanFactors() to use. See OPTIONS for details.


========= OPTIONS =======================

  Technically, no options are required for execution. However, if none
  of the options, --dic, --hyper, --parms, --probs, --rates, --time, 
  is set, then the program will not produce output even though
  each of these six options is not required individually.

  --data [DATA_FILE]   (equivalent: -i [DATA_FILE], -I [DATA_FILE])

    Specifies the RNA-seq data file, [DATA_FILE], to run through 
    the main program. Must be a flat file with no header and no row names. 
    The sole content is a G X N matrix of integers, where G is the number 
    of genes and N is the number of samples. Entry (g, n) is the expression 
    of gene g in sample n. Entries within a row are delimited by spaces 
    and rows are delimited by the line feed character. 
    See data/test/smallData.txt for an example of an input data file with 
    8 genes and 4 samples. If this option is unset, the data file defaults
    to data/data.txt.
	
  --group [GROUP_FILE]    (equivalent: -g [GROUP_FILE], -G [GROUP_FILE])

    Specifies the group file, [GROUP_FILE] that corresponds to 
    the RNA-seq data file. The group file contains a single row that 
    specifies the assignment of RNA-seq samples to experimental groups. 
    For example, the example group file, data/test/smallGroup.txt
    (corresponding to data/test/smallData.txt) contains
  
	1 2 3 3
  
    Hence, in data/test/smallData.txt, the first sample (column) corresponds
    to group 1 (parent 1) the second sample corresponds to group 2 (offspring)
    and samples 3 and 4 correspond to group 3 (parent 3). To model heterosis,
    GROUP_FILE must specify exactly 3 treatment groups. If only 2 groups are 
    specified - for example,
  
	1 2 2 1
  
    then the program does not consider heterosis and only models
    differential expression between the two groups given. If this option is not set,
    then the program defaults to data/group.txt.
  
  --chains [NUM_CHAINS]    (equivalent: -c [NUM_CHAINS], -C [NUM_CHAINS])

    Number of MCMC chains to run. Defaults to 2. Note: for Gelman factors to be
    computed using sh/gelman-factors.sh, R/gelman-factors.r, or some other tool,
    there must be at least 2 chains.
  
  --iter [NUM_ITERATIONS]    (equivalent: -m [NUM_ITERATIONS], -M [NUM_ITERATIONS])

    Length of each chain (number of MCMC iterations). If unset, defaults to 10.

  --burnin [BURNIN]    (equivalent: -b [BURNIN], -B [BURNIN])

    Specify the number of iterations to ignore (in addition to initial parameters)
    when calculating gene-specific differential expression and heterosis probabilities, 
    acceptance rates for Metropolis steps, and the deviance information criterion 
    (DIC). If unset, defaults to half the chain length (number of iterations).

  --seed [SEED]    (equivalent: -s [SEED], -S [SEED])

    Set the seed for random number generation. If unset, defaults to 0.

  --joint    (equivalent: -j, -J)

    If set, each (phi_g, alpha_g, delta_g) triplet will be sampled jointly in a
    single Metropolis step. If unset, the phi_g's, alpha_g's, and delta_g's will
    all be sampled in individual Metropolis steps. See inst/doc/method/method.Rnw
    for details about the model.
  
  --verbose    (equivalent: -v, -V)

    If set, the progress of the MCMC will be printed to stdout. If unset, nothing
    will be printed to stdout.

  --out [OUTPUT_DIR]    (equivalent: -o [OUTPUT_DIR], -O [OUTPUT_DIR])

    Specifies the output directory, [OUTPUT_DIR]. All the program's output will
    be written to directories and files inside [OUTPUT_DIR]. If this directory 
    does not already exist, the program will create it. However, [OUTPUT_DIR]
    will not be created if there is no output to produce: i.e., if none of the 
    options, --dic, --hyper, --parms, --probs, --rates, --time, is set. 
    If --out is not set, [OUTPUT_DIR] defaults to out/.
  
  --dic

    If set, the program will output one estimate of the model's deviance information
    criterion (DIC) to dic.txt inside [OUTPUT_DIR]. Note: this will slow both 
    programs down and take away most of the performance advantage of the 
    GPU-accelerated version for all the iterations after burn-in.

  --hyper    (equivalent: -h, -H)

    If set, the program will output the hyperparameters, including initial values,
    for each chain: hyper-chain1.txt, hyper-chain2.txt, etc. inside [OUTPUT_DIR]. 

  --parms    (equivalent: -P)

    If set, the program will output the non-hyper parameters, including initial values,
    for each chain: parms-chain1.txt, parms-chain2.txt, etc. inside [OUTPUT_DIR]. 
    Note: this will slow both programs down and take away most of the
    performance advantage of the GPU-accelerated version.

  --probs    (equivalent: -p)

    If set, the program will output the gene-specific probabilities of 
    differential expression, high parent heterosis, low parent heterosis, 
    and mid parent heterosis for each chain: probs-chain1.txt, probs-chain2.txt,
    etc. inside [OUTPUT_DIR].

  --rates    (equivalent: -r, -R)

    If set, the program will print the acceptance rates of Metropolis steps 
    for each chain: rates-chain1.txt, rates-chain2.txt, etc. inside [OUTPUT_DIR].

  --time    (equivalent: -t, -T)

    If set, the program will print out the time (in milliseconds) that it 
    takes to sample each parameter for each iteration in each chain:
    time-chain1.txt, time-chain2.txt, etc. inside [OUTPUT_DIR].
  
  --phi-prior [CHOICE]

    Use this option to choose the default prior for the phi_g's. 
    The following options are currently available:

    --phi-prior 0: default choice, documented in section 2 of inst/doc/method.pdf.

  --alpha-prior [CHOICE]

    Use this option to choose the default prior for the alpha_g's. 
    The following options are currently available:

    --alpha-prior 0: default choice, documented in section 2 of inst/doc/method.pdf.

  --delta-prior [CHOICE]

    Use this option to choose the default prior for the delta_g's. 
    The following options are currently available:

    --delta-prior 0: default choice, documented in section 2 of inst/doc/method.pdf.

  --debug [DEBUG_LEVEL]

    Use this option to print out internal data structures during runtime.
    Choose among 3 levels of verbosity:

    1) Set [DEBUG_LEVEL] TO 1 (--debug 1) to print out only the Config object.
    2) Set [DEBUG_LEVEL] TO 2 (--debug 2) to print out the Config object
       at the beginning and the Chain object just before each chain is run.
    3) Set [DEBUG_LEVEL] TO 3 (--debug 3) to print out the Config object
       at the beginning and the Chain object just before each chain is run
       and again after each iteration of the current chain.

  The following options set initial constants. 
  See inst/doc/method/method.Rnw for details on the model.

  --sigma-c0 [VAL]        (default: 10)
  --d0 [VAL]              (default: 1000)
  --a-tau [VAL]           (default: 100)
  --a-alpha [VAL]         (default: 1)
  --a-delta [VAL]         (default: 1)
  --b-tau [VAL]           (default: 100)
  --b-alpha [VAL]         (default: 1)
  --b-delta [VAL]         (default: 1)
  --gamma-phi [VAL]       (default: 2)
  --gamma-alpha [VAL]     (default: 2)
  --gamma-delta [VAL]     (default: 2)
  --sigma-phi0 [VAL]      (default: 2)
  --sigma-alpha0 [VAL]    (default: 2)
  --sigma-delta0 [VAL]    (default: 2)


  The following options set hyperparameters. The hyperparameters
  set here will remain constant throughout the MCMC. The others 
  will be sampled from their prior distributions initially
  and then sampled from their full conditional distributions
  in Gibbs steps in the MCMC.

  --sigma-c [VAL]    
  --d [VAL]               (equivalent: -d)    
  --tau [VAL]    
  --theta-phi [VAL]    
  --theta-alpha [VAL]    
  --theta-delta [VAL]    
  --sigma-phi [VAL]    
  --sigma-alpha [VAL]    
  --sigma-delta [VAL]    
  --pi-alpha [VAL]    
  --pi-delta [VAL]  


========= TESTING =======================    

To test the code, a small test dataset, data/test/smallData.txt, is provided. Its group file, data/test/smallGroup.txt, is also included. The user may run either the cpu or the gpu version on this test data to verify that the software runs properly. After compilation, change to the bin/ directory with your Linux-style command line interface.

  $ cd bin/

Then, run the cpu version,

  $ ./mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt --iter 20 --chains 3 --hyper --parms --probs --rates --time --dic --verbose --out cpu-output 

or the gpu version,

  $ ./gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt --iter 20 --chains 3 --hyper --parms --probs --rates --time --dic --verbose --out gpu-output

For a test with minimal output, you may instead run

  $ ./mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt --hyper --out cpu-output
  $ ./gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt --hyper --out gpu-output

To compute Gelman diagnostics on the output in the newly created cpu-ouput/ directory, change from the bin/ directory to the sh/ directory.

  $ cd ../sh

Then, run

  $ ./gelman-factors.sh ../bin/cpu-output



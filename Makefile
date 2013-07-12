CC=gcc
NVCC=nvcc

CCFLAGS=-c -Wall -pedantic -Iinclude/cpu
NVCCFLAGS=-c -Iinclude/gpu -arch=sm_20

LDFLAGS=-lm

BINDIR=bin/
OBJDIR=obj/
SRCDIR=src/
OUTDIR=out/

CCOBJDIR=$(OBJDIR)cpu/
CCSRCDIR=$(SRCDIR)cpu/

NVCCOBJDIR=$(OBJDIR)gpu/
NVCCSRCDIR=$(SRCDIR)gpu/

CCDEP=printArrays config getopts printConfig freeConfig
CCDEP+=mySampleInt readGrp readData
CCDEP+=allocChain newChain printChain freeChain chains
CCDEP+=mu runiform rnormal rgamma rbeta
CCDEP+=c sigC eta d tau eps
CCDEP+=phi alp del phiAlpDelJoint phiAlpDel
CCDEP+=thePhi theAlp theDel
CCDEP+=sigPhi sigAlp sigDel
CCDEP+=piAlp piDel
CCDEP+=runChain oneChain summarizeChain
CCDEP+=printProbs printRates printHyper printParms
CCDEP+=main

NVCCDEP=printArrays config getopts printConfig freeConfig
NVCCDEP+=mySampleInt readGrp readData
NVCCDEP+=allocChain newChain chainDeviceToHost printChain freeChain chains
NVCCDEP+=runiform rnormal rgamma rbeta
NVCCDEP+=c sigC eta d tau eps 
NVCCDEP+=phiAlpDel phi alp del phiAlpDelJoint
NVCCDEP+=thePhi theAlp theDel
NVCCDEP+=sigPhi sigAlp sigDel
NVCCDEP+=piAlp piDel
NVCCDEP+=runChain oneChain summarizeChain
NVCCDEP+=printProbs printRates printHyper printParms
NVCCDEP+=main

CCOBJ=$(foreach name, $(CCDEP), $(CCOBJDIR)$(name).o)
NVCCOBJ=$(foreach name, $(NVCCDEP), $(NVCCOBJDIR)$(name).o)

all: cpu
	
cpu: $(CCOBJ) $(OUTDIR)
	mkdir -p $(BINDIR)
	$(CC) $(CCOBJ) -o $(BINDIR)mcmc $(LDFLAGS)

$(CCOBJDIR)%.o: $(CCSRCDIR)%.c
	mkdir -p $(CCOBJDIR)
	$(CC) $(CCFLAGS) $< -o $@ 
	
gpu: $(NVCCOBJ) $(OUTDIR)
	mkdir -p $(BINDIR)
	$(NVCC) $(NVCCOBJ) -o $(BINDIR)gpumcmc $(LDFLAGS)

$(NVCCOBJDIR)%.o: $(NVCCSRCDIR)%.cu
	mkdir -p $(NVCCOBJDIR)
	$(NVCC) $(NVCCFLAGS) $< -o $@ 

$(OUTDIR):
	mkdir -p $(OUTDIR)

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	rm -rf $(OUTDIR)
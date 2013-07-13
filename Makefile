CC=gcc
NVCC=nvcc

CCFLAGS=-c -Wall -pedantic -Iinclude/cpu
NVCCFLAGS=-c -Iinclude/gpu -arch=sm_20

LDFLAGS=-lm

BINDIR=bin/
OBJDIR=obj/
SRCDIR=src/
OUTDIR=out/

CPUBIN=$(BINDIR)mcmc
GPUBIN=$(BINDIR)gpu-mcmc

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

cpu: $(CPUBIN)
	
$(CPUBIN): $(CCOBJ) 
	$(CC) $(CCOBJ) $(LDFLAGS) -o $(CPUBIN) 

$(CCOBJDIR)%.o: $(CCSRCDIR)%.c cpudirs
	$(CC) $(CCFLAGS) $< -o $@ 

gpu: $(GPUBIN)
	
$(GPUBIN): $(NVCCOBJ) 
	$(NVCC) $(NVCCOBJ) $(LDFLAGS) -o $(GPUBIN) 

$(NVCCOBJDIR)%.o: $(NVCCSRCDIR)%.cu
	$(NVCC) $(NVCCFLAGS) $< -o $@ 

.INTERMEDIATE: cpudirs gpudirs dirs

cpudirs: dirs
	mkdir -p $(CCOBJDIR)

gpudirs: dirs
	mkdir -p $(NVCCOBJDIR)

dirs:
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	rm -rf $(OUTDIR)
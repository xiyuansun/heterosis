CC=gcc
NVCC=nvcc

CCFLAGS=-c -Wall -pedantic -Iinclude/cpu
NVCCFLAGS=-c -Iinclude/gpu -arch=sm_20

LDFLAGS= -lm

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

DEP+=getopts readGrp readData config  
DEP+=printArrays printConfig printChain printHeaders
DEP+=allocChain newChain runChain resetChain freeChain
DEP+=runiform rnormal rgamma rbeta
DEP+=phi alp del phiAlpDelJoint phiAlpDel
DEP+=c sigC eta d tau eps
DEP+=thePhi theAlp theDel
DEP+=sigPhi sigAlp sigDel
DEP+=piAlp piDel
DEP+=intermResults summarizeChain
DEP+=updateDICprep dic
DEP+=mcmc main

CCDEP=$(DEP) mu logLik
NVCCDEP=$(DEP) chainDeviceToHost 

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

$(NVCCOBJDIR)%.o: $(NVCCSRCDIR)%.cu gpudirs
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
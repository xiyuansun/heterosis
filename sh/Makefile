CC=gcc
CCFLAGS=-c -Wall -pedantic -I../include/cpu

NVCC=nvcc
NVCCFLAGS=-c -I../include/gpu -arch=sm_20

LDFLAGS=-lm

CCDEP=printArrays config getopts printConfig freeConfig
CCDEP+=mySampleInt readGrp readData
CCDEP+=allocChain newChain printChain freeChain
CCDEP+=mu runiform rnormal rgamma rbeta
CCDEP+=c sigC eta d tau eps
CCDEP+=phi alp del phiAlpDelJoint phiAlpDel
CCDEP+=thePhi theAlp theDel
CCDEP+=sigPhi sigAlp sigDel
CCDEP+=piAlp piDel
CCDEP+=runChain oneChain summarizeChain
CCDEP+=printProbs printRates printHyper printParms
CCDEP+=main

CCOBJ=$(foreach name, $(CCDEP), ../obj/cpu/$(name).o)

NVCCDEP=printArrays config getopts printConfig freeConfig
NVCCDEP+=mySampleInt readGrp readData
NVCCDEP+=allocChain newChain chainDeviceToHost printChain freeChain
NVCCDEP+=runiform rnormal rgamma rbeta
NVCCDEP+=c sigC eta d tau # eps
# NVCCDEP+=phi alp del phiAlpDelJoint phiAlpDel
NVCCDEP+=thePhi theAlp theDel
NVCCDEP+=sigPhi sigAlp sigDel
NVCCDEP+=piAlp piDel
NVCCDEP+=runChain oneChain summarizeChain
NVCCDEP+=printProbs printRates printHyper printParms
NVCCDEP+=main

NVCCOBJ=$(foreach name, $(NVCCDEP), ../obj/gpu/$(name).o)

all: cpu
	
cpu: $(CCOBJ)
	mkdir -p ../bin
	$(CC) $(CCOBJ) -o ../bin/mcmc $(LDFLAGS)

$(CCOBJ): ../src/cpu/%.c
	mkdir -p ../obj/cpu
	$(CC) $(CCFLAGS) $< -o $@ 

gpu: $(NVCCOBJ)
	mkdir -p ../bin
	$(NVCC) $(NVCCOBJ) -o ../bin/gpumcmc $(LDFLAGS)

$(NVCCOBJ): ../src/gpu/%.c
	mkdir -p ../obj/gpu
	$(NVCC) $(NVCCFLAGS) $< -o $@ 

clean:
	rm -rf ../obj
	rm -rf ../bin
	
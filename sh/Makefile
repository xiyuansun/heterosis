CC=gcc
CCFLAGS=-c -Wall -pedantic -I../include/cpu

NVCC=nvcc
NVCCFLAGS=-c -I../include/gpu -arch=sm_20

LDFLAGS=-lm

DEP=printArrays config getopts printConfig freeConfig
DEP+=mySampleInt readGrp readData
DEP+=allocChain newChain printChain freeChain
DEP+=mu runiform rnormal rgamma rbeta
DEP+=c sigC eta d tau # eps
#DEP+=phi alp del phiAlpDelJoint phiAlpDel
DEP+=thePhi theAlp theDel
DEP+=sigPhi sigAlp sigDel
DEP+=piAlp piDel
DEP+=runChain oneChain summarizeChain
DEP+=printProbs printRates printHyper printParms
DEP+=main

CCOBJ=$(foreach name, $(DEP), ../obj/cpu/$(name).o)
NVCCOBJ=$(foreach name, $(DEP), ../obj/gpu/$(name).o)

all: cpu
	
cpu: $(CCOBJ)
	mkdir -p ../bin
	$(CC) $(CCOBJ) -o ../bin/mcmc $(LDFLAGS)

../obj/cpu/%.o: ../src/cpu/%.c
	mkdir -p ../obj/cpu
	$(CC) $(CCFLAGS) $< -o $@ 

clean:
	rm -rf ../obj
	rm -rf ../bin

gpu: $(NVCCOBJ)
	mkdir -p ../bin
	$(NVCC) $(NVCCOBJ) -o ../bin/gpumcmc $(LDFLAGS)

../obj/gpu/%.o: ../src/gpu/%.c
	mkdir -p ../obj/gpu
	$(NVCC) $(NVCCFLAGS) $< -o $@ 

clean:
	rm -rf ../obj
	rm -rf ../bin
	
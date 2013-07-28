CC=gcc
CCFLAGS=-c -Wall -pedantic -Iinclude/cpu
LDFLAGS= -lm

BINDIR=bin/
OBJDIR=obj/
SRCDIR=src/

CPUBIN=$(BINDIR)mcmc

CCOBJDIR=$(OBJDIR)cpu/
CCSRCDIR=$(SRCDIR)cpu/

DEP+=getopts readGrp readData config  
DEP+=printArrays printConfig printChain printHeaders
DEP+=allocChain newChain runChain resetChain freeChain
DEP+=runiform rnormal rgamma rbeta
DEP+=phi alp del phiAlpDelJoint phiAlpDel
DEP+=c sigC eta d tau eps
DEP+=thePhi theAlp theDel
DEP+=sigPhi sigAlp sigDel
DEP+=piAlp piDel
DEP+=interimResults summarizeChain
DEP+=updateDICprep dic
DEP+=mcmc main
DEP+=mu logLik

CCOBJ=$(foreach name, $(DEP), $(CCOBJDIR)$(name).o)

all: cpu

cpu: $(CPUBIN)
	
$(CPUBIN): $(CCOBJ) 
	$(CC) $(CCOBJ) $(LDFLAGS) -o $(CPUBIN) 

$(CCOBJDIR)%.o: $(CCSRCDIR)%.c cpudirs
	$(CC) $(CCFLAGS) $< -o $@ 

.INTERMEDIATE: cpudirs dirs

cpudirs: dirs
	mkdir -p $(CCOBJDIR)

dirs:
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
CC=gcc
CCFLAGS=-c -Wall -pedantic -I../include/cpu -shared -fPIC
LDFLAGS=-lm
SOFLAGS=-lm -shared

SRCDIR=cpu/
OBJDIR=../obj/
BINDIR=../bin/
BIN=$(BINDIR)heterosis-mcmc
SO=heterosis.so

DEP=getopts readGrp readData config  
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

OBJ=$(foreach name, $(DEP), $(OBJDIR)$(name).o)

all: $(BIN) $(SO)
	
$(BIN): $(OBJ) 
	$(CC) $(OBJ) $(LDFLAGS) -o $(BIN) 
	cp $(BIN) ../inst/doc/heterosis-mcmc

$(SO): $(OBJ)
	$(CC) $(OBJ) $(SOFLAGS) -o $(SO) 

$(OBJDIR)%.o: $(SRCDIR)%.c dirs
	$(CC) $(CCFLAGS) $< -o $@ 

.INTERMEDIATE: dirs

dirs:
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	rm $(SO)
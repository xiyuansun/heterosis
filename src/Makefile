CC=gcc
CFLAGS=-Wall -ansi -pedantic -lm -I../include

CPUSRC=betaHost.c config.c freeConfig.c gammaHost.c getopts.c normalHost.c printConfig.c test.c  uniformHost.c
CPUOBJ=betaHost.o config.o freeConfig.o gammaHost.o getopts.o normalHost.o printConfig.o test.o  uniformHost.o

all: test

test: $(CPUOBJ)
	$(CC) $(CFLAGS) $? -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $? -o $@

clean:
	rm *.o
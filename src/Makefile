all: test

test: betaHost.o config.o freeConfig.o gammaHost.o getopts.o normalHost.o printConfig.o test.o uniformHost.o
	gcc betaHost.o config.o freeConfig.o gammaHost.o getopts.o normalHost.o printConfig.o test.o uniformHost.o -o test

betaHost.o: betaHost.c
	gcc -c betaHost.c -o betaHost.o

config.o: config.c
	gcc -c config.c -o config.o

freeConfig.o: freeConfig.c
	gcc -c freeConfig.c -o freeConfig.o

gammaHost.o: gammaHost.c
	gcc -c gammaHost.c -o gammaHost.o

getopts.o: getopts.c
	gcc -c getopts.c -o getopts.o

normalHost.o: normalHost.c
	gcc -c normalHost.c -o normalHost.o

printConfig.o: printConfig.c
	gcc -c printConfig.c -o printConfig.o

test.o: test.c
	gcc -c test.c -o test.o

uniformHost.o: uniformHost.c
	gcc -c uniformHost.c -o uniformHost.o

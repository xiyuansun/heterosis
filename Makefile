SRCMK=src/Makefile

all: $(SRCMK)
	make -C src

$(SRCMK):
	./configure

clean:
	rm -rf obj
	rm -rf bin
	rm -rf inst/bin
	rm -f  $(SRCMK)
CFLAGS = -Wall -O1 -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY:all
all: out.txt
	cat out.txt

out.txt: matrix
	./$< > out.txt

.PHONEY:clean
clean:
	$(RM) *.o matrix *.txt

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)

cc = gcc
CFLAGS = -O -std=gnu11
LDLIBS = -lm

default: out.txt 
	cat out.txt

out.txt: mfk
	./mfk>out.txt

mfk: mfk.o
	$(cc) -o mfk mfk.o $(LDLIBS)

mfk.o: mfk.c
	$(cc) $(CFLAGS) -c mfk.c


CC=gcc
CFLAGS=-O -std=gnu11
LDLIBS=-lm

default: out.txt
	cat ./$<

out.txt: test
	./$< > $@


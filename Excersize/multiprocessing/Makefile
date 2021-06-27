CC = gcc
CFLAGS = -O -Wall  -std=gnu11
CFLAGS += -fopenmp
CFLAGS += -pthread 
LDLIBS = -lm 

default: out.mt.txt out.op.txt
	cat ./out.mt.txt
	cat ./out.op.txt


out.mt.txt: MonteCarlo
	./$< > $@

out.op.txt: MonteCarlo_op
	./$< > $@

.PHONEY:clean
clean:
	$(RM) *.o test *.txt MonteCarlo MonteCarlo_op 



CFLAGS = -Wall -Werror -O -std=gnu11 #-DNDEBUG
LDLIBS = -lm 
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)




.PHONEY: default
default: out.num.txt 
	cat  out.num.txt


out.num.txt: num 
	./$< > $@


num: num.o num_funcs.o




.PHONEY: clean

clean:
	$(RM) out* num *.o *.png *.log

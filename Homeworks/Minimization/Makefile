CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: out.mini.txt Exc.B.png
	cat out.mini.txt


out.mini.txt: mini
	./$< > $@
	$(RM) out.psitest.txt 

out.bwfit.txt: mini


mini: mini.o mini_funcs.o 



Exc.B.png: data.txt out.bwfit.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "Energy [GeV]";\
	set ylabel "Cross section [arb. unit]";\
	set key top;\
	set xrange [100:160];\
	set yrange [-3:5];\
	plot "$<" using 1:2:3 w yerr title "data"\
	, "out.bwfit.txt" using 1:2 w l title "fit";\
	' | gnuplot 



.PHONEY: clean

clean:
	$(RM) out* mini *.o *.png *.log


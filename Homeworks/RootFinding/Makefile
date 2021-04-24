CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: out.rf.txt Psi.png
	cat out.rf.txt


out.rf.txt: rf
	./$< > $@

rf: rf.o rf_funcs.o lin_funcs.o ode_funcs.o



Psi.png: out.Psi.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "r (a0) ";\
	set ylabel "Psi_0";\
	set key top;\
	set xrange [0:8];\
	set yrange [0:20];\
	plot "$<" using 1:2 w l title "wavefunction";\
	' | gnuplot 

.PHONEY: clean

clean:
	$(RM) out* rf *.o *.png *.log


CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: out.rf.txt Psi_0_boundary.png Psi_good_boundary.png
	cat out.rf.txt


out.rf.txt: rf
	./$< > $@
	$(RM) out.psitest.txt ##The ode driver prints the wavefunction to a file, which it also does during the rootfinding.
				#the out.Psi.test.txt is that file, and as a result the data in it is useless (writing on top of prewritten txt files is not a good idea...)
out.Psi5.txt: rf

out.Psi8.txt: rf

out.Psi65.txt: rf


rf: rf.o rf_funcs.o lin_funcs.o ode_funcs.o



Psi_0_boundary.png: out.Psi5.txt out.Psi8.txt out.Psi65.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "r (a0) ";\
	set ylabel "Psi_0";\
	set key top;\
	set xrange [0:8];\
	set yrange [0:0.5];\
	f(x)=x*exp(-x);\
	plot "$<" using 1:2 w l title "wavefunction f(x)=0 for x=5"\
	, "out.Psi8.txt" using 1:2 w l title "wavefunction f(x)=0 for x=8"\
	, "out.Psi65.txt" using 1:2 w l title "wavefunction f(x)=0 for x=10"\
	, f(x) w l title "analytical solution";\
	' | gnuplot

Psi_good_boundary.png: out.Psi5_2.txt out.Psi8_2.txt out.Psi65_2.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "r (a0) ";\
	set ylabel "Psi_0";\
	set key top;\
	set xrange [0:8];\
	set yrange [0:0.5];\
	f(x)=x*exp(-x);\
	plot "$<" using 1:2 w l title "wavefunction f(x)=x*exp(-x) for x=5"\
	, "out.Psi8_2.txt" using 1:2 w l title "wavefunction f(x)=x*exp(-x) for x=8"\
	, "out.Psi65_2.txt" using 1:2 w l title "wavefunction f(x)=x*exp(-x) for x=10"\
	, f(x) w l title "analytical solution";\
	' | gnuplot




.PHONEY: clean

clean:
	$(RM) out* rf *.o *.png *.log


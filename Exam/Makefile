CFLAGS = -Wall -Werror -O -std=gnu11 #-DNDEBUG
LDLIBS = -lm
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)




.PHONEY: default
default: out.ode1.txt out.ode2.txt Harmonic_oscillator1.png Harmonic_oscillator2.png Harmonic_oscillator3.png


out.ode1.txt: ode
	./$< > $@



out.ode2.txt: ode



ode: ode.o ode_funcs.o


Harmonic_oscillator1.png: out.osc1.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "time";\
	set ylabel "extension";\
	set title "Harmonic oscillator1";\
	set key top;\
	set xrange [0:20];\
	set yrange [-3:3];\
	fr(x)=sin(x);\
	fi(x)=0;\
	plot "$<" using 1:3 w l title "Numerical real part"\
	, "$<" using 1:4 w l title "Numerical imaginary part"\
	, fr(x) w l title "Excact real part"\
	, fi(x) w l title "Excact imaginary part";\
	' | gnuplot



Harmonic_oscillator2.png: out.osc2.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "time";\
	set ylabel "extension";\
	set title "Harmonic oscillator2";\
	set key top;\
	set xrange [0:5];\
	set yrange [-1:75];\
	fi(x)=cosh(x);\
	fr(x)=0;\
	plot "$<" using 2:3 w l title "Numerical real part"\
	, "$<" using 2:4 w l title "Numerical imaginary part"\
	, fr(x) w l title "Excact real part"\
	, fi(x) w l title "Excact imaginary part";\
	' | gnuplot


Harmonic_oscillator3.png: out.osc3.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "time";\
	set ylabel "extension";\
	set title "Harmonic oscillator3";\
	set key top;\
	set xrange [0:5];\
	set yrange [-75:30];\
	fr(x)=cosh(x)*sin(x);\
	fi(x)=cos(x)*sinh(x);\
	plot "$<" using 2:3 w l title "Numerical real part"\
	, "$<" using 2:4 w l title "Numerical imaginary part"\
	, fr(x) w l title "Excact real part"\
	, fi(x) w l title "Excact imaginary part";\
	' | gnuplot


.PHONEY: clean
clean:
	$(RM) out* ode *.o *.png *.log

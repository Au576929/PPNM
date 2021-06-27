CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
LDLIBS = -lm 




.PHONEY: default
default: out.ode.txt out.ode.grav.txt  Pandemic.png grav1.png grav2.png grav3.png Harmonic_oscillator.png
	cat  out.ode.txt


out.ode.txt: ode 
	./$< > $@

out.ode.grav.txt: ode_c
	./$< > $@

ode: ode.o ode_funcs.o

ode_c: ode_c.o ode_funcs.o


Harmonic_oscillator.png: out.osc.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "time";\
	set ylabel "extension";\
	set key top;\
	set xrange [0:20];\
	set yrange [-10:10];\
	plot "$<" using 1:2 w l title "result from stepper";\
	' | gnuplot





Pandemic.png: out.Tc1.txt out.Tc2.txt out.Tc3.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "time (days) ";\
	set ylabel "Infected";\
	set key top;\
	set xrange [0:100];\
	set yrange [0:5500000];\
	plot "$<" using 1:3 w l title "Infected Tc 5"\
	, "out.Tc2.txt" using 1:3 w l title "Infected Tc 2.5"\
	, "out.Tc3.txt" using 1:3 w l title "Infected Tc 1.25"\
	, "out.Tc4.txt" using 1:3 w l title "Infected Tc 0.625"\
	, "out.Tc5.txt" using 1:3 w l title "Infected Tc 0.3125";\
	' | gnuplot 




grav1.png: out.gravity.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x-coordinate";\
	set ylabel "y-coordinate";\
	set key top;\
	set xrange [-2:2];\
	set yrange [-1:1];\
	plot "$<" using 2:3 w l title "particle 1"\
        ' | gnuplot 

grav2.png: out.gravity.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x-coordinate";\
	set ylabel "y-coordinate";\
	set key top;\
	set xrange [-2:2];\
	set yrange [-1:1];\
	plot "$<" using 6:7 w l title "particle 2";\
	' | gnuplot 

grav3.png: out.gravity.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x-coordinate";\
	set ylabel "y-coordinate";\
	set key top;\
	set xrange [-2:2];\
	set yrange [-1:1];\
	plot "$<" using 10:11 w l title "particle 3";\
	' | gnuplot 





.PHONEY: clean

clean:
	$(RM) out* ode *.o *.png *.log

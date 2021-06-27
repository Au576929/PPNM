CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: out.ls.txt Decay.png
 	


out.ls.txt: ls
	./$< > $@
	cat out.ls.txt

ls: ls.o ls_funcs.o


Decay.png: out.fit.txt out.data.txt Makefile
	echo '\
	set terminal png;\
	set title "Dacay of Radium";\
	set output "$@";\
	set xlabel "Time (days)";\
	set ylabel "Counts";\
	set key top;\
	plot "$<" index 0 using 1:2 w l title "best fit"\
	, "$<" index 1 using 1:2 w l title "fit + Delta c"\
	, "$<" index 2 using 1:2 w l title "fit - Delta c"\
	, "out.data.txt" using 1:2:3 with errorbars title "Data";\
	' | gnuplot 



gsl.timing.png: out.gsl.times.txt Makefile
	echo '\
	set terminal png;\
	set title "QR decomposition gsl algoritm";\
	set output "$@";\
	set xlabel "matrix size";\
	set ylabel "time sec";\
	set key top;\
	a=1;\
	f(x)=a*x**3;\
	fit f(x) "$<" via a;\
	plot "$<" using 1:2 w lp title "gsl algorithm timing"\
	, f(x) with lines title "const*N^3";\
	' | gnuplot


.PHONEY: clean

clean:
	$(RM) out* ls *.o *.png *.log

CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: out.lin.txt Makefile my.timing.png gsl.timing.png both.timing.png
	cat out.lin.txt


out.lin.txt: lin
	./$< > $@

lin: lin.o lin_funcs.o


both.timing.png: out.my.times.txt out.gsl.times.txt Makefile
	echo '\
	set terminal png;\
	set title "QR decomposition my algorithm";\
	set output "$@";\
	set xlabel "matrix size";\
	set ylabel "time sec";\
	set key top;\
	set yrange [0:1];\
	a=1;\
	f(x)=a*x**3;\
	fit f(x) "$<" via a;\
	b=1;\
	g(x)=b*x**3;\
	fit g(x) "out.gsl.times.txt" via b;\
	plot "$<" using 1:2 w lp title "my algorithm timing"\
	, f(x) w l title "const*N³"\
	, g(x) w l title "const*N³"\
	, "out.gsl.times.txt" w lp title "gsl algoritm timing";\
	' | gnuplot


my.timing.png: out.my.times.txt Makefile
	echo '\
	set terminal png;\
	set title "QR decomposition my algorithm";\
	set output "$@";\
	set xlabel "matrix size";\
	set ylabel "time sec";\
	set key top;\
	a=1;\
	f(x)=a*x**3;\
	fit f(x) "$<" via a;\
	plot "$<" using 1:2 w lp title "my algorithm timing"\
	, f(x) w l title "const*N^3";\
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
	$(RM) out* lin *.o *.png *.log

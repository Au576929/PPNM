CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: out.eig.txt States.png times.png 
	cat  out.eig.txt


out.eig.txt: eig
	./$< > $@

eig: eig.o eig_funcs.o

States.png: out.states.txt Makefile
	echo '\
	set terminal png;\
	set title "States of infinite well";\
	set output "$@";\
	set xlabel "x";\
	set ylabel "State";\
	set key bottom right;\
	Pi=3.1415926535897932348;\
	a1=1;\
	a2=1;\
	a3=1;\
	S1(x)=a1*sin(Pi*x);\
	S2(x)=a2*sin(Pi*2*x);\
	S3(x)=a3*sin(Pi*3*x);\
	fit S1(x) "$<" index 0 via a1;\
	fit S2(x) "$<" index 1 via a2;\
	fit S3(x) "$<" index 2 via a3;\
	plot "$<" index 0 using 1:2 w lp title "numerical state 1"\
	, "$<" index 1 using 1:2 w lp title "numerical state 2"\
	, "$<" index 2 using 1:2 w lp title "numerical state 3"\
	, S1(x) w l title "exact state 1 norm fitted"\
	, S2(x) w l title "exact state 2 norm fitted"\
	, S3(x) w l title "exact state 3 norm fitted";\
	' | gnuplot 

times.png: out.my.times1.txt out.my.times2.txt out.gsl.times.txt Makefile
	echo '\
	set terminal png;\
	set title "times";\
	set output "$@";\
	set xlabel "Matrix size";\
	set ylabel "time (sec)";\
	set key top left;\
	a1=1;\
	a2=1;\
	a3=1;\
	f1(x)=a1*x**3;\
	f2(x)=a2*x**3;\
	f3(x)=a3*x**3;\
	fit f1(x) "$<" via a1;\
	fit f2(x) "out.my.times2.txt" via a2;\
	fit f3(x) "out.gsl.times.txt" via a3;\
	plot "$<" using 1:2 w lp title "My times 1"\
	, "out.my.times2.txt" using 1:2 w lp title "My times 2"\
	, "out.gsl.times.txt" using 1:2 w lp title "gsltimes"\
	, f1(x) w l title "x^3 my times 1 fit"\
	, f2(x) w l title "x^3 my times 2 fit"\
	, f3(x) w l title "x^3 gsl times fit";\
	' | gnuplot 


.PHONEY: clean

clean:
	$(RM) out* eig *.o *.png *.log

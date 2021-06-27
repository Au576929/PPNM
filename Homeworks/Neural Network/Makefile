CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS += $(shell /usr/bin/gsl-config --libs)

.PHONEY: default
default: out.NN.txt interpolation.png derivative.png integral.png partc.png
	


out.NN.txt: NN NN_c
	./$< > $@
	./NN_c >> $@
	cat $@

NN: NN.o NN_funcs.o mini_funcs.o print_funcs.o

NN_c: NN_c.o NN_funcs_partc.o mini_funcs.o print_funcs.o NN_funcs.o num_funcs.o


out.harm.osc.txt: NN_c


interpolation.png: out.interp.txt out.tab.func.txt  Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set key top;\
	set xrange [0:10];\
	set yrange [-5:5];\
	plot "out.tab.func.txt" using 1:2 w p title "data"\
	, "out.interp.txt" using 1:2 w l title "interpolation";\
	' | gnuplot

derivative.png: out.interp.deriv.txt out.tab.func.deriv.txt  Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set key top;\
	set xrange [0:10];\
	set yrange [-5:5];\
	plot "out.tab.func.deriv.txt" using 1:2 w l title "excact derivative"\
	, "out.interp.deriv.txt" using 1:2 w l title "derivative from interpolation";\
	' | gnuplot


integral.png: out.interp.integ.txt out.tab.func.integ.txt  Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set key top;\
	set xrange [0:10];\
	set yrange [-5:5];\
	plot "out.tab.func.integ.txt" using 1:2 w l title "excact integral"\
	, "out.interp.integ.txt" using 1:2 w l title "integral from interpolation";\
	' | gnuplot

partc.png: out.harm.osc.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set key top;\
	set xrange [0:10];\
	set yrange [-5:5];\
	f(x)=sin(x);\
	plot "$<" using 1:2 w l title "neural network result"\
	, f(x) w l title "excact solution";\
	' | gnuplot




.PHONEY: clean

clean:
	$(RM) out* NN tab_func *.o *.png *.log


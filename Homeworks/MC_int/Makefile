CFLAGS = -Wall -Werror  -O -std=gnu11 #-DNDEBUG
LDLIBS = -lm 




.PHONEY: default
default: out.mc.txt err.png
	cat  "$<"


out.mc.txt: mc
	./$< > $@


mc: mc.o mc_funcs.o



err.png: out.err.txt Makefile 
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "Number of points";\
	set ylabel "Error";\
	set key top;\
	set format x "%.0e";\
	set format y "%.0e";\
	a=1; b=-0.5; f(x)=a*x**b;\
	c=1; d=-1; g(x)=c*x**d;\
	fit f(x) "$<" using 1:2 via a,b;\
	fit g(x) "$<" using 1:3 via c,d;\
	set log y;\
	plot "$<" using 1:2 w p title "quasi monte carlo"\
	, "$<" using 1:3 w p title "plain monte carlo"\
	, f(x) w l title "fit1"\
	, g(x) w l title "fit2";\
	' | gnuplot 




.PHONEY: clean

clean:
	$(RM) out* mc *.o *.png *.log

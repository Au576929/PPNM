CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs) -lm

.PHONEY: default
default: gamma.png erf.png lngamma.png


gamma.png:gamma.txt Makefile tab_gamma.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set yrange [-20:20];\
	set key top left;\
	plot "$<" using 1:2 w l title "approx"\
	,  "$<" using 1:3 w l title "math.h"\
	,  "$<" using 1:4 w l title "gsl"\
	,  "tab_gamma.txt" using 1:2 w p title "tabulated data"\
	' | gnuplot 

lngamma.png:lngamma.txt Makefile tab_lngamma.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set yrange [-20:20];\
	set key top left;\
	plot "$<" using 1:2 w l title "approx"\
	,  "$<" using 1:3 w l title "math.h"\
	,  "$<" using 1:4 w l title "gsl"\
	,  "tab_lngamma.txt" using 1:2 w p title "tabulated data"\
	' | gnuplot 





erf.png:erf.txt Makefile tab_erf.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set yrange [-20:20];\
	set key top left;\
	plot "$<" using 1:2 w l title "approx"\
	,  "$<" using 1:3 w l title "math.h"\
	,  "$<" using 1:4 w l title "gsl"\
	,  "tab_erf.txt" using 1:2 w p title "tabulated data"\
	' | gnuplot 

erf.txt: erf
	./$< > $@


gamma.txt: gamma
	./$< > $@

lngamma.txt: lngamma
	./$< > $@

#data.txt: main
#	./$< > $@

.PHONEY: clean
clean:
	rm *.o gamma erf lngamma gamma.txt erf.txt lngamma.txt *.png  

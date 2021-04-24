CC = gcc
CFLAGS = -O  -std=gnu11 
LDLIBS = -lm 

.PHONEY: default
default: out.txt exp.png ms.pdf


exp.png: out.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set yrange [0:150];\
	set key top left;\
	plot "$<" using 1:2 w l title "approx"\
	,"$<" using 1:3 w l title "math.h"\
	' | gnuplot 

out.txt: main
	./$< > $@


ms.pdf: Repport.tex exp.png
	pdflatex $<
	pdflatex $<


.PHONEY: clean
clean:
	$(RM) main *.txt *.o exp.png

.PHONEY: test
test:
	echo $(LDLIBS)
	echo $(CC)
	echo $(RM)


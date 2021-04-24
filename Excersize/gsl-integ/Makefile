
CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.Phoney: default
default: out.mygslinteg.txt  out.myBessel.png  out.myerf.png


out.myerf.png: out.Erf.txt Makefile
	graph \
	--output-format png --bitmap-size 1024x768 \
	--x-label "x" --y-label "y" \
	--top-label "Erf function" \
	--line-mode 1 $< \
	> $@

out.myBessel.png: out.Bessel.txt Makefile
	graph \
	--output-format png --bitmap-size 1024x768 \
	--x-label "x" --y-label "y" \
	--top-label "Bessel function" \
	--line-mode 1 $< \
	> $@


out.mygslinteg.txt: gsl-integ
	./$< > $@

out.Erf.txt: Erf
	./$< > $@
out.Bessel.txt: Bessel
	./$< > $@


.PHONEY:clean
clean:
	$(RM) out* *.o gslinteg Erf Bessel


.PHONEY:test
test:
	echo CFLAGS
	echo $(CFLAGS)
	echo $(LDLIBs)

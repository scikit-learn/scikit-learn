CC	= colorgcc
LD	= gcc

PYREX	= python2.4-pyrexc

PYTHONINC	= -I/usr/include/python2.4
NUMPYINC	= -I/home/david/local/lib/python2.4/site-packages/numpy/core/include
OPTIMS		= -O3 -funroll-all-loops -march=pentium4 -msse2
WARN		= -Wall -W -Winline -Wstrict-prototypes -Wmissing-prototypes \
		-Waggregate-return -Wcast-align -Wcast-qual -Wnested-externs \
		-Wshadow -Wbad-function-cast -Wwrite-strings

CFLAGS	= $(PYTHONINC) $(NUMPYINC) $(OPTIMS) $(WARN)

targets: c_gden.so _rawden.so

c_gden.so: c_gden.o
	$(LD) -shared -o $@ $< -Wl,-soname,$@

_rawden.so: pure_den.o
	$(LD) -shared -o $@ $< -Wl,-soname,$@

c_gden.o: c_gden.c
	$(CC) -c $(CFLAGS) -fPIC $<

pure_den.o: pure_den.c
	$(CC) -c $(CFLAGS) -fPIC $<

#c_gmm.so: c_gmm.o
#	$(LD) -shared -o $@ $< -Wl,-soname,$@
#
#c_gmm.o: c_gmm.c
#	$(CC) -c $(CFLAGS) -fPIC $<
#
#c_gmm.c: c_gmm.pyx c_numpy.pxd c_python.pxd
#	$(PYREX) $<

clean:
	rm -f c_gmm.c
	rm -f *.o
	rm -f *.so
	rm -f *.pyc

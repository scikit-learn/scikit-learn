CC	= colorgcc
LD	= gcc

PYREX	= python2.4-pyrexc

PYTHONINC	= -I/usr/include/python2.4
NUMPYINC	= -I/usr/lib/python2.4/site-packages/numpy/core/include
OPTIMS		= -O2 -ffast-math -march=pentium4
WARN		= -W -Wall

CFLAGS	= $(PYTHONINC) $(NUMPYINC) $(OPTIMS) $(WARN)

libgden.so: c_gden.o
	$(LD) -shared -o $@ $< -Wl,-soname,$@

c_gden.o: c_gden.c
	$(CC) -c $(CFLAGS) -fPIC $<

c_gmm.so: c_gmm.o
	$(LD) -shared -o $@ $< -Wl,-soname,$@

c_gmm.o: c_gmm.c
	$(CC) -c $(CFLAGS) -fPIC $<

c_gmm.c: c_gmm.pyx c_numpy.pxd c_python.pxd
	$(PYREX) $<

clean:
	rm -f c_gmm.c
	rm -f *.o
	rm -f *.so
	rm -f *.pyc

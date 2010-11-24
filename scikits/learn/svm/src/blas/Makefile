AR     = ar rcv
RANLIB = ranlib 

HEADERS = blas.h blas.h blasp.h
FILES = dnrm2.o daxpy.o ddot.o dscal.o 

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

blas: $(FILES) $(HEADERS)
	$(AR) blas.a $(FILES)  
	$(RANLIB) blas.a

clean:
	- rm -f *.o
	- rm -f *.a
	- rm -f *~

.c.o:
	$(CC) $(CFLAGS) -c $*.c



CC =gcc
ILINK =-I/usr/local/include
#-I/home/paul/gsl/include
LLINK =-L/usr/local/lib
#-L/home/paul/gsl/lib
OBJ =  main.o qrwrapper.o sim.o lattice.o random.o
DEPS = qrwrapper.h lattice.h sim.h random.h includes.h

%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) $(ILINK) -g  -c -o $@ $<

main : $(OBJ)
	gcc $(CFLAGS) -g -o $@ $^ -lgsl -lgslcblas -lm

.PHONY : clean
clean :
	rm -f $(OBJ) main

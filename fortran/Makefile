CC = gfortran
LDFLAGS = -L/opt/local/lib
LIBS = -lnetcdff -llapack

OBJS = read_netcdf.o
FILES = ensvsa-TE grads-ensvsa-TE plot-ensvsa-TE

.SUFFIXES: .f90 .o

.f90.o:
	$(CC) $(LDFLAGS) $(LIBS) -c $<

all:	$(FILES)

ensvsa-TE:	$(OBJS) ensvsa-TE.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^
grads-ensvsa-TE:	$(OBJS) grads-ensvsa-TE.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^
plot-ensvsa-TE:		$(OBJS) plot-ensvsa-TE.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^

ensvsa-TE.o:	ensvsa-TE.f90
grads-ensvsa-TE.o:	grads-ensvsa-TE.f90
plot-ensvsa-TE.o:	plot-ensvsa-TE.f90

clean:;		rm -f *.o ensvsa-TE grads-ensvsa-TE plot-ensvsa-TE

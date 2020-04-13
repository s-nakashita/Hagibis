CC = gfortran
LDFLAGS = -L/opt/local/lib
LIBS = -lnetcdff -llapack
OBJS = read_netcdf.f90

all:	ensvsa-TE grads-ensvsa-TE plot-ensvsa-TE

ensvsa-TE:	$(OBJS) ensvsa-TE.f90
	$(CC) $(OBJS) $(LDFLAGS) $(LIBS) ensvsa-TE.f90 -o ensvsa-TE.o
grads-ensvsa-TE:	$(OBJS) grads-ensvsa-TE.f90
	$(CC) $(OBJS) $(LDFLAGS) $(LIBS) grads-ensvsa-TE.f90 -o grads-ensvsa-TE.o
plot-ensvsa-TE:		$(OBJS) plot-ensvsa-TE.f90
	$(CC) $(OBJS) $(LDFLAGS) $(LIBS) plot-ensvsa-TE.f90 -o plot-ensvsa-TE.o

clean:;		rm -f *.o ensvsa-TE grads-ensvsa-TE plot-ensvsa-TE

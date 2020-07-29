CC = gfortran
LDFLAGS = -L/opt/local/lib
LIBS = -lnetcdff -llapack

OBJS = read_netcdf.o
FILES = ensvsa-TE grads-ensvsa-TE plot-ensvsa-TE

.SUFFIXES: .f90 .o

.f90.o:
	$(CC) $(LDFLAGS) $(LIBS) -c $<

all:	$(FILES)
dTE:		ensvsa-dTE grads-ensvsa-dTE
KE:		ensvsa-KE grads-ensvsa-KE
reg:	reg-prtb reg-mul-prtb

ensvsa-TE:	$(OBJS) ensvsa-TE.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^
grads-ensvsa-TE:	$(OBJS) grads-ensvsa-TE.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^
plot-ensvsa-TE:		$(OBJS) plot-ensvsa-TE.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^
ensvsa-dTE:	$(OBJS) ensvsa-dTE.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^
grads-ensvsa-dTE:	$(OBJS) grads-ensvsa-dTE.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^
ensvsa-KE:	$(OBJS) ensvsa-KE.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^
grads-ensvsa-KE:	$(OBJS) grads-ensvsa-KE.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^
reg-prtb:		$(OBJS) stat.o reg-prtb.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^
reg-mul-prtb:		$(OBJS) stat.o reg-mul-prtb.o
	$(CC) $(LDFLAGS) $(LIBS) -o $@ $^

$(OBJS):	read_netcdf.f90
ensvsa-TE.o:	ensvsa-TE.f90
grads-ensvsa-TE.o:	grads-ensvsa-TE.f90
plot-ensvsa-TE.o:	plot-ensvsa-TE.f90
ensvsa-dTE.o:	ensvsa-dTE.f90
grads-ensvsa-dTE.o:	grads-ensvsa-dTE.f90
ensvsa-KE.o:	ensvsa-KE.f90
grads-ensvsa-KE.o:	grads-ensvsa-KE.f90
stat.o:		stat.f90
reg-prtb.o:		reg-prtb.f90
reg-mul-prtb.o:		reg-mul-prtb.f90

clean:;		rm -f *.o *.mod ensvsa-TE grads-ensvsa-TE plot-ensvsa-TE ensvsa-KE grads-ensvsa-KE reg-prtb reg-mul-prtb
FC = gfortran
FFLAGS = -O2
CPP = gfortran
CPPFLAGS = -cpp -E -P#-Dmoist
LD = $(FC)
LDFLAGS = -L/opt/local/lib
LIBS = -lnetcdff -llapack

OBJS = read_netcdf.o
FILES = ensvsa grads-ensvsa tevol-ensvsa #plot-ensvsa-TE

PARENT_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
WKD = $(PARENT_DIR)tmp
BINDIR = $(PARENT_DIR)bin

.SUFFIXES:
.SUFFIXES: .f90 .o .F90

.f90.o:
	$(FC) $(FFLAGS) -c $<

.F90.o:
	$(CPP) $(CPPFLAGS) $< > CPPOUT/`basename $*.f90`
	$(FC) $(FFLAGS) -c CPPOUT/`basename $*.f90`

.o:
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

ensvsa:	$(OBJS) ensvsa.o
grads-ensvsa:	$(OBJS) grads-ensvsa.o
tevol-ensvsa:	$(OBJS) tevol-ensvsa.o

ensvsa-TE:	$(OBJS) ensvsa-TE.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
grads-ensvsa-TE:	$(OBJS) grads-ensvsa-TE.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
plot-ensvsa-TE:		$(OBJS) plot-ensvsa-TE.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
filter-ensvsa-TE:	filter-ensvsa-TE.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
tevol-ensvsa-TE:	$(OBJS) tevol-ensvsa-TE.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
ensvsa-dTE:	$(OBJS) ensvsa-dTE.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
grads-ensvsa-dTE:	$(OBJS) grads-ensvsa-dTE.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
ensvsa-KE:	$(OBJS) ensvsa-KE.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
grads-ensvsa-KE:	$(OBJS) grads-ensvsa-KE.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
reg-prtb:		$(OBJS) stat.o reg-prtb.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^
reg-mul-prtb:		$(OBJS) stat.o reg-mul-prtb.o
#	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

#read_netcdf.o:	read_netcdf.f90
ensvsa.o:	$(OBJS) #ensvsa.F90
grads-ensvsa.o:	$(OBJS) #grads-ensvsa.F90
tevol-ensvsa.o:	$(OBJS) #tevol-ensvsa.F90
ensvsa-TE.o:	$(OBJS) ensvsa-TE.F90
grads-ensvsa-TE.o:	$(OBJS) grads-ensvsa-TE.F90
plot-ensvsa-TE.o:	$(OBJS) plot-ensvsa-TE.f90
filter-ensvsa-TE.o:	filter-ensvsa-TE.f90
tevol-ensvsa-TE.o:	$(OBJS) tevol-ensvsa-TE.f90
ensvsa-dTE.o:	$(OBJS) ensvsa-dTE.f90
grads-ensvsa-dTE.o:	$(OBJS) grads-ensvsa-dTE.f90
ensvsa-KE.o:	$(OBJS) ensvsa-KE.f90
grads-ensvsa-KE.o:	$(OBJS) grads-ensvsa-KE.f90
stat.o:		stat.f90
reg-prtb.o:		reg-prtb.f90
reg-mul-prtb.o:		reg-mul-prtb.f90

#all:	$(FILES)
TE:	#	ensvsa-TE grads-ensvsa-TE tevol-ensvsa-TE
	rm -r $(WKD) || echo
	mkdir -p $(WKD) || echo
	mkdir -p $(WKD)/CPPOUT || echo
	mkdir -p $(BINDIR) || echo
	cd $(WKD);\
	$(MAKE) ensvsa -f $(PARENT_DIR)Makefile CPPFLAGS="$(CPPFLAGS) -Dmoist" VPATH=$(PARENT_DIR);\
	$(MAKE) grads-ensvsa -f $(PARENT_DIR)Makefile CPPFLAGS="$(CPPFLAGS) -Dmoist" VPATH=$(PARENT_DIR);\
	$(MAKE) tevol-ensvsa -f $(PARENT_DIR)Makefile CPPFLAGS="$(CPPFLAGS) -Dmoist" VPATH=$(PARENT_DIR);\
	mv ensvsa grads-ensvsa tevol-ensvsa $(BINDIR);\
	ls -l $(BINDIR)
	@echo END
dTE:	#ensvsa-dTE grads-ensvsa-dTE
	rm -r $(WKD) || echo
	mkdir -p $(WKD) || echo
	mkdir -p $(WKD)/CPPOUT || echo
	mkdir -p $(BINDIR) || echo
	cd $(WKD);\
	$(MAKE) ensvsa -f $(PARENT_DIR)Makefile VPATH=$(PARENT_DIR);\
	$(MAKE) grads-ensvsa -f $(PARENT_DIR)Makefile VPATH=$(PARENT_DIR);\
	$(MAKE) tevol-ensvsa -f $(PARENT_DIR)Makefile VPATH=$(PARENT_DIR);\
	mv ensvsa grads-ensvsa tevol-ensvsa $(BINDIR);\
	ls -l $(BINDIR)
	@echo END
KE:		ensvsa-KE grads-ensvsa-KE
reg:	reg-prtb reg-mul-prtb

clean:
	rm -f *.o *.mod ensvsa-TE grads-ensvsa-TE plot-ensvsa-TE tevol-ensvsa-TE filter-ensvsa-TE ensvsa-KE grads-ensvsa-KE reg-prtb reg-mul-prtb
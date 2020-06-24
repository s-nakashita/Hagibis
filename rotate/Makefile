CC = pgf90
LDFLAGS = -L/opt/local/lib 
LIBS = -lnetcdff 
SPLIBS = -l spherepack
OBJS = read_netcdf.o
FILES = test_sp

.SUFFIXES: .f90 .o

.f90.o:
	$(CC) -c $< $(LDFLAGS) $(LIBS) $(SPLIBS)

all:	$(FILES)

test_sp:	$(OBJS) test_sp.o
	$(CC) -o $@ $^ $(LDFLAGS) $(LIBS) $(SPLIBS)

$(OBJS):	read_netcdf.f90
test_sp.o:	test_sp.f90

clean:;		rm -f *.o *.mod test_sp
# Use this Makefile with make

# Executable name
CMD = dgradq.exe
# -------- description of DFLAGS ---------------


# -------- Define environmental variable C_COMPILER -----------
# Make sure it is defined
 FC = icpc -openmp -I$(MKLROOT)/include
OFLAGS =  # optimization
F95ROOT = $(MKLROOT)
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)
DFLAGS = -DVERSION=\"$(GIT_VERSION)\"


#Intel Linkers
#LINKERFLAGS =  -L$(MKLROOT)/lib/em64t $(F95ROOT)/lib/em64t/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#Intel parallel openmp (only w/icpc compiler)
#LINKERFLAGS =  -L$(MKLROOT)/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
LINKERFLAGS =  -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
# MAC OS linkers
#LINKERFLAGS = -lm -framework Accelerate


OFLAGS =  # optimization



#
# Implicit rules for handling src files
#  ( uses static pattern rules, see info for make )
.c.o:
	$(FC) -c -g $(DFLAGS) -Wimplicit $<
.cpp.o:
	$(FC) -c -g $(DFLAGS) $<

OBJECTS = main.o pTable.o stringtools.o mem.o icoord.o print.o

$(CMD) : $(OBJECTS)
	$(FC) $(DEBUG_FLAGS) $(OFLAGS) $(OBJECTS) $(LINKERFLAGS)   -o ./$(CMD)

clean:
	/bin/rm -f *.o *.i *.mod *.exe a.out make.log

cleano:
	rm -f *.o *.i

depend :
	g++ -MM *.cpp *.c >> Makefile 

# DO NOT DELETE created with g++ -MM *.cpp *.c
mem.o: mem.cpp icoord.h 
print.o: print.cpp icoord.h
icoord.o: icoord.cpp icoord.h stringtools.h pTable.h 
main.o: main.cpp icoord.h 
pTable.o: pTable.cpp pTable.h
stringtools.o: stringtools.cpp stringtools.h

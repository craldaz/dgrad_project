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

OBJECTS = main.o pTable.o stringtools.o mem.o icoord.o print.o grad.o bmat.o utils.o

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
grad.o: grad.cpp grad.h
bmat.o: bmat.cpp icoord.h 
icoord.o: icoord.cpp icoord.h stringtools.h pTable.h utils.h
main.o: main.cpp icoord.h 
pTable.o: pTable.cpp pTable.h
stringtools.o: stringtools.cpp stringtools.h
utils.o: utils.cpp utils.h constants.h
bmat.o: bmat.cpp icoord.h stringtools.h pTable.h grad.h utils.h \
 constants.h /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_types.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_blas.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_trans.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_cblas.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_spblas.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_lapack.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_lapacke.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_solver.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_dss.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_pardiso.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_rci.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_service.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vml.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vml_defines.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vml_types.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vml_functions.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vsl.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vsl_defines.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vsl_functions.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vsl_types.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_df.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_df_defines.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_df_functions.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_df_types.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_dfti.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_trig_transforms.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_poisson.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_solvers_ee.h
grad.o: grad.cpp grad.h stringtools.h pTable.h
icoord.o: icoord.cpp icoord.h stringtools.h pTable.h grad.h utils.h \
 constants.h
main.o: main.cpp icoord.h stringtools.h pTable.h grad.h
mem.o: mem.cpp icoord.h stringtools.h pTable.h grad.h
print.o: print.cpp icoord.h stringtools.h pTable.h grad.h constants.h
pTable.o: pTable.cpp pTable.h
stringtools.o: stringtools.cpp stringtools.h
utils.o: utils.cpp icoord.h stringtools.h pTable.h grad.h utils.h \
 constants.h /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_types.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_blas.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_trans.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_cblas.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_spblas.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_lapack.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_lapacke.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_solver.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_dss.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_pardiso.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_rci.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_service.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vml.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vml_defines.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vml_types.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vml_functions.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vsl.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vsl_defines.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vsl_functions.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_vsl_types.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_df.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_df_defines.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_df_functions.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_df_types.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_dfti.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_trig_transforms.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_poisson.h \
 /export/apps/Intel/composer_xe_2013.4.183/mkl/include/mkl_solvers_ee.h

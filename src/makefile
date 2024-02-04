DEBUG ?= 1

ifeq ($(DEBUG), 1)
   CFLAGS=-Wall -g
   CxxFLAGS=-Wall -std=c++17 -g -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(SLEPC_DIR)/$(PETSC_ARCH)/lib
else
   CFLAGS=-Wall -O3 
   CxxFLAGS=-Wall -std=c++17 -O3 
endif


CCxx=mpicxx
CxxINCS=-I../../OpenParEMCommon/src -I$(MFEM_DIR) -I$(MFEM_DIR)/linalg -I$(HYPRE_DIR)/src/hypre/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(EIGEN_DIR)
CxxLDIR=-L../../OpenParEMCommon/src -L$(MFEM_DIR) -L$(HYPRE_DIR)/src/hypre/lib -L$(METIS_DIR) -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/9
CxxLIBS=-lOpenParEMCommon -lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lflapack -lfblas -lptesmumps -lptscotchparmetisv3 -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lm -lX11 -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lmfem -lHYPRE -lmetis -lrt -lslepcpep -lslepcsys -lslepceps -lslepclme -lslepcnep -lslepcmfn -lslepcsvd /usr/lib/x86_64-linux-gnu/liblapacke64.a /usr/lib/x86_64-linux-gnu/liblapack64.a -lpetscmat -lpetscsnes -lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lflapack -lfblas -lslepcpep -lslepcsys -lslepceps -lslepclme -lslepcnep -lslepcmfn -lslepcsvd -lgfortran -lc  -lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp

CC=mpicc
CINCS=-I../../OpenParEMCommon/src -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(EIGEN_DIR)
CLDIR=-L../../OpenParEMCommon/src 
CLIBS=-lm 


builder.o: builder.cpp builder.hpp 
	$(CCxx) $(CxxFLAGS) -c builder.cpp $(CxxINCS)

convergence.o: convergence.cpp convergence.hpp project.h
	$(CCxx) $(CxxFLAGS) -c convergence.cpp $(CxxINCS)

eigensolve.o: eigensolve.c eigensolve.h project.h Hsolve.h
	$(CC) $(CFLAGS) -c eigensolve.c $(CINCS)

fem2D.o: fem2D.cpp fem2D.hpp modes.hpp project.h convergence.hpp fieldPoints.hpp results.hpp Hsolve.h
	$(CCxx) $(CxxFLAGS) -c fem2D.cpp $(CxxINCS)

fieldPoints.o: fieldPoints.cpp fieldPoints.hpp 
	$(CCxx) $(CxxFLAGS) -c fieldPoints.cpp $(CxxINCS)

findPoints.o: findPoints.cpp 
	$(CCxx) $(CxxFLAGS) -c findPoints.cpp $(CxxINCS)

Hsolve.o: Hsolve.c Hsolve.h
	$(CC) $(CFLAGS) -c Hsolve.c $(CINCS)

modes.o: modes.cpp modes.hpp
	$(CCxx) $(CxxFLAGS) -c modes.cpp $(CxxINCS)

OpenParEM2D.o: OpenParEM2D.cpp convergence.hpp fem2D.hpp fieldPoints.hpp modes.hpp project.h results.hpp 
	$(CCxx) $(CxxFLAGS) -c OpenParEM2D.cpp $(CxxINCS)

process.o: process.cpp process.hpp project.h
	$(CCxx) $(CxxFLAGS) -c process.cpp $(CxxINCS)

project.o: project.c project.h
	$(CC) $(CFLAGS) -c project.c $(CLDIR) $(CINCS)

results.o: results.cpp results.hpp
	$(CCxx) $(CxxFLAGS) -c results.cpp $(CxxINCS)

simplify.o: simplify.cpp simplify.hpp project.h
	$(CCxx) $(CxxFLAGS) -c simplify.cpp $(CxxINCS)

waveguide: waveguide.o
	$(CC) $(CFLAGS) waveguide.o -o waveguide $(CLDIR) $(CLIBS)

Zsolve.o: Zsolve.c Zsolve.h
	$(CC) $(CFLAGS) -c Zsolve.c $(CINCS)

simplify: simplify.o project.o
	$(CCxx) $(CxxFLAGS) simplify.o project.o -o simplify $(CxxLDIR) $(CxxLIBS)

process: process.o project.o
	$(CCxx) $(CxxFLAGS) process.o project.o -o process $(CxxLDIR) $(CxxLIBS)

builder: builder.o 
	$(CCxx) $(CxxFLAGS) builder.o -o builder $(CxxLDIR) $(CxxLIBS)

OpenParEM2D: OpenParEM2D.o convergence.o eigensolve.o fem2D.o fieldPoints.o findPoints.o Hsolve.o modes.o project.o results.o Zsolve.o
	$(CCxx) $(CxxFLAGS) OpenParEM2D.o convergence.o eigensolve.o fem2D.o fieldPoints.o findPoints.o Hsolve.o modes.o project.o results.o Zsolve.o -o OpenParEM2D  $(CxxLDIR) $(CxxLIBS)


.PHONY: all clean

all: OpenParEM2D waveguide simplify process builder 

install:
	rm -f ../bin/OpenParEM2D
	rm -f ../bin/process
	rm -f ../bin/simplify
	rm -f ../bin/builder
	rm -f $(HOME)/bin/OpenParEM2D
	rm -f $(HOME)/bin/process
	rm -f $(HOME)/bin/simplify
	rm -f $(HOME)/bin/builder
	rm -f $(HOME)/bin/process.sh
	rm -f $(HOME)/bin/proj_search
	rm -f $(HOME)/bin/regression.sh
	cp OpenParEM2D $(HOME)/bin
	cp process $(HOME)/bin
	cp simplify $(HOME)/bin
	cp builder $(HOME)/bin
	mv OpenParEM2D ../bin
	mv process ../bin
	mv simplify ../bin
	mv builder ../bin
	cp ../scripts/process.sh $(HOME)/bin/process.sh
	cp ../scripts/proj_search $(HOME)/bin/proj_search
	cp ../scripts/regression.sh $(HOME)/bin/regression.sh
	cp ../scripts/checkErrors.sh $(HOME)/bin/checkErrors.sh
	cp ../scripts/listErrors $(HOME)/bin/listErrors
	cp ../scripts/checkBuilder.sh $(HOME)/bin/checkBuilder.sh
	cp ../scripts/OpenParEM2D_save.py $(HOME)/.FreeCAD/Macro/OpenParEM2D_save.py

clean:
	rm -f OpenParEM2D
	rm -f process
	rm -f builder
	rm -f simplify
	rm -f waveguide
	rm -f builder.o
	rm -f convergence.o
	rm -f eigensolve.o
	rm -f fem2D.o
	rm -f fieldPoints.o
	rm -f findPoints.o
	rm -f Hsolve.o
	rm -f modes.o
	rm -f OpenParEM2D.o
	rm -f process.o
	rm -f project.o
	rm -f results.o
	rm -f simplify.o
	rm -f waveguide.o
	rm -f Zsolve.o

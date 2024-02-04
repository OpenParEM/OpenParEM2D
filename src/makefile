DEBUG ?= 0

ifeq ($(DEBUG), 1)
   CFLAGS=-Wall -g
   CxxFLAGS=-Wall -std=c++17 -g -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(SLEPC_DIR)/$(PETSC_ARCH)/lib
else
   CFLAGS=-Wall -O3 
   CxxFLAGS=-Wall -std=c++17 -O3 
endif


CCxx=mpicxx
CxxINCS=-I$(MFEM_DIR) -I$(MFEM_DIR)/linalg -I$(HYPRE_DIR)/src/hypre/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(EIGEN_DIR)
CxxLDIR=-L$(MFEM_DIR) -L$(HYPRE_DIR)/src/hypre/lib -L$(METIS_DIR) -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/gcc/x86_64-linux-gnu/9
CxxLIBS=-lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lflapack -lfblas -lptesmumps -lptscotchparmetis -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lm -lX11 -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lmfem -lHYPRE -lmetis -lrt -lslepcpep -lslepcsys -lslepceps -lslepclme -lslepcnep -lslepcmfn -lslepcsvd /usr/lib/x86_64-linux-gnu/liblapacke64.a /usr/lib/x86_64-linux-gnu/liblapack64.a -lpetscmat -lpetscsnes -lpetscmat -lpetscsnes -lpetscdm -lpetscvec -lpetscts -lpetsctao -lpetscsys -lpetscksp -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lflapack -lfblas -lslepcpep -lslepcsys -lslepceps -lslepclme -lslepcnep -lslepcmfn -lslepcsvd

CC=mpicc
CINCS=-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(EIGEN_DIR)
CLDIR=
CLIBS=-lm

OpenParEM2D: OpenParEM2D.o  modes.o results.o convergence.o frequencyPlan.o findPoints.o fem2D.o OpenParEMmaterials.o mesh.o fieldPoints.o license.o project.o eigensolve.o Hsolve.o Zsolve.o triplet.o 
	$(CCxx) $(CxxFLAGS) OpenParEM2D.o modes.o results.o convergence.o frequencyPlan.o findPoints.o fem2D.o OpenParEMmaterials.o mesh.o fieldPoints.o license.o project.o eigensolve.o Hsolve.o Zsolve.o triplet.o -o OpenParEM2D  $(CxxLDIR) $(CxxLIBS)

OpenParEM2D.o: OpenParEM2D.cpp
	$(CCxx) $(CxxFLAGS) -c OpenParEM2D.cpp $(CxxINCS)

modes.o: modes.cpp modes.hpp
	$(CCxx) $(CxxFLAGS) -c modes.cpp $(CxxINCS)

results.o: results.cpp results.hpp
	$(CCxx) $(CxxFLAGS) -c results.cpp $(CxxINCS)

convergence.o: convergence.cpp convergence.hpp
	$(CCxx) $(CxxFLAGS) -c convergence.cpp $(CxxINCS)

frequencyPlan.o: frequencyPlan.cpp frequencyPlan.hpp
	$(CCxx) $(CxxFLAGS) -c frequencyPlan.cpp $(CxxINCS)

findPoints.o: findPoints.cpp fem2D.hpp
	$(CCxx) $(CxxFLAGS) -c findPoints.cpp $(CxxINCS)

fem2D.o: fem2D.cpp fem2D.hpp
	$(CCxx) $(CxxFLAGS) -c fem2D.cpp $(CxxINCS)

OpenParEMmaterials.o: OpenParEMmaterials.cpp OpenParEMmaterials.hpp
	$(CCxx) $(CxxFLAGS) -c OpenParEMmaterials.cpp $(CxxINCS)

mesh.o: mesh.cpp mesh.hpp
	$(CCxx) $(CxxFLAGS) -c mesh.cpp $(CxxINCS)

fieldPoints.o: fieldPoints.cpp fieldPoints.hpp
	$(CCxx) $(CxxFLAGS) -c fieldPoints.cpp $(CxxINCS)

license.o: license.cpp
	$(CCxx) $(CxxFLAGS) -c license.cpp $(CxxINCS)

project.o: project.c project.h
	$(CC) $(CFLAGS) -c project.c $(CINCS)

eigensolve.o: eigensolve.c eigensolve.h
	$(CC) $(CFLAGS) -c eigensolve.c $(CINCS) 

Hsolve.o: Hsolve.c Hsolve.h
	$(CC) $(CFLAGS) -c Hsolve.c $(CINCS)

Zsolve.o: Zsolve.c Zsolve.h
	$(CC) $(CFLAGS) -c Zsolve.c $(CINCS)

triplet.o: triplet.c triplet.h
	$(CC) $(CFLAGS) -c triplet.c $(CINCS)

waveguide: waveguide.o
	$(CC) $(CFLAGS) waveguide.o -o waveguide $(CLDIR) $(CLIBS)

waveguide.o: waveguide.c waveguide.h
	$(CC) $(CFLAGS) -c waveguide.c $(CINCS)

process: process.o project.o mesh.o OpenParEMmaterials.o
	$(CCxx) process.o project.o mesh.o OpenParEMmaterials.o -o process $(CxxLDIR) $(CxxLIBS)

process.o: process.cpp process.hpp
	$(CCxx) $(CxxFLAGS) -c process.cpp $(CxxINCS)

builder: builder.o OpenParEMmaterials.o mesh.o
	$(CCxx) $(CxxFLAGS) builder.o mesh.o OpenParEMmaterials.o -o builder $(CxxLDIR) $(CxxLIBS)

builder.o: builder.cpp builder.hpp
	$(CCxx) $(CxxFLAGS) -c builder.cpp $(CxxINCS)

simplify: simplify.o project.o 
	$(CCxx) $(CxxFLAGS) simplify.o project.o -o simplify $(CxxLDIR) $(CxxLIBS)

simplify.o: simplify.cpp simplify.hpp
	$(CCxx) $(CxxFLAGS) -c simplify.cpp $(CxxINCS)

license_process: license_process.o
	$(CCxx) $(CxxFLAGS) license_process.o -o license_process $(CxxLDIR) $(CxxLIBS)

license_process.o: license_process.cpp
	$(CCxx) $(CxxFLAGS) -c license_process.cpp $(CxxINCS)

.PHONY: all clean

all: OpenParEM2D waveguide process builder simplify license_process

install:
	rm -f ../bin/OpenParEM2D
	cp OpenParEM2D ../bin/OpenParEM2D
	rm -f $(HOME)/bin/OpenParEM2D
	mv OpenParEM2D $(HOME)/bin/OpenParEM2D

	rm -f ../bin/process
	cp process ../bin/process
	rm -f $(HOME)/bin/process
	mv process $(HOME)/bin/process

	rm -f ../bin/simplify
	cp simplify ../bin/simplify
	rm -f $(HOME)/bin/simplify
	mv simplify $(HOME)/bin/simplify

	rm -f ../bin/builder
	cp builder ../bin/builder
	rm -f $(HOME)/bin/builder
	mv builder $(HOME)/bin/builder

	rm -f $(HOME)/bin/process.sh
	cp ../scripts/process.sh $(HOME)/bin/process.sh

	rm -f $(HOME)/bin/regression.sh
	cp ../scripts/regression.sh $(HOME)/bin/regression.sh

	rm -f $(HOME)/bin/checkBuilder.sh
	cp ../scripts/checkBuilder.sh $(HOME)/bin/checkBuilder.sh

clean:
	rm -f OpenParEM2D
	rm -f process
	rm -f builder
	rm -f simplify
	rm -f license_process
	rm -f waveguide
	rm -f OpenParEM2D.o
	rm -f OpenParEMmaterials.o
	rm -f eigensolve.o
	rm -f Hsolve.o
	rm -f Zsolve.o
	rm -f triplet.o
	rm -f waveguide.o
	rm -f mesh.o
	rm -f project.o
	rm -f convergence.o
	rm -f fem2D.o
	rm -f frequencyPlan.o
	rm -f modes.o
	rm -f results.o
	rm -f fieldPoints.o
	rm -f findPoints.o
	rm -f process.o
	rm -f simplify.o
	rm -f builder.o
	rm -f license_process.o
	rm -f license.o

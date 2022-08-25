CCxx=mpicxx -Wall -std=c++17 -g -Wl,-rpath,$(PETSC_DIR)/$(PETSC_ARCH)/lib -Wl,-rpath,$(SLEPC_DIR)/$(PETSC_ARCH)/lib
CxxFLAGS=-I$(MFEM_DIR) -I$(MFEM_DIR)/linalg -I$(HYPRE_DIR)/src/hypre/include -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(EIGEN_DIR)
CxxLDIR=-L$(MFEM_DIR) -L$(HYPRE_DIR)/src/hypre/lib -L$(METIS_DIR) -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L/usr/lib/x86_64-linux-gnu
CxxLIBS=-lmfem -lHYPRE -lmetis -lrt -lpetsc -lslepc -lflapack -llapacke64 -lfblas -lm 

CC=mpicc -Wall -g 
CFLAGS=-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -I$(EIGEN_DIR)
CLDIR=-L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L/usr/lib/x86_64-linux-gnu
CLIBS=-lrt -lpetsc -lslepc -lflapack -llapacke64 -lfblas -lm


OpenParEM2D: OpenParEM2D.o  modes.o results.o convergence.o frequencyPlan.o findPoints.o fem2D.o OpenParEMmaterials.o mesh.o fieldPoints.o license.o project.o eigensolve.o Hsolve.o Zsolve.o triplet.o 
	$(CCxx) OpenParEM2D.o modes.o results.o convergence.o frequencyPlan.o findPoints.o fem2D.o OpenParEMmaterials.o mesh.o fieldPoints.o license.o project.o eigensolve.o Hsolve.o Zsolve.o triplet.o -o OpenParEM2D  $(CxxLDIR) $(CxxLIBS)

OpenParEM2D.o: OpenParEM2D.cpp
	$(CCxx) -c OpenParEM2D.cpp $(CxxFLAGS)

modes.o: modes.cpp modes.hpp
	$(CCxx) -c modes.cpp $(CxxFLAGS)

results.o: results.cpp results.hpp
	$(CCxx) -c results.cpp $(CxxFLAGS)

convergence.o: convergence.cpp convergence.hpp
	$(CCxx) -c convergence.cpp $(CxxFLAGS)

frequencyPlan.o: frequencyPlan.cpp frequencyPlan.hpp
	$(CCxx) -c frequencyPlan.cpp $(CxxFLAGS)

findPoints.o: findPoints.cpp fem2D.hpp
	$(CCxx) -c findPoints.cpp $(CxxFLAGS)

fem2D.o: fem2D.cpp fem2D.hpp
	$(CCxx) -c fem2D.cpp $(CxxFLAGS)

OpenParEMmaterials.o: OpenParEMmaterials.cpp OpenParEMmaterials.hpp
	$(CCxx) -c OpenParEMmaterials.cpp $(CxxFLAGS)

mesh.o: mesh.cpp mesh.hpp
	$(CCxx) -c mesh.cpp $(CxxFLAGS)

fieldPoints.o: fieldPoints.cpp fieldPoints.hpp
	$(CCxx) -c fieldPoints.cpp $(CxxFLAGS)

license.o: license.cpp
	$(CCxx) -c license.cpp $(CxxFLAGS)

project.o: project.c project.h
	$(CC) -c project.c $(CFLAGS)

eigensolve.o: eigensolve.c eigensolve.h
	$(CC) -c eigensolve.c $(CFLAGS) 

Hsolve.o: Hsolve.c Hsolve.h
	$(CC) -c Hsolve.c $(CFLAGS)

Zsolve.o: Zsolve.c Zsolve.h
	$(CC) -c Zsolve.c $(CFLAGS)

triplet.o: triplet.c triplet.h
	$(CC) -c triplet.c $(CFLAGS)

waveguide: waveguide.o
	$(CC) waveguide.o -o waveguide $(CLDIR) $(CLIBS)

waveguide.o: waveguide.c waveguide.h
	$(CC) -c waveguide.c $(CFLAGS)

process: process.o project.o mesh.o OpenParEMmaterials.o
	$(CCxx) process.o project.o mesh.o OpenParEMmaterials.o -o process $(CxxLDIR) $(CxxLIBS)

process.o: process.cpp process.hpp
	$(CCxx) -c process.cpp $(CxxFLAGS)

builder: builder.o OpenParEMmaterials.o mesh.o
	$(CCxx) builder.o mesh.o OpenParEMmaterials.o -o builder $(CxxLDIR) $(CxxLIBS)

builder.o: builder.cpp builder.hpp
	$(CCxx) -c builder.cpp $(CxxFLAGS)

simplify: simplify.o project.o 
	$(CCxx) simplify.o project.o -o simplify $(CxxLDIR) $(CxxLIBS)

simplify.o: simplify.cpp simplify.hpp
	$(CCxx) -c simplify.cpp $(CxxFLAGS)

license_process: license_process.o
	$(CCxx) license_process.o -o license_process $(CxxLDIR) $(CxxLIBS)

license_process.o: license_process.cpp
	$(CCxx) -c license_process.cpp $(CxxFLAGS)

.PHONY: all clean

all: OpenParEM2D waveguide process builder simplify license_process

install:
	cp OpenParEM2D $(HOME)/bin
	cp process $(HOME)/bin
	cp simplify $(HOME)/bin
	cp license_process $(HOME)/bin
	cp builder $(HOME)/bin

clean:
	rm -f OpenParEM2D
	rm -f OpenParEM2D.o
	rm -f OpenParEMmaterials.o
	rm -f eigensolve.o
	rm -f Hsolve.o
	rm -f Zsolve.o
	rm -f triplet.o
	rm -f waveguide
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
	rm -f process
	rm -f process.o
	rm -f simplify
	rm -f simplify.o
	rm -f builder
	rm -f builder.o
	rm -f license_process
	rm -f license_process.o
	rm -f license.o

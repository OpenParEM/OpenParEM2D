AC_PREREQ([2.69])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE([-Wall foreign])
AC_CONFIG_SRCDIR([fem2D.cpp])
AC_CONFIG_HEADERS([config.h])
AC_PROG_RANLIB()
AM_PROG_AR()

# Checks for programs.
AC_PROG_CXX(mpicxx)
AC_PROG_CC(mpicc)
AC_PROG_INSTALL

# environment variables for required packages
AC_ARG_VAR(MFEM_DIR,[libmfem.a in $MEM_DIR; includes in $MFEM_DIR and $MFEM_DIR/linalg])
AC_ARG_VAR(HYPRE_DIR,[libHYPRE.a in $HYPRE_DIR/src/hypre/lib; includes in $HYPRE_DIR/src/hypre/include])
AC_ARG_VAR(PETSC_DIR,[libpetscsys.a + others in $PETSC_DIR/$PETSC_ARCH/lib; includes in $PETSC_DIR/include and $PETSC_DIR/$PETSC_ARCH/include])
AC_ARG_VAR(PETSC_ARCH)
AC_ARG_VAR(SLEPC_DIR,[libslepcsys.a + others in $SLEPC_DIR/$PETSC_ARCH/lib; includes in $SLEPC_DIR/include $SLEPC_DIR/$PETSC_ARCH/include])
AC_ARG_VAR(EIGEN_DIR,[includes for blas.h and others in $EIGEN_DIR])
AC_ARG_VAR(METIS_DIR,[libmetis.a in $METIS_DIR])

# Checks for libraries.
AC_CHECK_LIB([c], [malloc])
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([quadmath], [acosq])
AC_CHECK_LIB([X11], [XCreateSimpleWindow])
AC_CHECK_LIB([rt], [clock_getres])
AC_CHECK_LIB([pthread], [pthread_mutex_init])
AC_CHECK_LIB([dl], [dlvsym])

CXXFLAGS="$CXXFLAGS -I/usr/lib/x86_64-linux-gnu/openmpi/include"
CFLAGS="$CFLAGS -I/usr/lib/x86_64-linux-gnu/openmpi/include"
AC_CHECK_LIB([mpi], [MPI_Init])
AC_CHECK_LIB([mpi_mpifh], [mpi_init_])
AC_CHECK_LIB([mpi++], [mpi_init_])
AC_CHECK_LIB([mpi_usempif08], [mpi_init_])
AC_CHECK_LIB([mpi_usempi_ignore_tkr], [mpi_init_])

LDFLAGS="$LDFLAGS -L/usr/lib/x86_64-linux-gnu"
CXXFLAGS="$CXXFLAGS -I$EIGEN_DIR"
CFLAGS="$CFLAGS -I$EIGEN_DIR"
AC_CHECK_LIB([lapack64], [ztgsyl_])
AC_CHECK_LIB([lapacke64], [LAPACKE_zlagsy])

LDFLAGS="$LDFLAGS -L$HYPRE_DIR/src/hypre/lib"
CXXFLAGS="$CXXFLAGS -I$HYPRE_DIR/src/hypre/include"
CFLAGS="$CFLAGS -I$HYPRE_DIR/src/hypre/include"
AC_CHECK_LIB([HYPRE], [hypre_ParCSRMatrixCreate])

LDFLAGS="$LDFLAGS -L$METIS_DIR"
AC_CHECK_LIB([metis], [METIS_PartGraphRecursive])

LDFLAGS="$LDFLAGS -L/usr/lib/gcc/x86_64-linux-gnu/9/"
AC_CHECK_LIB([stdc++], [_ZNSt9bad_allocD2Ev])
AC_CHECK_LIB([gfortran], [_gfortran_rand])
AC_CHECK_LIB([gcc_s], [__clzdi2])

LDFLAGS="$LDFLAGS -L$PETSC_DIR/$PETSC_ARCH/lib"
CXXFLAGS="$CXXFLAGS -I$PETSC_DIR/include -I$PETSC_DIR/$PETSC_ARCH/include"
CFLAGS="$CFLAGS -I$PETSC_DIR/include -I$PETSC_DIR/$PETSC_ARCH/include"
AC_CHECK_LIB([fblas], [cswap_])
AC_CHECK_LIB([flapack], [cungqr_])
AC_CHECK_LIB([scalapack], [pbchkvect])
AC_CHECK_LIB([scotcherr], [SCOTCH_errorPrint])
AC_CHECK_LIB([ptscotcherr], [SCOTCH_errorPrint])
AC_CHECK_LIB([scotch], [scotchfstratinit_])
AC_CHECK_LIB([ptscotch], [scotchfdgraphorderinit_])
AC_CHECK_LIB([ptscotchparmetisv3], [scotch_parmetis_v3_partkway])
AC_CHECK_LIB([pord], [printFactorMtx])
AC_CHECK_LIB([esmumps], [esmumps])
AC_CHECK_LIB([ptesmumps], [esmumpsv])
AC_CHECK_LIB([mumps_common], [mumps_abort_])
AC_CHECK_LIB([cmumps], [cmumps_root_solve_])
AC_CHECK_LIB([dmumps], [dmumps_gather_root_])
AC_CHECK_LIB([smumps], [smumps_gather_root_])
AC_CHECK_LIB([zmumps], [zmumps_gather_root_])
AC_CHECK_LIB([petscsys], [PetscObjectComm])
AC_CHECK_LIB([petscvec], [vecsum_])
AC_CHECK_LIB([petscmat], [matcreateaij_])
AC_CHECK_LIB([petsctao], [MatSMFResetRowColumn])
AC_CHECK_LIB([petscsnes], [SNESLineSearchCreate_Shell])
AC_CHECK_LIB([petscdm], [DMFieldCreateDS])
AC_CHECK_LIB([petscksp], [PCFieldSplitSetOffDiagUseAmat])
AC_CHECK_LIB([petscts], [dmplexlandauijacobian_])

LDFLAGS="$LDFLAGS -L$SLEPC_DIR/$PETSC_ARCH/lib"
CXXFLAGS="$CXXFLAGS -I$SLEPC_DIR/include -I$SLEPC_DIR/$PETSC_ARCH/include"
CFLAGS="$CFLAGS -I$SLEPC_DIR/include -I$SLEPC_DIR/$PETSC_ARCH/include"
AC_CHECK_LIB([slepcsys], [STPostSolve_Sinvert])
AC_CHECK_LIB([slepclme], [LMECreate_Krylov])
AC_CHECK_LIB([slepceps], [EPSCreate_KrylovSchur])
AC_CHECK_LIB([slepcmfn], [MFNCreate_Krylov])
AC_CHECK_LIB([slepcpep], [PEPSTOARSetDimensions])
AC_CHECK_LIB([slepcnep], [NEPSolve_NArnoldi])
AC_CHECK_LIB([slepcsvd], [SVDSetUp_LAPACK])

LDFLAGS="$LDFLAGS -L$MFEM_DIR"
CXXFLAGS="$CXXFLAGS -I$MFEM_DIR -I$MFEM_DIR/linalg"
AC_CHECK_LIB([mfem], [_ZN4mfem11DenseMatrix3NegEv])

LDFLAGS="$LDFLAGS -L$HOME/lib"
CXXFLAGS="$CXXFLAGS -I../../OpenParEMCommon/src"
CFLAGS="$CXXFLAGS -I../../OpenParEMCommon/src"
AC_CHECK_LIB([OpenParEMCommon], [prefix])

# Checks for header files.
AC_CHECK_HEADERS([stddef.h stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt])

# force use of c++17
m4_include([ax_cxx_compile_stdcxx.m4])
AX_CXX_COMPILE_STDCXX([17],[noext],[mandatory])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

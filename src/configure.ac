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
AC_CHECK_LIB([c], [malloc], ,AC_MSG_ERROR([cannot link to libc.a; check config.log]))
AC_CHECK_LIB([m], [cos], ,AC_MSG_ERROR([cannot link to libm.a; check config.log]))
AC_CHECK_LIB([quadmath], [acosq], ,AC_MSG_ERROR([cannot link to libquadmath.a; check config.log]))
AC_CHECK_LIB([X11], [XCreateSimpleWindow], ,AC_MSG_ERROR([cannot link to libX11.a; check config.log]))
AC_CHECK_LIB([rt], [clock_getres], ,AC_MSG_ERROR([cannot link to librt.a; check config.log]))
AC_CHECK_LIB([pthread], [pthread_mutex_init], ,AC_MSG_ERROR([cannot link to libpthread.a; check config.log]))
AC_CHECK_LIB([dl], [dlvsym], ,AC_MSG_ERROR([cannot link to libdl.a; check config.log]))

CXXFLAGS="$CXXFLAGS -I/usr/lib/x86_64-linux-gnu/openmpi/include"
CFLAGS="$CFLAGS -I/usr/lib/x86_64-linux-gnu/openmpi/include"
AC_CHECK_LIB([mpi], [MPI_Init], ,AC_MSG_ERROR([cannot link to libmpi.a; check config.log]))
AC_CHECK_LIB([mpi_mpifh], [mpi_init_], ,AC_MSG_ERROR([cannot link to libmpi_mpifh.a; check config.log]))
AC_CHECK_LIB([mpi++], [mpi_init_], ,AC_MSG_ERROR([cannot link to libmpi++.a; check config.log]))
AC_CHECK_LIB([mpi_usempif08], [mpi_init_], ,AC_MSG_ERROR([cannot link to libmpi_usempif08.a; check config.log]))
AC_CHECK_LIB([mpi_usempi_ignore_tkr], [mpi_init_], ,AC_MSG_ERROR([cannot link to libmpi_usempi_ignore_tkr.a; check config.log]))

LDFLAGS="$LDFLAGS -L/usr/lib/x86_64-linux-gnu"
CXXFLAGS="$CXXFLAGS -I$EIGEN_DIR"
CFLAGS="$CFLAGS -I$EIGEN_DIR"
AC_CHECK_LIB([lapack64], [ztgsyl_], ,AC_MSG_ERROR([cannot link to liblapack64.a; check config.log]))
AC_CHECK_LIB([lapacke64], [LAPACKE_zlagsy], ,AC_MSG_ERROR([cannot link to liblapacke64.a; check config.log]))

LDFLAGS="$LDFLAGS -L$HYPRE_DIR/src/hypre/lib"
CXXFLAGS="$CXXFLAGS -I$HYPRE_DIR/src/hypre/include"
CFLAGS="$CFLAGS -I$HYPRE_DIR/src/hypre/include"
AC_CHECK_LIB([HYPRE], [hypre_ParCSRMatrixCreate], ,AC_MSG_ERROR([cannot link to libHYPRE.a; check config.log]))

LDFLAGS="$LDFLAGS -L$METIS_DIR"
AC_CHECK_LIB([metis], [METIS_PartGraphRecursive], ,AC_MSG_ERROR([cannot link to libmetis.a; check config.log]))

LDFLAGS="$LDFLAGS -L/usr/lib/gcc/x86_64-linux-gnu/9/"
AC_CHECK_LIB([gfortran], [_gfortran_rand], ,AC_MSG_ERROR([cannot link to libgfortran.a; check config.log]))
AC_CHECK_LIB([gcc_s], [__clzdi2], ,AC_MSG_ERROR([cannot link to libgcc_s.a; check config.log]))

AC_LANG([C++])
LIBS="$LIBS -lstdc++"
AC_RUN_IFELSE([AC_LANG_PROGRAM([[
    using namespace std; ]], [[ double *a=new double ]])],
       [       AC_MSG_RESULT([checking for new in -lstdc++... yes])
               AC_DEFINE_UNQUOTED([HAVE_LIBSTDC__], [1], [Set to 1 when libstdc++ is found])
       ],
       [       AC_MSG_RESULT([checking for new in -lstdc++... no])
               AC_MSG_ERROR([cannot link to libstdc++.a; check config.log])
       ],
       []
)
AC_LANG([C])

LDFLAGS="$LDFLAGS -L$PETSC_DIR/$PETSC_ARCH/lib"
CXXFLAGS="$CXXFLAGS -I$PETSC_DIR/include -I$PETSC_DIR/$PETSC_ARCH/include"
CFLAGS="$CFLAGS -I$PETSC_DIR/include -I$PETSC_DIR/$PETSC_ARCH/include"
AC_CHECK_LIB([fblas], [cswap_], ,AC_MSG_ERROR([cannot link to libfblas.a; check config.log]))
AC_CHECK_LIB([flapack], [cungqr_], ,AC_MSG_ERROR([cannot link to libflapack.a; check config.log]))
AC_CHECK_LIB([scalapack], [pbchkvect], ,AC_MSG_ERROR([cannot link to libscalapack.a; check config.log]))
AC_CHECK_LIB([scotcherr], [SCOTCH_errorPrint], ,AC_MSG_ERROR([cannot link to libscotcherr.a; check config.log]))
AC_CHECK_LIB([ptscotcherr], [SCOTCH_errorPrint], ,AC_MSG_ERROR([cannot link to libptscotcherr.a; check config.log]))
AC_CHECK_LIB([scotch], [scotchfstratinit_], ,AC_MSG_ERROR([cannot link to libscotch.a; check config.log]))
AC_CHECK_LIB([ptscotch], [scotchfdgraphorderinit_], ,AC_MSG_ERROR([cannot link to libptscotch.a; check config.log]))
AC_CHECK_LIB([ptscotchparmetisv3], [scotch_parmetis_v3_partkway], ,AC_MSG_ERROR([cannot link to libptscotchparmetisv3.a; check config.log]))
AC_CHECK_LIB([pord], [printFactorMtx], ,AC_MSG_ERROR([cannot link to libpord.a; check config.log]))
AC_CHECK_LIB([esmumps], [esmumps], ,AC_MSG_ERROR([cannot link to libesmumps.a; check config.log]))
AC_CHECK_LIB([ptesmumps], [esmumpsv], ,AC_MSG_ERROR([cannot link to libptesmumps.a; check config.log]))
AC_CHECK_LIB([mumps_common], [mumps_abort_], ,AC_MSG_ERROR([cannot link to libmumps_common.a; check config.log]))
AC_CHECK_LIB([cmumps], [cmumps_root_solve_], ,AC_MSG_ERROR([cannot link to libcmumps.a; check config.log]))
AC_CHECK_LIB([dmumps], [dmumps_gather_root_], ,AC_MSG_ERROR([cannot link to libdmumps.a; check config.log]))
AC_CHECK_LIB([smumps], [smumps_gather_root_], ,AC_MSG_ERROR([cannot link to libsmumps.a; check config.log]))
AC_CHECK_LIB([zmumps], [zmumps_gather_root_], ,AC_MSG_ERROR([cannot link to libzmumps.a; check config.log]))
AC_CHECK_LIB([petscsys], [PetscObjectComm], ,AC_MSG_ERROR([cannot link to libpetscsys.a; check config.log]))
AC_CHECK_LIB([petscvec], [vecsum_], ,AC_MSG_ERROR([cannot link to libpetscvec.a; check config.log]))
AC_CHECK_LIB([petscmat], [matcreateaij_], ,AC_MSG_ERROR([cannot link to libpetscmat.a; check config.log]))
AC_CHECK_LIB([petsctao], [MatSMFResetRowColumn], ,AC_MSG_ERROR([cannot link to libpetsctao.a; check config.log]))
AC_CHECK_LIB([petscsnes], [SNESLineSearchCreate_Shell], ,AC_MSG_ERROR([cannot link to libpetscsnes.a; check config.log]))
AC_CHECK_LIB([petscdm], [DMFieldCreateDS], ,AC_MSG_ERROR([cannot link to libpetscdm.a; check config.log]))
AC_CHECK_LIB([petscksp], [PCFieldSplitSetOffDiagUseAmat], ,AC_MSG_ERROR([cannot link to libpetscksp.a; check config.log]))
AC_CHECK_LIB([petscts], [dmplexlandauijacobian_], ,AC_MSG_ERROR([cannot link to libpetscts.a; check config.log]))

LDFLAGS="$LDFLAGS -L$SLEPC_DIR/$PETSC_ARCH/lib"
CXXFLAGS="$CXXFLAGS -I$SLEPC_DIR/include -I$SLEPC_DIR/$PETSC_ARCH/include"
CFLAGS="$CFLAGS -I$SLEPC_DIR/include -I$SLEPC_DIR/$PETSC_ARCH/include"
AC_CHECK_LIB([slepcsys], [STPostSolve_Sinvert], ,AC_MSG_ERROR([cannot link to libslepcsys.a; check config.log]))
AC_CHECK_LIB([slepclme], [LMECreate_Krylov], ,AC_MSG_ERROR([cannot link to libslepclme.a; check config.log]))
AC_CHECK_LIB([slepceps], [EPSCreate_KrylovSchur], ,AC_MSG_ERROR([cannot link to libslepceps.a; check config.log]))
AC_CHECK_LIB([slepcmfn], [MFNCreate_Krylov], ,AC_MSG_ERROR([cannot link to libslepcmfn.a; check config.log]))
AC_CHECK_LIB([slepcpep], [PEPSTOARSetDimensions], ,AC_MSG_ERROR([cannot link to libslepcpep.a; check config.log]))
AC_CHECK_LIB([slepcnep], [NEPSolve_NArnoldi], ,AC_MSG_ERROR([cannot link to libslepcnep.a; check config.log]))
AC_CHECK_LIB([slepcsvd], [SVDSetUp_LAPACK], ,AC_MSG_ERROR([cannot link to libslepcsvd.a; check config.log]))

AC_LANG([C++])
LDFLAGS="$LDFLAGS -L$MFEM_DIR"
CXXFLAGS="$CXXFLAGS -I$MFEM_DIR -I$MFEM_DIR/linalg"
LIBS="$LIBS -lmfem -lHYPRE -lmetis"
AC_RUN_IFELSE([AC_LANG_PROGRAM([[
    #include "mfem.hpp"
    using namespace mfem; ]], [[ HypreParMatrix *mat; ]])],
       [       AC_MSG_RESULT([checking for HypreParMatrix in -lmfem... yes])
               AC_DEFINE_UNQUOTED([HAVE_LIBMFEM], [1], [Set to 1 when libmfem.a is found])
       ],
       [       AC_MSG_RESULT([checking for HypreParMatrix in -lmfem... no])
               AC_MSG_ERROR([cannot link to libmfem.a; check config.log])
       ],
       []
)
AC_LANG([C])

LDFLAGS="$LDFLAGS -L$HOME/lib"
CXXFLAGS="$CXXFLAGS -I../../OpenParEMCommon/src"
CFLAGS="$CXXFLAGS -I../../OpenParEMCommon/src"
AC_CHECK_LIB([OpenParEMCommon], [prefix], ,AC_MSG_ERROR([cannot link to libOpenParEMCommon.a; check config.log]))

# Checks for header files.
AC_CHECK_HEADERS([stddef.h stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt getcwd strdup strtoull])

# force use of c++17
m4_include([ax_cxx_compile_stdcxx.m4])
AX_CXX_COMPILE_STDCXX([17],[noext],[mandatory])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

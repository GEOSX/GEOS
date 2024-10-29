/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BlasLapackFunctions.h
 */

#ifndef GEOS_DENSELINEARALGEBRA_INTERFACES_BLASLAPACK_BLASLAPACKFUNCTIONS_HPP_
#define GEOS_DENSELINEARALGEBRA_INTERFACES_BLASLAPACK_BLASLAPACKFUNCTIONS_HPP_

/// This macro provide a flexible interface for Fortran naming convetion for compiled objects
#ifdef FORTRAN_MANGLE_NO_UNDERSCORE
#define FORTRAN_MANGLE( name ) name
#else
#define FORTRAN_MANGLE( name ) name ## _
#endif

#ifdef __cplusplus
extern "C"
{

#define GEOS_dasum FORTRAN_MANGLE( dasum )
double GEOS_dasum( int const * N,
                   double const * DX,
                   int const * INCX );

#define GEOS_daxpy FORTRAN_MANGLE( daxpy )
void GEOS_daxpy( int const * N,
                 double const * DA,
                 double const * DX,
                 int const * INCX,
                 double * DY,
                 int const * INCY );

#define GEOS_dcopy FORTRAN_MANGLE( dcopy )
void GEOS_dcopy( int const * N,
                 double const * DX,
                 int const * INCX,
                 double * DY,
                 int const * INCY );

#define GEOS_ddot FORTRAN_MANGLE( ddot )
double GEOS_ddot( int const * N,
                  double const * DX,
                  int const * INCX,
                  double const * DY,
                  int const * INCY );

#define GEOS_idamax FORTRAN_MANGLE( idamax )
int GEOS_idamax( int const * N,
                 double const * DX,
                 int const * INCX );

#define GEOS_dgemm FORTRAN_MANGLE( dgemm )
void GEOS_dgemm( char const * TRANSA,
                 char const * TRANSB,
                 int const * M,
                 int const * N,
                 int const * K,
                 double const * ALPHA,
                 double const * A,
                 int const * LDA,
                 double const * B,
                 int const * LDB,
                 double const * BETA,
                 double * C,
                 int const * LDC );

#define GEOS_dgetrf FORTRAN_MANGLE( dgetrf )
void GEOS_dgetrf( int const * M,
                  int const * N,
                  double * A,
                  int const * LDA,
                  int * IPIV,
                  int * INFO );

#define GEOS_dgetri FORTRAN_MANGLE( dgetri )
void GEOS_dgetri( int const * N,
                  double * A,
                  int const * LDA,
                  int const * IPIV,
                  double * WORK,
                  int const * LWORK,
                  int * INFO );

#define GEOS_dlange FORTRAN_MANGLE( dlange )
double GEOS_dlange( char const * NORM,
                    int const * M,
                    int const * N,
                    double const * A,
                    int const * LDA,
                    double * WORK );

#define GEOS_dlarnv FORTRAN_MANGLE( dlarnv )
void GEOS_dlarnv( int const * IDIST,
                  int * ISEED,
                  int const * N,
                  double * X );

#define GEOS_dnrm2 FORTRAN_MANGLE( dnrm2 )
double GEOS_dnrm2( int const * N,
                   double const * X,
                   int const * INCX );

#define GEOS_dscal FORTRAN_MANGLE( dscal )
void GEOS_dscal( int const * N,
                 double const * DA,
                 double * DX,
                 int const * INCX );

#define GEOS_dgesvd FORTRAN_MANGLE( dgesvd )
void GEOS_dgesvd( char const * JOBU,
                  char const * JOBVT,
                  int const * M,
                  int const * N,
                  double * A,
                  int const * LDA,
                  double * S,
                  double * U,
                  int const * LDU,
                  double * VT,
                  int const * LDVT,
                  double * WKOPT,
                  int const * LWORK,
                  int * INFO );

#define GEOS_dgeev FORTRAN_MANGLE( dgeev )
void GEOS_dgeev( char const * JOBVL,
                 char const * JOBVR,
                 int const * N,
                 double * A,
                 int const * LDA,
                 double * WR,
                 double * WI,
                 double * VL,
                 int const * LDVL,
                 double * VR,
                 int const * LDVR,
                 double * WORK,
                 int * LWORK,
                 int * INFO );

#define GEOS_dgetrs FORTRAN_MANGLE( dgetrs )
void GEOS_dgetrs( char const * TRANS,
                  int const * N,
                  int const * NRHS,
                  double * A,
                  int const * LDA,
                  int const * IPIV,
                  double * B,
                  int const * LDB,
                  int * INFO );

#define GEOS_dgels FORTRAN_MANGLE( dgels )
void GEOS_dgels( char const * TRANS,
                 int const * M,
                 int const * N,
                 int const * NRHS,
                 double * A,
                 int const * LDA,
                 double * B,
                 int const * LDB,
                 double * WORK,
                 int const * LWORK,
                 int * INFO );

}
#endif

#endif //GEOS_DENSELINEARALGEBRA_INTERFACES_BLASLAPACKFUNCTIONS_HPP_

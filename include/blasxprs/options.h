// set some defaults

#ifdef BLASXPRS_BLAS
	// set BLASXPRS_WITH_XYZ_BLAS according to BLASXPRS_BLAS in {openblas, mkl, acml, atlas, netlib}
#endif

// handle AMD ACML constraints
#ifdef BLASXPRS_WITH_AMD_ACML
  #ifdef BLASXPRS_MTRX_ROW_MAJOR
	#undef BLASXPRS_MTRX_ROW_MAJOR
	#ifndef BLASXPRS_MTRX_COLUMN_MAJOR
	  #define BLASXPRS_MTRX_COLUMN_MAJOR
	#endif
  #endif
#endif

#if defined( BLASXPRS_WITH_AMD_ACML ) || defined( BLASXPRS_WITH_ATLAS )
  #ifdef BLASXPRS_USE_F77BLAS
    #undef BLASXPRS_USE_F77BLAS
	#ifndef BLASXPRS_USE_CBLAS
      #define BLASXPRS_USE_CBLAS
	#endif
  #endif
#endif

// Use CBLAS
#if defined( BLASXPRS_USE_CBLAS )
  // Define CBLAS flags (layout, transpose, etc.)
  // OpenBLAS & Intel MKL
  #ifdef BLASXPRS_WITH_AMD_ACML
    // Layout is fixed to be column-major for AMD ACML
    // CBLAS_TRANSPOSE
	#define BLASXPRS_CBLAS_TRANSPOSE char
    #define BLASXPRS_CBLAS_NOTRANS 'N'
    #define BLASXPRS_CBLAS_TRANS 'T'
  #else 
    // enum CBLAS_LAYOUT
	#ifdef BLASXPRS_WITH_ATLAS
      #define BLASXPRS_CBLAS_LAYOUT CBLAS_ORDER
	#else
	  #define BLASXPRS_CBLAS_LAYOUT CBLAS_LAYOUT
	#endif
	#define BLASXPRS_CBLAS_ROW_MAJOR CblasRowMajor
	#define BLASXPRS_CBLAS_COLUMN_MAJOR CblasColMajor
    // enum CBLAS_TRANSPOSE
	#define BLASXPRS_CBLAS_TRANSPOSE CBLAS_TRANSPOSE
    #define BLASXPRS_CBLAS_NOTRANS CblasNoTrans
    #define BLASXPRS_CBLAS_TRANS CblasTrans
  // AMD ACML
  #endif
  // Include CBLAS headers
  #if defined( BLASXPRS_WITH_OPENBLAS ) || defined( BLASXPRS_WITH_NETLIB_BLAS ) 
    #include <cblas.h>
  #elif defined( BLASXPRS_WITH_ATLAS )
    extern "C" {
		#include <cblas.h>
	}
  #elif defined( BLASXPRS_WITH_INTEL_MKL )
	#include <mkl_cblas.h>
  #elif defined( BLASXPRS_WITH_AMD_ACML )
    #include <acml.h>
  #endif
// Use Fortran77 BLAS directly
#elif defined( BLASXPRS_USE_F77BLAS )
  // Restrict to column major representation
  #ifdef BLASXPRS_MTRX_ROW_MAJOR
	#undef BLASXPRS_MTRX_ROW_MAJOR
  #endif
  #ifndef BLASXPRS_MTRX_COLUMN_MAJOR
	#define BLASXPRS_MTRX_COLUMN_MAJOR
  #endif
  // Include Fortran77 C headers
  #if defined( BLASXPRS_WITH_OPENBLAS )
  	#include <f77blas.h>
  #elif defined( BLASXPRS_WITH_NETLIB_BLAS )
    #include <cblas_mangling.h>
    #include <cblas_f77.h>
  #elif defined( BLASXPRS_WITH_INTEL_MKL )
	#include <mkl_blas.h>
  #endif
  #ifdef BLASXPRS_WITH_INTEL_MKL
    #define F77BLAS_CALL(routine)     routine
  #else
    #define F77BLAS_CALL(routine)     routine##_
  #endif
#endif

#ifdef BLASXPRS_WITH_ARMADILLO_UTIL
  #include <armadillo>
#endif

#ifdef BLASXPRS_WITH_BLAZE_UTIL
  #include <blaze/Math.h>
#endif

#ifdef BLASXPRS_WITH_EIGEN_UTIL
  #include <Eigen/Dense>
#endif

// Defines for layout and layout-dependent access 
#if defined( BLASXPRS_MTRX_ROW_MAJOR )
  #define BLASXPRS_MTRX_LAYOUT                      BLASXPRS_CBLAS_ROW_MAJOR
  #define BLASXPRS_MTRX_ROW(i,M,ncols)              M + i*ncols 
  #define BLASXPRS_MTRX_ENTRY(i,j,M,ncols)          *(BLASXPRS_MTRX_ROW(i,M,ncols) + j)
#elif defined( BLASXPRS_MTRX_COLUMN_MAJOR )
  #define BLASXPRS_MTRX_LAYOUT                      BLASXPRS_CBLAS_COLUMN_MAJOR
  #define BLASXPRS_MTRX_COLUMN(j,M,nrows)           M + j*nrows
  #define BLASXPRS_MTRX_ENTRY(i,j,M,nrows)          *(BLASXPRS_MTRX_COLUMN(j,M,nrows) + i)
#endif


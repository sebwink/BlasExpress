#ifndef BLASXPRS_REAL_XGEMM_H
#define BLASXPRS_REAL_XGEMM_H

#include <cassert>
#include <iostream>

#include <blasxprs/matrix.h>

#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif

using namespace blasxprs;

using std::cout;

template <typename ValueType>
void fill_with_ixj(Matrix<ValueType>& M) {
	for (std::size_t i = 0; i < M.num_rows(); ++i)
		for(std::size_t j = 0; j < M.num_cols(); ++j) 
			M(i,j) = (i+1)*(j+1);
}

template <typename ValueType>
void xgemm_test() {
	Matrix<ValueType> I { Identity<ValueType>(5) };

	//ValueType a = 2.0;
	//ValueType b = -2.0;

    Matrix<ValueType> A(3,3);
	fill_with_ixj(A);
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	arma::Mat<ValueType> AA = A.arma();
	assert( A == AA );
	assert( AA == A );
#endif
	
    Matrix<ValueType> B(3,3);
	fill_with_ixj(B); 
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	arma::Mat<ValueType> BA = B.arma();
	assert( B == BA );
	assert( BA == B );
#endif
    
	Matrix<ValueType> C(3,3);
    fill_with_ixj(C);
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	arma::Mat<ValueType> CA = C.arma();
	assert( C == CA );
	assert( CA == C );
#endif
/*
	Matrix<ValueType> case1 = A*B;
	Matrix<ValueType> case2 = A.t()*B;
	Matrix<ValueType> case3 = A*B.t();
	Matrix<ValueType> case4 = A.t()*B.t();
	Matrix<ValueType> case5 = a*A*B;
	Matrix<ValueType> case6 = a*A.t()*B;
	Matrix<ValueType> case7 = a*A*B.t();
	Matrix<ValueType> case8 = a*A.t()*B.t();
	Matrix<ValueType> case9( A*B );
	Matrix<ValueType> case10( A.t()*B );
	Matrix<ValueType> case11( A*B.t() );
	Matrix<ValueType> case12( A.t()*B.t() );
	Matrix<ValueType> case13( a*A*B );
	Matrix<ValueType> case14( a*A.t()*B );
	Matrix<ValueType> case15( a*A*B.t() );
	Matrix<ValueType> case16( a*A.t()*B.t() );
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	arma::Mat<ValueType> case1A = AA*BA;
	cassert( case1 == case1A )
	arma::Mat<Type> case2 = AA.t()*BA;
	cassert( case2 == case2A )
	arma::Mat<ValueType> case3 = AA*BA.t();
	cassert( case3 == case3A )
	arma::Mat<ValueType> case4 = AA.t()*BA.t();
	cassert( case5 == case5A )
	arma::Mat<ValueType> case5 = a*AA*BA;
	cassert( case6 == case6A )
	arma::Mat<ValueType> case6 = a*AA.t()*BA;
	cassert( case7 == case7A )
	arma::Mat<ValueType> case7 = a*AA*BA.t();
	cassert( case7 == case7A )
	arma::Mat<ValueType> case8 = a*AA.t()*BA.t();
	cassert( case8 == case8A )
	arma::Mat<ValueType> case9( AA*BA );
	cassert( case9 == case9A )
	arma::Mat<ValueType> case10( AA.t()*BA );
	cassert( case10 == case10A )
	arma::Mat<ValueType> case11( AA*BA.t() );
	cassert( case11 == case11A )
	arma::Mat<ValueType> case12( AA.t()*BA.t() );
	cassert( case12 == case12A )
	arma::Mat<ValueType> case13( a*AA*BA );
	cassert( case13 == case13A )
	arma::Mat<ValueType> case14( a*AA.t()*BA );
	cassert( case14 == case14A )
	arma::Mat<ValueType> case15( a*AA*BA.t() );
	cassert( case15 == case15A )
	arma::Mat<ValueType> case16( a*AA.t()*BA.t() );
	cassert( case16 == case16A )
#endif
*/
/*
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
	Matrix<ValueType> case1 = A*B;
	Matrix<ValueType> case2 = A.t()*B;
	Matrix<ValueType> case3 = A*B.t();
	Matrix<ValueType> case4 = A.t()*B.t();
	Matrix<ValueType> case5 = a*A*B;
	Matrix<ValueType> case6 = a*A.t()*B;
	Matrix<ValueType> case7 = a*A*B.t();
	Matrix<ValueType> case8 = a*A.t()*B.t();
	Matrix<ValueType> case9( A*B );
	Matrix<ValueType> case10( A.t()*B );
	Matrix<ValueType> case11( A*B.t() );
	Matrix<ValueType> case12( A.t()*B.t() );
	Matrix<ValueType> case13( a*A*B );
	Matrix<ValueType> case14( a*A.t()*B );
	Matrix<ValueType> case15( a*A*B.t() );
	Matrix<ValueType> case16( a*A.t()*B.t() );
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
	Matrix<ValueType> case1 = A*B;
	Matrix<ValueType> case2 = A.t()*B;
	Matrix<ValueType> case3 = A*B.t();
	Matrix<ValueType> case4 = A.t()*B.t();
	Matrix<ValueType> case5 = a*A*B;
	Matrix<ValueType> case6 = a*A.t()*B;
	Matrix<ValueType> case7 = a*A*B.t();
	Matrix<ValueType> case8 = a*A.t()*B.t();
	Matrix<ValueType> case9( A*B );
	Matrix<ValueType> case10( A.t()*B );
	Matrix<ValueType> case11( A*B.t() );
	Matrix<ValueType> case12( A.t()*B.t() );
	Matrix<ValueType> case13( a*A*B );
	Matrix<ValueType> case14( a*A.t()*B );
	Matrix<ValueType> case15( a*A*B.t() );
	Matrix<ValueType> case16( a*A.t()*B.t() );
#endif
	C = A*B + C;           // case 17
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C += A*B;              // case 18
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A.t()*B + C;       // case 19
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C += A.t()*B;          // case 20
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A*B.t() + C;       // case 21
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C += A*B.t();          // case 22
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A.t()*B.t() + C;   // case 23
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C += A.t()*B.t();      // case 24
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A*B + C;         // case 25
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C += a*A*B;            // case 26
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A.t()*B + C;       // case 27 
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C += a*A.t()*B;          // case 28
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A*B.t() + C;       // case 29
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C += a*A*B.t();          // case 30
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A.t()*B.t() + C;   // case 31
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C += a*A.t()*B.t();      // case 32
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A*B - C;           // case 33
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C -= A*B;              // case 34
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A.t()*B - C;       // case 35
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C -= A.t()*B;          // case 36
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A*B.t() - C;       // case 37
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C -= A*B.t();          // case 38
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A.t()*B.t() - C;   // case 39
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C -= A.t()*B.t();      // case 40
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A*B - C;         // case 41
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C -= a*A*B;            // case 41
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A.t()*B - C;       // case 42 
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C -= a*A.t()*B;          // case 43
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A*B.t() - C;       // case 44
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C -= a*A*B.t();          // case 45
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A.t()*B.t() - C;   // case 46
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C -= a*A.t()*B.t();      // case 47
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A*B + b*C;           // case 48
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A.t()*B + b*C;       // case 49
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A*B.t() + b*C;       // case 50
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A.t()*B.t() + b*C;   // case 51
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A*B + b*C;         // case 52
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A.t()*B + b*C;       // case 53 
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A*B.t() + b*C;       // case 54
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A.t()*B.t() + b*C;   // case 55
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A*B - b*C;           // case 56
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A.t()*B - b*C;       // case 57
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A*B.t() - b*C;       // case 58
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = A.t()*B.t() - b*C;   // case 59
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A*B - b*C;         // case 60
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A.t()*B - b*C;       // case 61 
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A*B.t() - b*C;       // case 62
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = a*A.t()*B.t() - b*C;   // case 63
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	Matrix<ValueType> case64 = -A*B;
	Matrix<ValueType> case65 = -A.t()*B;
	Matrix<ValueType> case66 = -A*B.t();
	Matrix<ValueType> case67 = -A.t()*B.t();
	Matrix<ValueType> case68 = -a*A*B;
	Matrix<ValueType> case69 = -a*A.t()*B;
	Matrix<ValueType> case70 = -a*A*B.t();
	Matrix<ValueType> case71 = -a*A.t()*B.t();
	Matrix<ValueType> case72( -A*B );
	Matrix<ValueType> case73( -A.t()*B );
	Matrix<ValueType> case74( -A*B.t() );
	Matrix<ValueType> case75( -A.t()*B.t() );
	Matrix<ValueType> case76( -a*A*B );
	Matrix<ValueType> case77( -a*A.t()*B );
	Matrix<ValueType> case78( -a*A*B.t() );
	Matrix<ValueType> case79( -a*A.t()*B.t() );
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	Matrix<ValueType> case64 = -A*B;
	Matrix<ValueType> case65 = -A.t()*B;
	Matrix<ValueType> case66 = -A*B.t();
	Matrix<ValueType> case67 = -A.t()*B.t();
	Matrix<ValueType> case68 = -a*A*B;
	Matrix<ValueType> case69 = -a*A.t()*B;
	Matrix<ValueType> case70 = -a*A*B.t();
	Matrix<ValueType> case71 = -a*A.t()*B.t();
	Matrix<ValueType> case72( -A*B );
	Matrix<ValueType> case73( -A.t()*B );
	Matrix<ValueType> case74( -A*B.t() );
	Matrix<ValueType> case75( -A.t()*B.t() );
	Matrix<ValueType> case76( -a*A*B );
	Matrix<ValueType> case77( -a*A.t()*B );
	Matrix<ValueType> case78( -a*A*B.t() );
	Matrix<ValueType> case79( -a*A.t()*B.t() );
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
	Matrix<ValueType> case64 = -A*B;
	Matrix<ValueType> case65 = -A.t()*B;
	Matrix<ValueType> case66 = -A*B.t();
	Matrix<ValueType> case67 = -A.t()*B.t();
	Matrix<ValueType> case68 = -a*A*B;
	Matrix<ValueType> case69 = -a*A.t()*B;
	Matrix<ValueType> case70 = -a*A*B.t();
	Matrix<ValueType> case71 = -a*A.t()*B.t();
	Matrix<ValueType> case72( -A*B );
	Matrix<ValueType> case73( -A.t()*B );
	Matrix<ValueType> case74( -A*B.t() );
	Matrix<ValueType> case75( -A.t()*B.t() );
	Matrix<ValueType> case76( -a*A*B );
	Matrix<ValueType> case77( -a*A.t()*B );
	Matrix<ValueType> case78( -a*A*B.t() );
	Matrix<ValueType> case79( -a*A.t()*B.t() );
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
	Matrix<ValueType> case64 = -A*B;
	Matrix<ValueType> case65 = -A.t()*B;
	Matrix<ValueType> case66 = -A*B.t();
	Matrix<ValueType> case67 = -A.t()*B.t();
	Matrix<ValueType> case68 = -a*A*B;
	Matrix<ValueType> case69 = -a*A.t()*B;
	Matrix<ValueType> case70 = -a*A*B.t();
	Matrix<ValueType> case71 = -a*A.t()*B.t();
	Matrix<ValueType> case72( -A*B );
	Matrix<ValueType> case73( -A.t()*B );
	Matrix<ValueType> case74( -A*B.t() );
	Matrix<ValueType> case75( -A.t()*B.t() );
	Matrix<ValueType> case76( -a*A*B );
	Matrix<ValueType> case77( -a*A.t()*B );
	Matrix<ValueType> case78( -a*A*B.t() );
	Matrix<ValueType> case79( -a*A.t()*B.t() );
#endif
	C = -A*B + C;           // case 80
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A.t()*B + C;       // case 81
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A*B.t() + C;       // case 82
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A.t()*B.t() + C;   // case 83
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A*B + C;         // case 84
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A.t()*B + C;       // case 85 
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A*B.t() + C;       // case 86
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A.t()*B.t() + C;   // case 87
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A*B - C;           // case 88
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A.t()*B - C;       // case 89
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A*B.t() - C;       // case 90
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A.t()*B.t() - C;   // case 91
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A*B - C;         // case 92
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A.t()*B - C;       // case 93 
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A*B.t() - C;       // case 94
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A.t()*B.t() - C;   // case 95
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A*B + b*C;           // case 96
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A.t()*B + b*C;       // case 97
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A*B.t() + b*C;       // case 98
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A.t()*B.t() + b*C;   // case 99
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A*B + b*C;         // case 100
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A.t()*B + b*C;       // case 101 
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A*B.t() + b*C;       // case 102
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A.t()*B.t() + b*C;   // case 103
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A*B - b*C;           // case 104
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A.t()*B - b*C;       // case 105
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A*B.t() - b*C;       // case 106
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -A.t()*B.t() - b*C;   // case 107
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A*B - b*C;         // case 108
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A.t()*B - b*C;       // case 109 
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
  #include <armadillo>
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A*B.t() - b*C;       // case 110
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
	cassert ( C == CA );
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif
	C = -a*A.t()*B.t() - b*C;   // case 111
#if defined( BLASXPRS_WITH_ARMADILLO_UTIL ) && defined( BLASXPRS_TEST_AGAINST_ARMADILLO )
    CA = -a*AA.t()*BA.t() - b*CA;
	cassert ( C == CA );
#endif
#if defined( BLASXPRS_WITH_BLAZE_UTIL ) && defined( BLASXPRS_TEST_AGAINST_BLAZE )
  #include <blaze/Math.h>
#endif
#if defined( BLASXPRS_WITH_EIGEN_UTIL ) && defined( BLASXPRS_TEST_AGAINST_EIGEN )
  #include <Eigen/Dense>
#endif

*/
}



#endif    //    BLASXPRS_REAL_XGEMM_H

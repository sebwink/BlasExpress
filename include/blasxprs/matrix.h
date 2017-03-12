#ifndef BLASXPRS_MATRIX_H
#define BLASXPRS_MATRIX_H

#include <cstdlib>
#include <cstring>

#include <iostream>

#include <blasxprs/options.h>

namespace blasxprs {

#if defined( BLASXPRS_USE_CBLAS )
  using blasxprs_size_t = std::size_t; 
#elif defined( BLASXPRS_USE_F77BLAS )
  using blasxprs_size_t = int;
#endif

template <typename ValueType>
class Matrix {

private:

	blasxprs_size_t nrows;
	blasxprs_size_t ncols;
	ValueType* mtrx;

	blasxprs_size_t _stride;

public:

	using type = ValueType;

	Matrix(blasxprs_size_t num_rows, blasxprs_size_t num_cols) 
	: nrows { num_rows },
	  ncols { num_cols },
	  mtrx { (ValueType*) std::calloc(num_rows * num_cols, sizeof(ValueType)) }
	{ 
#if defined( BLASXPRS_MTRX_ROW_MAJOR )
		_stride = ncols;
#elif defined( BLASXPRS_MTRX_COLUMN_MAJOR )
		_stride = nrows;
#endif
	}

	template <typename T>
	Matrix(T* array, blasxprs_size_t num_rows, blasxprs_size_t num_cols, bool row_major) 
	: nrows { num_rows },
	  ncols { num_cols }
	{
#if defined( BLASXPRS_MTRX_ROW_MAJOR )
		_stride = ncols;
#elif defined( BLASXPRS_MTRX_COLUMN_MAJOR )
		_stride = nrows;
#endif
		mtrx = (ValueType*) malloc(nrows * ncols * sizeof(ValueType));
		blasxprs_size_t array_stride;
		if (row_major) array_stride = ncols;
		else array_stride = nrows;
		for (blasxprs_size_t i = 0; i < nrows; ++i)
			for (blasxprs_size_t j = 0; j < ncols; ++j)
				BLASXPRS_MTRX_ENTRY(i, j, mtrx, _stride) = BLASXPRS_MTRX_ENTRY(i, j, array, array_stride);
	}
	
	template <typename T>
	Matrix(Matrix<T>& other) 
	: nrows { other.num_rows() },
	  ncols { other.num_cols() },
	  mtrx { (ValueType*) std::calloc(nrows * ncols, sizeof(ValueType)) }
	{
		cpyctor(other);
    }
  	
	template<typename T>
	Matrix<ValueType>& operator=(Matrix<T>& other) {
		return cpyassign(other);
	}

    /*
	template <typename T>
	Matrix(Matrix<T>&& other) {

	}
	*/

	~Matrix() {
		free(mtrx);
	}

	ValueType& operator()(blasxprs_size_t i, blasxprs_size_t j) {
		return BLASXPRS_MTRX_ENTRY(i, j, mtrx, _stride);
	}

	const ValueType& operator()(blasxprs_size_t i, blasxprs_size_t j) const {
		return BLASXPRS_MTRX_ENTRY(i, j, mtrx, _stride);
	}

	void print(std::ostream& out) const {
		for (blasxprs_size_t i = 0; i < nrows; ++i) {
			for (blasxprs_size_t j = 0; j < ncols; ++j)
				out << BLASXPRS_MTRX_ENTRY(i, j, mtrx, _stride) << ' ';
			out << '\n';
		}
	}

	blasxprs_size_t num_rows() const{
		return nrows;
	}

	blasxprs_size_t num_cols() const {
		return ncols;
	}

	ValueType* raw() {
		return mtrx;
	}

	blasxprs_size_t stride() {
		return _stride;
	}

	// Armadillo, Blaze and Eigen utilities

	template <typename MatrixType>
	void as_matrix_type(MatrixType& M) { 
		for (blasxprs_size_t i = 0; i < nrows; ++i) 
			for (blasxprs_size_t j = 0; j < ncols; ++j)
				M(i,j) = BLASXPRS_MTRX_ENTRY(i, j, mtrx, _stride);
	}

	template <typename MatrixType>
	void from_matrix_type(const MatrixType& M) {
		for (blasxprs_size_t i = 0; i < nrows; ++i) 
			for (blasxprs_size_t j = 0; j < ncols; ++j)
				BLASXPRS_MTRX_ENTRY(i, j, mtrx, _stride) = M(i,j);
	}

#ifdef BLASXPRS_WITH_ARMADILLO_UTIL

	template<typename VT, typename T>
	friend bool operator==(const Matrix<ValueType>&, const arma::Mat<T>&);

	template<typename VT, typename T>
	friend bool operator==(const arma::Mat<T>&, const Matrix<ValueType>&);

	template <typename T>
	Matrix(const arma::Mat<T>& M) 
	: nrows { M.n_rows },
	  ncols { M.n_cols },
	  mtrx { (ValueType*) std::calloc(nrows * ncols, sizeof(ValueType)) }
	{
#if defined( BLASXPRS_MTRX_ROW_MAJOR )
		_stride = ncols;
#elif defined( BLASXPRS_MTRX_COLUMN_MAJOR )
		_stride = nrows;
#endif
		from_matrix_type(M);
	}

	template <typename T>
	Matrix<ValueType>& operator=(const arma::Mat<T> M) {
		nrows = M.n_rows();
		ncols = M.n_cols();
#if defined( BLASXPRS_MTRX_ROW_MAJOR )
		_stride = ncols;
#elif defined( BLASXPRS_MTRX_COLUMN_MAJOR )
		_stride = nrows;
#endif
		mtrx = (ValueType*) std::realloc(mtrx, nrows * ncols * sizeof(ValueType));
		from_matrix_type(M);
		return *this;
	}

	arma::Mat<ValueType> arma() {
		arma::Mat<ValueType> M(nrows, ncols);
		as_matrix_type(M);
		return M;
	}
#endif

#ifdef BLASXPRS_WITH_BLAZE_UTIL
	template <typename T>
	blaze::DynamicMatrix<T> blaze() {
		return as matrix_type<blaze::DynamicMatrix<T>>();
	}
#endif

#ifdef BLASXPRS_WITH_EIGEN_UTIL
	template <typename T>
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigen() {
		return as_matrix_type<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>();
	}
#endif

private:

	template <typename T>
	void cpyctor(Matrix<T>& other) {
#if defined( BLASXPRS_MTRX_ROW_MAJOR )
		_stride = ncols;
#elif defined( BLASXPRS_MTRX_COLUMN_MAJOR )
		_stride = nrows;
#endif
		T* other_raw = other.raw();
		for (blasxprs_size_t i = 0; i < nrows; ++i)
			for (blasxprs_size_t j = 0; j < ncols; ++j)
				BLASXPRS_MTRX_ENTRY(i, j, mtrx, _stride) = BLASXPRS_MTRX_ENTRY(i, j, other_raw, _stride);
	}

	void cpyctor(Matrix<ValueType>& other) {
#if defined( BLASXPRS_MTRX_ROW_MAJOR )
		_stride = ncols;
#elif defined( BLASXPRS_MTRX_COLUMN_MAJOR )
		_stride = nrows;
#endif
		std::memcpy(mtrx, other.raw(), nrows * ncols * sizeof(ValueType));
	}

	template <typename T>
	Matrix<ValueType>& cpyassign(Matrix<T>& other) {
		nrows = other.num_rows();
		ncols = other.num_cols();
		_stride = other.stride();
		T* other_raw = other.raw();
		mtrx = (ValueType*) std::realloc(mtrx, nrows * ncols * sizeof(ValueType));
		for (blasxprs_size_t i = 0; i < nrows; ++i)
			for (blasxprs_size_t j = 0; j < ncols; ++j)
				BLASXPRS_MTRX_ENTRY(i, j, mtrx, _stride) = BLASXPRS_MTRX_ENTRY(i, j, other_raw, _stride);
		return *this;
	}

	Matrix<ValueType>& cpyassign(Matrix<ValueType>& other) {
		nrows = other.num_rows();
		ncols = other.num_cols();
		_stride = other.stride();
		mtrx = (ValueType*) std::realloc(mtrx, nrows * ncols * sizeof(ValueType));
		std::memcpy(mtrx, other.raw(), nrows * ncols * sizeof(ValueType));
		return *this;
    }

};

#ifdef BLASXPRS_WITH_ARMADILLO_UTIL
	template<typename VT, typename T>
	bool operator==(const Matrix<VT>& M, const arma::Mat<T>& MA) {
		for (blasxprs_size_t i = 0; i < M.num_rows(); ++i)
			for (blasxprs_size_t j = 0; j < M.num_cols(); ++j)
				if ( M(i,j) != MA(i,j) )
					return false;
		return true;
	}

	template<typename VT, typename T>
	bool operator==(const arma::Mat<T>& MA, const Matrix<VT>& M) {
		return M == MA;
	}
#endif

template <typename ValueType>
void print(const Matrix<ValueType>& matrix, std::ostream& out) {
	matrix.print(out);
}

template<typename ValueType>
Matrix<ValueType> Identity(blasxprs_size_t dim) {
	Matrix<ValueType> I(dim, dim);
	for (blasxprs_size_t i = 0; i < dim; ++i)
		I(i,i) = 1.0;
	return I;         // not very efficient, implement move semantics
}

template <typename ValueType>
Matrix<ValueType> gemm(Matrix<ValueType>* Ap, 
				       Matrix<ValueType>* Bp, 
		               ValueType* alpha = nullptr,
				       Matrix<ValueType>* Cp = nullptr,
                       ValueType* beta = nullptr) {
// naive implementation
	// gemm_check_size_mismatch(A,B,C);
	Matrix<ValueType>& A = *Ap;
	Matrix<ValueType>& B = *Bp;
    // Get dimensions
	blasxprs_size_t m = A.num_rows();
	blasxprs_size_t k = A.num_cols();
	blasxprs_size_t n = B.num_cols();
	// Define the value of the alpha multiplier
	ValueType alpha_val = 1.0;
	if (alpha)
		alpha_val = *alpha;
	// Calculate alpha*A*B + beta*C
	Matrix<ValueType> Result(m,n);
	for (blasxprs_size_t i = 0; i < m; ++i)
		for (blasxprs_size_t j = 0; j < n; ++j) {
			if (Cp && beta)
				Result(i,j) += (*beta) * (*Cp)(i,j);
			for (blasxprs_size_t p = 0; p < k; ++p)
				Result(i,j) += alpha_val * A(i,p) * B(p,j);
		}
	return Result;
}

Matrix<double> blas_dgemm(Matrix<double>* Ap,
                          Matrix<double>* Bp,
					      double alpha = 1.0,
					      Matrix<double>* Cp = nullptr,
					      double beta = 0.0,
						  bool transposeA = false,
						  bool transposeB = false) {
	// extend to full dgemm functionality
	Matrix<double>& A = *Ap;
	Matrix<double>& B = *Bp;
	// dimensions M, N, K
	blasxprs_size_t m { A.num_rows() };
	blasxprs_size_t k { A.num_cols() };
	blasxprs_size_t n { B.num_cols() };
	// C
	Matrix<double> C(m,n);
	if (Cp)	C = *Cp;
    // ldx
    blasxprs_size_t lda, ldb, ldc;
	if (!transposeA && !transposeB) {
		lda = A.stride();
		ldb = B.stride();
		ldc = ldb;
	}
#if defined( BLASXPRS_USE_CBLAS )
  #ifdef BLASXPRS_WITH_AMD_ACML 
    dgemm(BLASXPRS_CBLAS_NOTRANS, BLASXPRS_CBLAS_NOTRANS, 
	      m, n, k, alpha, A.raw(), lda, B.raw(), ldb, beta, C.raw(), ldc);
  #else
	cblas_dgemm(BLASXPRS_MTRX_LAYOUT, BLASXPRS_CBLAS_NOTRANS, BLASXPRS_CBLAS_NOTRANS, 
	            m, n, k, alpha, A.raw(), lda, B.raw(), ldb, beta, C.raw(), ldc);
  #endif
#elif defined( BLASXPRS_USE_F77BLAS )
	char T = 'N';
	F77BLAS_CALL(dgemm)(&T, &T, &m, &n, &k, &alpha, A.raw(), &lda, B.raw(), &ldb, &beta, C.raw(), &ldc);
#endif
	return C;
}


template <>
Matrix<double> gemm(Matrix<double>* A,
                    Matrix<double>* B,
					double* alpha,
					Matrix<double>* C,
					double* beta) {
	double alpha_val { 1.0 };
	if (alpha) alpha_val = *alpha;
	double beta_val { 0.0 };
	if (beta) beta_val = *beta;
	return blas_dgemm(A, B, alpha_val, C, beta_val);
}

template <typename ValueType>
Matrix<ValueType> operator*(Matrix<ValueType>& A, Matrix<ValueType>& B) {
	return gemm<ValueType>(&A, &B);
}

}    //    namespace blasxprs

#endif    //    BLASXPRS_MATRIX_H

#include <iostream>

#include <blasxprs/matrix.h>

using namespace blasxprs;

int main() {
	std::cout << "\nfloat:\n " << std::endl;
	Matrix<float> If { Identity<float>(5) };
    Matrix<float> Af(5,7);
	for (int i = 0; i < 5; ++i)
		for(int j = 0; j < 7; ++j) 
			Af(i,j) = (i+1)*(j+1);
	If.print(std::cout);
	std::cout << std::endl;
	Af.print(std::cout);
	std::cout << std::endl;
	print(If*Af, std::cout);
	std::cout << "\ndouble:\n " << std::endl;
	Matrix<double> I(If);
	Matrix<double> A = Af;
	I.print(std::cout);
	std::cout << std::endl;
	A.print(std::cout);
	std::cout << std::endl;
	print(I*A, std::cout);
	std::cout << std::endl;
	print(blas_dgemm(&I,&A,10), std::cout);
	return 0;
}

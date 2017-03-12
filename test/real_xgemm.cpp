#include "real_xgemm.h"

int main() {
	xgemm_test<float>();
	xgemm_test<double>();
	return 0;
}

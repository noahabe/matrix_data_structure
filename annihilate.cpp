#include "matrix.h"

int main(void) { 
	matrix<long double> m {
		{5,1,1},
		{1,5,-1},
		{1,-1,5}
	};
	//std::function<matrix<long double>(matrix<long double>)> f = [](matrix<long double>& m) {
	auto f = [](matrix<long double>& m) {
		return (m^2) + m*((long double)-9) + (identity<long double>(m.nrow()) * (long double)18);
	};
	
	// f(m) = 0. this means f annihilates the matrix m
	std::cout << m << f(m);
	return 0;

}


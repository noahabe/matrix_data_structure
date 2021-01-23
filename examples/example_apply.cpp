#include "matrix.h"

long double square(long double x) { 
	return x * x;
}

int main(void) { 
	
	matrix<long double> m { 
		{5, 3, 1},
		{1, 2, 3},
	};
	std::cout << m;
	m.apply(square);
	std::cout << m;
	
	return 0;
} 


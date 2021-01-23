#include "matrix.h"
#include <cmath>

char NEWLINE = '\n';

int main(void) { 
	matrix<long double> m { 
		{5, 3, 1},
		{1, 2, 3},
	};
	std::cout << m;

	matrix<long double> other = m.apply<long double>(
	[](long double x){
		return std::sin(x);
	},0);
		
	std::cout << NEWLINE << m 
		<< NEWLINE << other; 
	
	return 0;
} 


#include "matrix.h"

int main(void) { 
	matrix<int> m {
		{5,1,1},
		{1,5,-1},
		{1,-1,5}
	};
	std::cout << m << std::endl 
		<< m^10 << std::endl;
	return 0;

}


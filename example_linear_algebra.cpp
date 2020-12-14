/*
Author: Noah Abe
shows some usage of 'matrix.h' 
*/

#include "matrix.h"

int main(void)
{
	matrix<double> a {
		{2,2,-1},
		{4,5,2},
		{-2,1,2}
	};
	std::vector<double> b {0,-1,2};

	std::cout << a;	
	a.append_column(b);
	std::cout << a;	

	rref(a); // rref - Reduced Row Echelon Form
	std::cout << a;
	return 0;
}

/*
Author: Noah Abe
shows some usage of 'matrix.h' 
*/

#include "matrix.h"

int main(void)
{
	matrix<double> a {
		{5,1,1},
		{1,5,-1},
		{1,-1,5}
	};
	matrix<double> p { 
		{-1,1,1},
		{1,1,0},
		{1,0,1}
	};
	matrix<double> pp (p);
	matrix<double> p_inverse = pp.inverse();
	matrix<double> d {
		{3,0,0},
		{0,6,0},
		{0,0,6}
	};

	std::cout << a << (p * d * p_inverse);	
	return 0;
}

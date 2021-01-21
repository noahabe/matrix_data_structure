#include "matrix.h"

int main(void) { 
	matrix<long double> m {
		{5,1,1},
		{1,5,-1},
		{1,-1,5}
	};
	std::cout << m
	 << (m^10);

	const matrix<long double>& m_backup = m;
	matrix<long double> mi = m_backup.inverse();
	std::cout << mi << std::endl;
	std::cout << m * mi << std::endl;
	return 0;

}


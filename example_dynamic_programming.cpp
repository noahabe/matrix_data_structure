/*
Author: Noah Abe
a typical usage of 'matrix.h' for dynamic programming. 
	[in this case subset sum problem]

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-> Problem Statement: Given an array of non-negative integers
and a positive number X, determine if there exists a subset
of the elements of the array with sum equal to X.
Example:
	INPUT ARRAY: {5, 9, 1, 8, 2} & X = 17 
	OUPUT: True 
		because 17 equals 9 + 8
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

#include <iostream>
#include <vector> 
#include <algorithm>
#include <iterator>
#include "matrix.h"

/* returns true if there exists a subset of A whose sum equals x
** 	   false otherwise
*/
bool subset_sum(const std::vector<int>& A,int x){
	int n = A.size();
	matrix<bool> MAT (n,x+1);
	
	for (int i=0;i<=n-1;++i)
		MAT[i][0] = true;
	for (int j=1;j<=x;++j){
		if (j == A[0])
			MAT[0][j] = true;
		else
			MAT[0][j] = false;	
	}
	for (int i=1;i<=n-1;++i){
		for (int j=1;j<=x;++j) {
			MAT[i][j] = MAT[i-1][j];
			if (j-A[i] >= 0)
				MAT[i][j] = MAT[i][j] or MAT[i-1][j-A[i]];
		}
	}

#ifdef _DEBUG
	std::cerr << MAT;
#endif 

	return MAT[n-1][x];
}

int main(void)
{
	std::vector<int> nums;
	int x;

	int y;
	do{
		std::cin >> y;
		nums.push_back(y);
	} while (y != -1); nums.pop_back();

	std::cout << "X = ";
	std::cin >> x;
	
	std::cout << std::boolalpha;
	std::cout << subset_sum(nums,x) << std::endl;

	return 0;
}

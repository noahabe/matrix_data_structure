// Noah Abe (data53)
/*
 * should do fine for small projects. 
 * TODO: enable the following idiom
 *
 * 	ifstream f ('m.txt');
 *  	matrix<T> m;
 *      while (f >> m) {
 * 		// do sth
 * 	}
 * there is a bug that doesn't read the last matrix.
 * TODO 2: create a variadic function that accepts N arguments
 * and multiplies them in an optimal way. // Matrix-Chain multiplication clrs page 370
 */

#include <vector>
#include <exception>
#include <iostream>
#include <iterator>
#include <initializer_list>
#include <functional>
#include <algorithm>
#include <cmath>
#include <ctype.h>

class IncompatibleMatricesError : public std::exception {
public:
	const char* what() const noexcept override {
		return "The matrices are incompatible types";
	}
};

class NotASquareMatrixError : public IncompatibleMatricesError {
public:
	const char* what() const noexcept override {
		return "The matrix needs to be a square matrix to take the power of it.";
	}
};
	
template<typename T> class matrix;
 
template<typename T>
std::ostream& operator<<(std::ostream& os,const matrix<T>& m);

template<typename T>
std::istream& operator>>(std::istream& os, matrix<T>& m);

template<typename T>
matrix<T> gaussian_elimination(matrix<T>& m,matrix<T> id = matrix<T>());

template<typename T>
matrix<T> jordan_elimination(matrix<T>& m, matrix<T> id = matrix<T>());

template<typename T>
T determinant(matrix<T>& m);

template<typename T>
class matrix{
public:
	// empty matrix with in_nor rows and in_noc columns
	explicit matrix(size_t in_nor=0,size_t in_noc=0);

	// initializer_list constructor
	matrix(std::initializer_list<std::initializer_list<T>> m);

	// copy construction
	matrix(const matrix<T>&) = default;
	matrix<T>& operator=(const matrix<T>&) = default;

	// changes the current matrix in place	
	void apply(std::function<T(T&)> f);
	// returns a whole new matrix
	// accepts an unused int parameter	
	template<typename F> 
	matrix<F> apply(std::function<F(T&)> f,int);
	
	// Add a new column, the new column must have
	// a size equal to the matrix's nor.
	void append_column(const std::vector<T>& x);

	// Matrix Transposition - returns a new matrix
	matrix<T> transpose() const; 

	// inverse of a matrix.
	// Effects on (*this): changes it to its rref
	//  (this effect is for the non const version)
	// It returns the inverse of (*this)
	matrix<T> inverse();	
	matrix<T> inverse() const;

	// returns the ith row of the matrix
	std::vector<T>& operator[](size_t i);
	const std::vector<T>& operator[](size_t i) const;

	// matrix addition
	// throws IncompatibleMatricesError if *this and rhs
	// don't have the same dimensions.
	matrix<T> operator+(const matrix<T>& rhs) const;
	matrix<T>& operator+=(const matrix<T>& rhs);

	// matrix subtraction 
	matrix<T> operator-(const matrix<T>& rhs) const;
	matrix<T>& operator-=(const matrix<T>& rhs);

	// scalar multiplication
	matrix<T> operator*(const T& k) const;	
	matrix<T>& operator*=(const T& k); 

	// matrix multiplication
	matrix<T> operator*(const matrix<T>& rhs) const;
	matrix<T>& operator*=(const matrix<T>& rhs); 

	// matrix-vector multiplication
	// b is treated as a column vector.
	std::vector<T> operator*(const std::vector<T>& b) const;
	
	// powers of a matrix.	
	matrix<T> operator^(int x) const;
	matrix<T>& operator^=(int x);
	
	bool operator==(const matrix<T>& rhs);
	bool operator!=(const matrix<T>& rhs);

	// getters
	size_t nrow() const { return nor; }
	size_t ncol() const { return noc; }

	friend std::ostream& operator<< <T>(std::ostream& os,const matrix<T>& m);
	friend std::istream& operator>> <T>(std::istream& os, matrix<T>& m); 
	friend matrix<T> gaussian_elimination<T>(matrix<T>& m,matrix<T> id);
	friend matrix<T>   jordan_elimination<T>(matrix<T>& m,matrix<T> id);
	friend T determinant<T>(matrix<T>& m);

	// a decoder
	// used primarily to create objects at run time from stored files.
	// see the friend of the class operator>> for why it is necessary
	static std::function<T(const std::string&)> str_to_T;
private:
	void init();
	std::vector<std::vector<T>> M;
	size_t nor; // Number Of Rows
	size_t noc; // Number Of Columns 
};

template<typename T>
std::function<T(const std::string&)> matrix<T>::str_to_T;

template<typename T>
std::ostream& operator<<(std::ostream& os,const matrix<T>& m){
	os << m.nor << 'x' << m.noc << '\n';
	for (size_t i=0; i < m.nor; ++i){
		for (size_t j = 0; j < m.noc; ++j)
			os << m.M[i][j] << ',';
		os << '\n';
	}
	os << "--" << '\n';
	return os;
}

template<typename T>
std::istream& operator>>(std::istream& is, matrix<T>& m){
	// read m.nor and m.noc
	/* assuming the header is in the form 
		m.nor x m.noc \n
	   the spaces don't have to exist.
	*/
	std::string buf;
	char next;
	while (is.get(next)) {
		if (next == 'x') {		
			m.nor = std::stoi(buf);
			buf = "";
			continue;
		}
		if (next == '\n') {
			m.noc = std::stoi(buf);
			break;
		}
		buf += next;
	}
	m.init();
	
	/* reads the body of the matrix */
	for (int i = 0;i<m.nor;++i) {
		for (int j=0;j<m.noc;++j) {
			buf = "";
			while (is.get(next)) {
				if (next == ',') {
					m.M[i][j] = m.str_to_T(buf);
					break;
				}
				buf += next;
			}
			 			
			if (j == m.noc-1){
				// if it is the last column read until \n
				while (is.get() != '\n'); //continue
			}	
		}
	}
	
	// read the terminating characters -- (two dashes)
	if (is.get() != '-' or is.get() != '-')
		throw std::runtime_error("Incorrect or Corrupted File read."); 

	while (isspace(is.get()));
	is.unget();

	return is;
}


template<typename T>
void matrix<T>::init(){
	M.resize(nor);
	for(auto& row : M)
		row.resize(noc);
}

template<typename T>
matrix<T>::matrix(size_t in_nor,size_t in_noc)
	: nor(in_nor),noc(in_noc) { init(); }

template<typename T>	
matrix<T>::matrix(std::initializer_list<std::initializer_list<T>> m)
{
	nor = m.size();
	noc = (*std::begin(m)).size();
	init();
	size_t counter = 0;	
	for (auto& i : m){
		M[counter++] = i;
	}
}

template<typename T>
std::vector<T>& matrix<T>::operator[](size_t i){
	return M[i];
}

template<typename T>
const std::vector<T>& matrix<T>::operator[](size_t i) const {
	return M[i];
}

template<typename T>
matrix<T> matrix<T>::operator+(const matrix<T>& rhs) const {
	if (nor != rhs.nor or noc != rhs.noc) 
		throw IncompatibleMatricesError();
	matrix<T> retval (nor,noc);
	for (size_t i=0; i<nor;++i) 
		for (size_t j=0; j<noc;++j) 
			retval.M[i][j] = M[i][j] + rhs.M[i][j];
	return retval;
}

template<typename T>
matrix<T>& matrix<T>::operator+=(const matrix<T>& rhs){
	if (nor != rhs.nor or noc != rhs.noc) 
		throw IncompatibleMatricesError();
	for (size_t i=0; i<nor;++i) 
		for (size_t j=0; j<noc;++j) 
			M[i][j] += rhs.M[i][j];
	return *this;
}

template<typename T>
matrix<T> matrix<T>::operator-(const matrix<T>& rhs) const {
	if (nor != rhs.nor or noc != rhs.noc) 
		throw IncompatibleMatricesError();
	matrix<T> retval (nor,noc);
	for (size_t i=0; i<nor;++i) 
		for (size_t j=0; j<noc;++j) 
			retval.M[i][j] = M[i][j] - rhs.M[i][j];
	return retval;
}

template<typename T>
matrix<T>& matrix<T>::operator-=(const matrix<T>& rhs){
	if (nor != rhs.nor or noc != rhs.noc) 
		throw IncompatibleMatricesError();
	for (size_t i=0; i<nor;++i) 
		for (size_t j=0; j<noc;++j) 
			M[i][j] -= rhs.M[i][j];
	return *this;
}

template<typename T>
matrix<T> matrix<T>::operator*(const T& k) const{
	matrix<T> retval (this->nor,this->noc);
	for (size_t i=0;i<nor;++i)
		for (size_t j=0;j<noc;++j)
			retval.M[i][j] = k * M[i][j];
	return retval;
}

template<typename T>
matrix<T>& matrix<T>::operator*=(const T& k){
	for (size_t i=0;i<nor;++i)
		for (size_t j=0;j<noc;++j)
			M[i][j] = k * M[i][j];
	return *this;
} 

template<typename T>
matrix<T> matrix<T>::operator*(const matrix<T>& rhs) const {
	if (this->noc != rhs.nor)
		throw IncompatibleMatricesError(); 

	matrix<T> retval (this->nor,rhs.noc);
	for (size_t i=0; i<retval.nor;++i) 
		for (size_t j=0; j<retval.noc;++j) 
			for (size_t k=0; k<this->noc;++k)
				retval.M[i][j] += M[i][k] * rhs.M[k][j];

	return retval;
}

template<typename T>
matrix<T>& matrix<T>::operator*=(const matrix<T>& rhs){
	(*this) = ((*this) * rhs);
	return (*this);
} 

template<typename T>
std::vector<T> matrix<T>::operator*(const std::vector<T>& b) const {
	if (noc != b.size())
		throw IncompatibleMatricesError();
	std::vector<T> retval (nor);
	for (int i = 0;i < nor;++i)
		for (int j = 0;j< noc; ++j)
			retval[i] += M[i][j] * b[j];
	return retval;
}


// forward declaration
template<typename T>
matrix<T> identity(int n);

template<typename T>
matrix<T> matrix<T>::operator^(int x) const {
	if (nor != noc) throw NotASquareMatrixError(); 
	if (x == 0) return identity<T>(nor);
	if (x > 0) {
		matrix<T> retval = (*this); 
		int y = std::log2(x); x -= std::pow(2,y);
		// uses successive squaring 
		while (y > 0) {
			retval = retval * retval; y--;
		} 
		// for the remaining....
		while (x > 0) {
			retval *= (*this); x--;
		}	
		return retval;
	}

//	if (x < 0) {
		x = -1 * x; // basically absolute value of x
		matrix<T> retval;
		retval = this->inverse(); // the const version
		return retval ^ x;	
//	}

}

template<typename T>
matrix<T> matrix<T>::inverse() {
	if (nor != noc)
		throw NotASquareMatrixError(); 
	return jordan_elimination(*this,gaussian_elimination(*this,identity<T>(nor))); 
}

template<typename T>
matrix<T> matrix<T>::inverse() const {
	if (nor != noc)
		throw NotASquareMatrixError(); 
	matrix<T> b = *this;
	return b.inverse(); 
}

template<typename T>
bool matrix<T>::operator==(const matrix<T>& rhs){
	if (nor != rhs.nor or noc != rhs.noc) return false; 
	for (int i=0;i<nor;++i)
		for (int j=0;j<noc;++j)
			if (M[i][j] != rhs.M[i][j])
				return false;
	return true;
}
 
template<typename T>
bool matrix<T>::operator!=(const matrix<T>& rhs){
	return not (*this == rhs);
}


template<typename T>
void matrix<T>::apply(std::function<T(T&)> f){
	for (int i=0;i<nor;++i)
		for (int j=0;j<noc;++j)
			M[i][j] = f(M[i][j]);
}

template<typename T>
template<typename F> 
matrix<F> matrix<T>::apply(std::function<F(T&)> f,int){
	matrix<F> retval (nor,noc);
	for (int i=0;i<nor;++i)
		for (int j=0;j<noc;++j)
			retval[i][j] = f(M[i][j]);
	return retval;
}
 
template<typename T>
void matrix<T>::append_column(const std::vector<T>& x){
	if (x.size() != this->nrow()) 
		throw IncompatibleMatricesError();
	for (int i = 0;i < this->nrow(); ++i) {
		M[i].push_back(x[i]);
	}
	++noc;
}

template<typename T>
matrix<T> matrix<T>::transpose() const {
	matrix<T> retval (noc,nor);
	for (int i=0;i<nor;++i)
		for (int j=0;j<noc;++j)
			retval.M[j][i] = M[i][j];
	return retval;
}

// Vector OPerationS
namespace VOPS{

// multiplying a vector by a scalar
template<typename F,typename T>
auto operator*(const F& k,const std::vector<T>& v){
	std::vector<decltype(k * v.at(0))> retval;
	retval.resize(v.size());
	for (int i=0;i<v.size();++i)
		retval[i] = k * v[i];
	return retval;		
}

// add two vectors<T> 
// the two vectors need to have the same size,
// if they don't, behavior is undefined.
template<typename T>
std::vector<T> operator+(const std::vector<T>& a,
		const std::vector<T>& b){
	std::vector<T> retval;
	retval.resize(a.size());
	for (int i=0;i<a.size();++i)
		retval[i] = a[i] + b[i];
	return retval;		
}
}// namespace VOPS

// HELPER FUNCTIONS and MATRIX ALGORITHMS
template<typename T>
matrix<T> operator*(const T& x,const matrix<T>& m){
	return m * x;
}

/*
performs a gaussian elimination on a matrix<T>.
refer "A modern introduction to Linear Algebra, Ricardo page 91"


whatever elementary operation that is done on m, it is also done on id.
this is useful to calculate the inverse of a matrix, in which case id
would be the apropriate sized identity matrix. (id for identity...)
*/
template<typename T>
matrix<T> gaussian_elimination(matrix<T>& m,matrix<T> id)
{
	const bool apply_on_id = (id.nrow() > 0 and id.ncol() > 0);

	using namespace VOPS;		
	int current_row = 0;
	for (int j=0;j<m.ncol() and current_row<m.nrow();++j) {
		// find a nonzero entry in the jth column	
		// record the row that contains this element
		// swap the row with current_row
		int pivotal_row = current_row;
		for (; pivotal_row < m.nrow(); ++pivotal_row) {
			if (m.M[pivotal_row][j] != 0)
				break;
		}
		if (pivotal_row < m.nrow()) {
			std::swap(m.M[current_row],m.M[pivotal_row]);//elementary row operation
			if (apply_on_id) {
				std::swap(id.M[current_row],id.M[pivotal_row]);
			}
			// make all the entries below this pivot equal to 0
			for (int i=current_row+1;i<m.nrow();++i) {
				if (m.M[i][j] == 0) continue;
				// elementary row operation
				auto K = (-m.M[i][j]/m.M[current_row][j]);
				m.M[i] = m.M[i] + (K * m.M[current_row]);
				if (apply_on_id) {
					id.M[i] = id.M[i] + (K * id.M[current_row]);
				}
			}
			++current_row;	
		}
	}

	return std::move(id);
}


// assumes that the matrix is already in row-echelon form
template<typename T>
matrix<T> jordan_elimination(matrix<T>& m, matrix<T> id)
{
	const bool apply_on_id = (id.nrow() > 0 and id.ncol() > 0);

	using namespace VOPS;	
	for (int i = m.nrow()-1;i >= 0; --i) {
		int pivotal_column,j;
		for (j = 0; j<m.ncol(); ++j) {
			if (m.M[i][j] != 0) {
				pivotal_column = j;
				break;	
			}
		}
		if (j == m.ncol()) {
			// the iths row was the zero vector
			continue;
		} 
		// to change the leading entry to 1
		auto K = (1/m.M[i][pivotal_column]);
		m.M[i] = K * m.M[i]; 	
		if (apply_on_id)
			id.M[i] = K * id.M[i];
		for (int higher_rows=i-1;higher_rows>=0;--higher_rows) {
			if (m.M[higher_rows][pivotal_column] == 0)
				continue;
			K = -1 * m.M[higher_rows][pivotal_column];
			m.M[higher_rows] = m.M[higher_rows] + (K * m.M[i]);
			if (apply_on_id)
				id.M[higher_rows] = id.M[higher_rows] + (K * id.M[i]);
		}
	}		
	
	return std::move(id);
}

// performs the gauss-jordan eliminatiom 
// rref - reduced row echelon form
template<typename T>
void rref(matrix<T>& m){
	gaussian_elimination(m);
	jordan_elimination(m);
}

/* some minor editions to gaussian elimination */
/* changes m to reduced echelon form as a side effect */
template<typename T>
T determinant(matrix<T>& m)
{
	if (m.nrow() != m.ncol()) 
		throw NotASquareMatrixError(); 

	T retval;
	size_t count = 0; /* number of row interchanges */

	using namespace VOPS;		
	int current_row = 0;
	for (int j=0;j<m.ncol() and current_row<m.nrow();++j) {
		// find a nonzero entry in the jth column	
		// record the row that contains this element
		// swap the row with current_row
		int pivotal_row = current_row;
		for (; pivotal_row < m.nrow(); ++pivotal_row) {
			if (m.M[pivotal_row][j] != 0)
				break;
		}
		if (pivotal_row < m.nrow()) {
			std::swap(m.M[current_row],m.M[pivotal_row]);//elementary row operation
			if (current_row != pivotal_row) {
				++count;
			}
			// make all the entries below this pivot equal to 0
			for (int i=current_row+1;i<m.nrow();++i) {
				if (m.M[i][j] == 0) continue;
				// elementary row operation
				auto K = (-m.M[i][j]/m.M[current_row][j]);
				m.M[i] = m.M[i] + (K * m.M[current_row]);
			}
			++current_row;	
		}
	}

	if (count % 2 == 0) 
		retval = 1;
	else 
		retval = -1;
	
	for (size_t i = 0;i<m.nrow();++i)
		retval *= m[i][i];
	
	return retval;
}

// used to count the number of non-zero rows
// can be used to calculate the rank of a matrix
// rank = the number of non zero rows in rref of a matrix
template<typename T>
size_t count_non_zero_rows(const matrix<T>& m){
	size_t retval = 0;
	T default_obj {};
	for (int i = 0;i < m.nrow(); ++i) {
		for (int j = 0;j < m.ncol(); ++j) {
			if (m[i][j] != default_obj) {
				// a non-zero row
				++retval;
				break;
			}		
		}
	}
	return retval;
}


// assumes the identity element is 1
// creates an identity matrix 
template<typename T>
matrix<T> identity(int n) {
	matrix<T> retval(n,n);
	T _1 {1};
	for (int i = 0;i < n;++i)
		retval[i][i] = _1;
	return retval;
}

// returns the the vandermonde matrix. Each row
// is a geometric progression of element of v
template<typename T>
matrix<T> vandermonde(const std::vector<T>& v) {
	matrix<T> retval (v.size(),v.size());
	
	for (int i = 0;i < retval.nrow(); ++i) {
		for (int j = 0; j < retval.ncol(); ++j) {
			if (j == 0) 
				retval[i][j] = 1;
			else
				retval[i][j] = retval[i][j-1] * v[i];
		}
	}

	return retval;
}

/* calculates the determinant of a vandermonde matrix
   using the method given in CLRS 3rd edition. Appendix D
   page 1226. Problem D-1 */
template<typename T>
T det_vandermonde(const std::vector<T>& v) {
	T retval {1};

	for (int j = 0;j < v.size();++j)
		for (int k = j+1;k < v.size(); ++k)
			retval *= v[k] - v[j];
	return retval;
}
	

/* calculate the trace of a matrix.
   a trace is just the sum of the diagonal elements
   of a matrix.
*/
template<typename T>
T trace(const matrix<T>& m) {
	if (m.nrow() != m.ncol()) 
		throw NotASquareMatrixError(); 
	T result {0};
	for (int i = 0; i < m.nrow(); ++i) 
		result += m[i][i];
	return result;
}


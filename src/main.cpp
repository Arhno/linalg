#include <iostream>

#include "Linalg.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

int main(int argc, char* argv[]){
	double values[] = {1.0, 2.0, 3.0} ;
	linalg::Vector<double> X(3, values) ;

	linalg::print(X) ;
	linalg::print(2.0 * X.t()) ;
	linalg::print(2.0 * X.t() * X) ;
	linalg::print(2.0 * X * X.t()) ;

	linalg::Matrix<double> B(X*X.t()) ;
	linalg::print(B) ;

	linalg::Matrix<double> M(3,3) ;
	M(0,0) = 1 ; M(1,1) = 2 ; M(2,2) = 3 ;

	linalg::print(linalg::MeanSquareSolve(M,2.0*X)) ;

	double Avalues[] = {1,2,3,4,5,6,7,8,9} ;
	linalg::Matrix<double> A(3,3,Avalues) ;
	linalg::print(A) ;

	linalg::Matrix<double> Id(3,3) ;
	Id(0,0) = Id(1,1) = Id(2,2) = 1.0 ;
	linalg::print(Id) ;
	
	linalg::print(linalg::solve(A,Id)) ;
	linalg::print(linalg::solve(A,Id)*A) ;
	linalg::print(A*linalg::solve(A,Id)) ;
	return 0 ;
}
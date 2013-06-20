#ifndef LINALG_H
#define LINALG_H

#include <cassert>

#include "MatrixExpression.h"
#include "Matrix.h"

namespace linalg {

	template<class T>
	void swapRows(Matrix<T> &A, int row1, int row2){
		if(row1 == row2)
			return ;

		int nbCols = A.nbCols() ;
		for(int i=0 ; i<nbCols ; ++i){
			T temp = A(row1, i) ;
			A(row1, i) = A(row2, i) ;
			A(row2, i) = temp ;
		}
	}

	// solve A*X=B
	template<class E1, class E2, class T>
	Matrix<T> solve(const MatrixExpression<E1,T> &ma, const MatrixExpression<E2,T> &mb){
		Matrix<T> A(ma) ;
		Matrix<T> B(mb) ;
		int n = A.nbCols() ;
		int m = B.nbCols() ;
		assert(A.nbRows() == n) ;
		assert(n == B.nbRows()) ;

		// for each column of A
		for(int i=0 ; i<n ; ++i){
			// We look for the best pivot
			int bestPivotRow = i ;
			T bestPivot = A(i,i) ;
			for(int row=i+1 ; row<n ; ++row){
				if(A(row, i) > bestPivot){
					bestPivot = A(row, i) ;
					bestPivotRow = row ;
				}
			}
			// We place the best pivot row on the ith row in both A and B
			swapRows(A, i, bestPivotRow) ;
			swapRows(B, i, bestPivotRow) ;
			// We scale the ith row
			T pinv = 1.0/bestPivot ;
			for(int col=i+1 ; col<n ; ++col)
				A(i, col) *= pinv ;
			for(int col=0 ; col<m ; ++col)
				B(i, col) *= pinv ;
			// We reduce the other rows
			for(int row=0 ; row<n ; ++row){
				if(row != i){
					T coef = A(row, i) ;
					for(int col=i+1 ; col<n ; ++col)
						A(row, col) -= coef * A(i, col) ;
					for(int col=0 ; col<m ; ++col)
						B(row, col) -= coef * B(i, col) ;
				}
			}
		}

		return B ;
	}

	template<class E1, class E2, class T>
	Matrix<T> MeanSquareSolve(const MatrixExpression<E1,T> &A, const MatrixExpression<E2,T> &B){
		return solve(A.t()*A, A.t()*B) ;
	}
}

#endif //LINALG_H
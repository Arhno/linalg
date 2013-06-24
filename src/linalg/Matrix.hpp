#ifndef MATRIX_H
#define MATRIX_H

#include <cassert>
#include <iostream>

#include "MatrixExpression.hpp"

namespace linalg {

	template<class T> class Vector ;

	template<class T>
	class Matrix : public MatrixExpression<Matrix<T>, T> {

		friend class Vector<T>;

	private:
		T** values_ ;
		int nbRows_ ;
		int nbCols_ ;

	public:
		Matrix(int nbRows, int nbCols, T* values = nullptr)
			: nbRows_(nbRows),
			  nbCols_(nbCols),
			  values_(new T*[nbRows])
		{
			int size = nbRows_*nbCols_ ;
			*values_ = new T[size] ;

			for(int i=1 ; i<nbRows_ ; ++i){
				*(values_+i) = *(values_+(i-1)) + nbCols_ ;
			}

			if(values == nullptr)
				for(int i=0 ; i<nbRows_ ; ++i)
					for(int j=0 ; j<nbCols_ ; ++j)
						values_[i][j] = 0.0 ;
			else
				for(int i=0 ; i<nbRows_ ; ++i)
					for(int j=0 ; j<nbCols_ ; ++j)
						values_[i][j] = values[i*nbCols_+j] ;
		}

		Matrix(Matrix&& m)
			: nbRows_(m.nbRows_),
			  nbCols_(m.nbCols_),
			  values_(m.values_)
		{
			m.values_ = nullptr ;
		}

		Matrix(Vector<T>&& v)
			: nbRows_(v.nbRows()),
			  nbCols_(1),
			  values_(new T*[v.nbRows()])
		{
			*values_ = v.values_ ;
			v.values_ = nullptr ;
			for(int i=1 ; i<nbRows_ ; ++i){
				*(values_+i) = *(values_+(i-1)) + 1 ;
			}
		}

		template <typename E>
		Matrix(const MatrixExpression<E,T>& me)
			: nbRows_(me.nbRows()),
			  nbCols_(me.nbCols()),
			  values_(new T*[me.nbRows()])
		{
			int size = nbRows_*nbCols_ ;
			*values_ = new T[size] ;

			for(int i=1 ; i<nbRows_ ; ++i){
				*(values_+i) = *(values_+(i-1)) + nbCols_ ;
			}

			for(int i=0 ; i<nbRows_ ; ++i)
				for(int j=0 ; j<nbCols_ ; ++j)
					values_[i][j] = me(i,j) ;
		}

		~Matrix(void){
			if(values_ != nullptr)
				delete[] (*values_) ;
			delete[] (values_) ;
		}

		int nbRows(void) const { return nbRows_; }
		int nbCols(void) const { return nbCols_; }

		// We never know if this matrix is used in the MatrixExpression, which could be problematic due to possible matrix multiplications
		template <typename E>
		void operator=(const MatrixExpression<E,T>& me) {
			*this = Matrix<T>(me) ;
		}

		void operator=(const Matrix& m) {
			if(values_ != nullptr)
				delete[] (*values_) ;
			delete[] (values_) ;

			nbRows_ = m.nbRows() ;
			nbCols_ = m.nbCols() ;

			values_ = new T*[nbRows_] ;

			int size = nbRows_*nbCols_ ;
			*values_ = new T[size] ;

			for(int i=1 ; i<nbRows_ ; ++i){
				*(values_+i) = *(values_+(i-1)) + nbCols_ ;
			}

			for(int i=0 ; i<nbRows_ ; ++i)
				for(int j=0 ; j<nbCols_ ; ++j)
					values_[i][j] = m(i,j) ;
		}

		void operator=(Matrix&& m) {
			if(values_ != nullptr)
				delete[] (*values_) ;
			delete[] (values_) ;

			nbRows_ = m.nbRows_ ;
			nbCols_ = m.nbCols_ ;
			values_ = m.values_ ;

			m.values_ = nullptr ;
		}

		void operator=(const Vector<T>& v) {
			if(values_ != nullptr)
				delete[] (*values_) ;
			delete[] (values_) ;

			nbRows_ = v.nbRows() ;
			nbCols_ = v.nbCols() ;

			values_ = new T*[nbRows_] ;

			int size = nbRows_*nbCols_ ;
			*values_ = new T[size] ;

			for(int i=1 ; i<nbRows_ ; ++i){
				*(values_+i) = *(values_+(i-1)) + nbCols_ ;
			}

			for(int i=0 ; i<nbRows_ ; ++i)
				for(int j=0 ; j<nbCols_ ; ++j)
					values_[i][j] = v(i,j) ;
		}

		void operator=(Vector<T>&& v) {
			if(values_ != nullptr)
				delete[] (*values_) ;
			delete[] (values_) ;

			nbRows_ = v.nbRows_ ;
			nbCols_ = v.nbCols_ ;

			values_ = new T*[nbRows_] ;
			*values_ = v.values_ ;
			v.values_ = nullptr ;
			for(int i=1 ; i<nbRows_ ; ++i){
				*(values_+i) = *(values_+(i-1)) + 1 ;
			}
		}

		Matrix& operator+=(const Matrix& m) {
		    assert(nbRows_ == m.nbRows()) ;
		    assert(nbCols_ == m.nbCols()) ;
		    for(int i=0 ; i<nbRows_ ; ++i)
				for(int j=0 ; j<nbCols_ ; ++j)
					values_[i][j] += m(i,j) ;
			return *this ;
		}

		template <typename E>
		Matrix& operator+=(const MatrixExpression<E,T>& me) {
			return *this += Matrix<T>(me) ;
		}

		Matrix& operator-=(const Matrix& m) {
		    assert(nbRows_ == m.nbRows()) ;
		    assert(nbCols_ == m.nbCols()) ;
		    for(int i=0 ; i<nbRows_ ; ++i)
				for(int j=0 ; j<nbCols_ ; ++j)
					values_[i][j] -= m(i,j) ;
			return *this ;
		}

		template <typename E>
		Matrix& operator-=(const MatrixExpression<E,T>& me) {
			return *this -= Matrix<T>(me) ;
		}

		const T operator()(int row, int col) const {
			return values_[row][col] ;
		}

		T& operator()(int row, int col) {
			return values_[row][col] ;
		}

		void reshape(int nbRows, int nbCols){
			nbRows_ = nbRows ;
			nbCols_ = nbCols ;
			T* temps = *values_ ;
			delete[] (values_) ;
			values_ = new T[nbRows_](nullptr) ;
			*values_ = temps ;
			for(int i=1 ; i<nbRows_ ; ++i){
				*(values_+i) = *(values_+(i-1)) + nbCols_ ;
			}
		}
	};

}

#endif //MATRIX_H

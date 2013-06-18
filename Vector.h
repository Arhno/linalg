#ifndef VECTOR_H
#define VECTOR_H

#include <cassert>
#include <iostream>
#include <vector>

#include "MatrixExpression.h"

namespace linalg {

	template<class T> class Matrix ;

	template<class T>
	class Vector : public MatrixExpression<Vector<T>, T> {

		friend class Matrix<T>;

	private:
		T* values_ ;
		int nbRows_ ;

	public:
		Vector(int nbRows, T* values = nullptr)
			: nbRows_(nbRows),
			  values_(new T[nbRows])
		{
			if(values == nullptr)
				for(int i=0 ; i<nbRows_ ; ++i)
					values_[i] = 0.0 ;
			else
				for(int i=0 ; i<nbRows_ ; ++i)
					values_[i] = values[i] ;
		}

		Vector(std::vector<T> values)
			: nbRows_(values.size()),
			  values_(new T[values.size()])
		{
			for(int i=0 ; i<nbRows_ ; ++i)
				values_[i] = values[i] ;
		}

		/*Vector(const Matrix<T>& m)
			: nbRows_(m.nbRows())
			  values_(new T[m.nbRows()])
		{
			assert(m.nbCols() == 1) ;
			for(int i=0 ; i<nbRows_ ; ++i)
				values_[i] = m(i,0) ;
		}*/

		Vector(const Matrix<T>&& m)
			: nbRows_(m.nbRows_)
		{
			assert(m.nbCols() == 1) ;
			values_ = *(m.values_) ;
			*(m.values_) = nullptr ;
		}

		Vector(Vector&& v)
			: nbRows_(v.nbRows_),
			  values_(v.values_)
		{
			v.values_ = nullptr ;
		}

		template <typename E>
		Vector(const MatrixExpression<E,T>& me)
			: nbRows_(me.nbRows()),
			  values_(new T[me.nbRows()])
		{
			assert(me.nbCols() == 1) ;
			for(int i=0 ; i<nbRows_ ; ++i)
				values_[i] = me(i,0) ;
		}

		~Vector(void){
			delete[] (values_) ;
		}
		
		int nbRows(void) const { return nbRows_; }
		int nbCols(void) const { return 1; }

		template <typename E>
		void operator=(const MatrixExpression<E,T>& me) {
			/*assert(me.nbCols() == 1) ;
			delete[] (values_) ;
			
			nbRows_ = me.nbRows() ;
			
			values_ = new T[nbRows_] ;

			for(int i=0 ; i<nbRows_ ; ++i)
				values_[i] = me(i,0) ;*/
			*this = Vector<T>(me) ;
		}

		void operator=(const Vector& v) {
			delete[] (values_) ;
			
			nbRows_ = v.nbRows_ ;
			
			values_ = new T[nbRows_] ;

			for(int i=0 ; i<nbRows_ ; ++i)
				values_[i] = v.values_[i] ;
		}

		void operator=(Vector&& v) {
			delete[] (values_) ;
			
			nbRows_ = v.nbRows_ ;
			values_ = v.values_ ; 

			v.values_ = nullptr ;
		}

		void operator=(const Matrix<T>& m) {
			assert(m.nbCols_ == 1) ;
			delete[] (values_) ;
			
			nbRows_ = m.nbRows_ ;
			
			values_ = new T[nbRows_] ;

			for(int i=0 ; i<nbRows_ ; ++i)
				values_[i] = m(i,1) ;
		}

		void operator=(Matrix<T>&& m) {
			assert(m.nbCols_ == 1) ;
			delete[] (values_) ;
			
			nbRows_ = m.nbRows_ ;
			values_ = *(v.values_) ;
			*(v.values_) = nullptr ;
		}

		const T operator()(int row, int col) const {
			return values_[row] ;
		}

		T& operator()(int row, int col) {
			return values_[row] ;
		}

		const T operator[](int i) const {
			return values_[i] ;
		}

		T& operator[](int i) {
			return values_[i] ;
		}

		operator std::vector<T>(void) const {
			std::vector<T> res(nbRows_) ;
			for(int i=0 ; i<nbRows_ ; ++i)
				res[i] = values_[i] ;

			return res ;
		}
	};

}

#endif //VECTOR_H
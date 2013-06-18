#ifndef MATRIX_EXPRESSION_H
#define MATRIX_EXPRESSION_H

#include <cassert>
#include <iostream>

namespace linalg {
	template<class T, class Op> class MatrixTranspose ;

	/*
	 * Base class for any expression involving Matrix and vector
	 */
	template<class E, class T>
	struct MatrixExpression {
		int nbRows(void) const { return static_cast<E const&>(*this).nbRows(); }
		int nbCols(void) const { return static_cast<E const&>(*this).nbCols(); }
		const T operator()(int row, int col) const { return static_cast<E const&>(*this)(row,col); }

		MatrixTranspose<T,E> const t(void) {
		   return MatrixTranspose<T,E>(static_cast<E&>(*this));
		}
 
		operator E&()             { return static_cast<      E&>(*this); }
		operator E const&() const { return static_cast<const E&>(*this); }
	};

	/*
	 * Class representing the sum of two Expressions
	 */
	template<class T, class OpLeft, class OpRight>
	class MatrixPlus : public MatrixExpression<MatrixPlus<T, OpLeft, OpRight>, T> {
	private:
		const OpLeft& left_ ;
		const OpRight& right_ ;

	public:
		MatrixPlus(const MatrixExpression<OpLeft,T>& left, const MatrixExpression<OpRight,T>& right)
			: left_(left),
			  right_(right)
		{
			assert(left.nbRows() == right.nbRows()) ;
			assert(left.nbCols() == right.nbCols()) ;
		}

		const T operator()(int row, int col) const {
			return left_(row,col) + right_(row,col) ;
		}

		int nbRows(void) const {
			return left_.nbRows() ;
		}

		int nbCols(void) const {
			return left_.nbCols() ;
		}
	};

	/*
	 * Class representing the difference of two Expressions
	 */
	template<class T, class OpLeft, class OpRight>
	class MatrixMinus : public MatrixExpression<MatrixMinus<T, OpLeft, OpRight>, T> {
	private:
		const OpLeft& left_ ;
		const OpRight& right_ ;

	public:
		MatrixMinus(const MatrixExpression<OpLeft,T>& left, const MatrixExpression<OpRight,T>& right)
			: left_(left),
			  right_(right)
		{
			assert(left.nbRows() == right.nbRows()) ;
			assert(left.nbCols() == right.nbCols()) ;
		}

		const T operator()(int row, int col) const {
			return left_(row,col) - right_(row,col) ;
		}

		int nbRows(void) const {
			return left_.nbRows() ;
		}

		int nbCols(void) const {
			return left_.nbCols() ;
		}
	};

	/*
	 * Class representing the product of two Expressions
	 */
	template<class T, class OpLeft, class OpRight>
	class MatrixTimes : public MatrixExpression<MatrixTimes<T, OpLeft, OpRight>, T> {
	private:
		const OpLeft& left_ ;
		const OpRight& right_ ;

	public:
		MatrixTimes(const MatrixExpression<OpLeft,T>& left, const MatrixExpression<OpRight,T>& right)
			: left_(left),
			  right_(right)
		{
			assert(left.nbCols() == right.nbRows()) ;
		}

		const T operator()(int row, int col) const {
			int K = left_.nbCols() ;
			T res = 0.0 ;
			for(int k=0 ; k<K ; ++k)
				res += left_(row,k) * right_(k,col) ;

			return res ;
		}

		int nbRows(void) const {
			return left_.nbRows() ;
		}

		int nbCols(void) const {
			return right_.nbCols() ;
		}
	};

	/*
	 * Class representing the scaling of an Expressions
	 */
	template<class T, class Op>
	class MatrixScale : public MatrixExpression<MatrixScale<T, Op>, T> {
	private:
		const Op& op_ ;
		const T& scale_ ;

	public:
		MatrixScale(const T& scale , const MatrixExpression<Op,T>& op)
			: op_(op),
			  scale_(scale)
		{
		}

		const T operator()(int row, int col) const {
			return scale_*op_(row,col) ;
		}

		int nbRows(void) const {
			return op_.nbRows() ;
		}

		int nbCols(void) const {
			return op_.nbCols() ;
		}
	};

	/*
	 * Class representing the transpose of an Expressions
	 */
	template<class T, class Op>
	class MatrixTranspose : public MatrixExpression<MatrixTranspose<T, Op>, T> {
	private:
		const Op& op_ ;

	public:
		MatrixTranspose(const MatrixExpression<Op,T>& op)
			: op_(op)
		{
		}

		const T operator()(int row, int col) const {
			return op_(col,row) ;
		}

		int nbRows(void) const {
			return op_.nbCols() ;
		}

		int nbCols(void) const {
			return op_.nbRows() ;
		}
	};

	/*
	 * Overload of the arithmetical operators to agregate MatixExpression's
	 */
	template <class T, class E1, class E2>
	MatrixPlus<T,E1,E2> const operator+(MatrixExpression<E1,T> const& left, MatrixExpression<E2,T> const& right) {
	   return MatrixPlus<T,E1,E2>(left,right);
	}

	template <class T, class E1, class E2>
	MatrixMinus<T,E1,E2> const operator-(MatrixExpression<E1,T> const& left, MatrixExpression<E2,T> const& right) {
	   return MatrixMinus<T,E1,E2>(left,right);
	}

	template <class T, class E1, class E2>
	MatrixTimes<T,E1,E2> const operator*(MatrixExpression<E1,T> const& left, MatrixExpression<E2,T> const& right) {
	   return MatrixTimes<T,E1,E2>(left,right);
	}

	template <class T, class E>
	MatrixScale<T,E> const operator*(T const& scale, MatrixExpression<E,T> const& op) {
	   return MatrixScale<T,E>(scale,op);
	}

	template <class T, class E>
	MatrixScale<T,E> const operator*(MatrixExpression<E,T> const& op, T const& scale) {
	   return MatrixScale<T,E>(scale,op);
	}

	template <class T, class E>
	void print(MatrixExpression<E,T> const& me, std::ostream& os = std::cout){
		for(int i=0 ; i<me.nbRows() ; ++i){
			for(int j=0 ; j<me.nbCols() ; ++j)
				os << me(i,j) << " " ;
			os << std::endl ;
		}
		os << std::endl ;
	}

}

#endif //MATRIX_EXPRESSION_H
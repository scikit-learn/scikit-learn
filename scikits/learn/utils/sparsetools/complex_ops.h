#ifndef COMPLEX_OPS_H
#define COMPLEX_OPS_H

/*
 *  Functions to handle arithmetic operations on NumPy complex values
 */

#include <numpy/arrayobject.h>
#include <iostream>

template <class c_type, class npy_type>
class complex_wrapper : public npy_type {
    template<class x, class y> 
        friend std::ostream& operator<<(std::ostream&, const complex_wrapper<x,y>& );

    public:
        complex_wrapper( const c_type r = c_type(0), const c_type i = c_type(0) ){
            npy_type::real = r;
            npy_type::imag = i;
        }
        complex_wrapper operator-() const {
          return complex_wrapper(-npy_type::real,-npy_type::imag);
        }
        complex_wrapper operator+(const complex_wrapper& B) const {
          return complex_wrapper(npy_type::real + B.real, npy_type::imag + B.imag);
        }
        complex_wrapper operator-(const complex_wrapper& B) const {
          return complex_wrapper(npy_type::real - B.real, npy_type::imag - B.imag);
        }
        complex_wrapper operator*(const complex_wrapper& B) const {
          return complex_wrapper(npy_type::real * B.real - npy_type::imag * B.imag, 
                                 npy_type::real * B.imag + npy_type::imag * B.real);
        }
        complex_wrapper operator/(const complex_wrapper& B) const {
          complex_wrapper result;
          c_type denom = 1.0 / (B.real * B.real + B.imag * B.imag);
          result.real = (npy_type::real * B.real + npy_type::imag * B.imag) * denom;
          result.imag = (npy_type::imag * B.real - npy_type::real * B.imag) * denom;
          return result;
        }
        complex_wrapper& operator+=(const complex_wrapper & B){
          npy_type::real += B.real;
          npy_type::imag += B.imag;
          return (*this);
        }
        complex_wrapper& operator-=(const complex_wrapper & B){
          npy_type::real -= B.real;
          npy_type::imag -= B.imag;
          return (*this);
        }
        complex_wrapper& operator*=(const complex_wrapper & B){
          c_type temp    = npy_type::real * B.real - npy_type::imag * B.imag;
          npy_type::imag = npy_type::real * B.imag + npy_type::imag * B.real;
          npy_type::real = temp;
          return (*this);
        }
        complex_wrapper& operator/=(const complex_wrapper & B){
          c_type denom   = 1.0 / (B.real * B.real + B.imag * B.imag);
          c_type temp    = (npy_type::real * B.real + npy_type::imag * B.imag) * denom; 
          npy_type::imag = (npy_type::imag * B.real - npy_type::real * B.imag) * denom;
          npy_type::real = temp;
          return (*this);
        }
        bool operator==(const complex_wrapper& B) const{
          return npy_type::real == B.real && npy_type::imag == B.imag;
        }
        bool operator!=(const complex_wrapper& B) const{
          return npy_type::real != B.real || npy_type::imag != B.imag;
        }
        bool operator==(const c_type& B) const{
          return npy_type::real == B && npy_type::imag == c_type(0);
        }
        bool operator!=(const c_type& B) const{
          return npy_type::real != B || npy_type::imag != c_type(0);
        }
        complex_wrapper& operator=(const complex_wrapper& B){
          npy_type::real = B.real;
          npy_type::imag = B.imag;
          return (*this);
        }
        complex_wrapper& operator=(const c_type& B){
          npy_type::real = B;
          npy_type::imag = c_type(0);
          return (*this);
        }
};

template<class x, class y> 
std::ostream& operator<<(std::ostream& out, const complex_wrapper<x,y>& cw){
    return out << cw.real << " " << cw.imag;
}

typedef complex_wrapper<float, npy_cfloat> npy_cfloat_wrapper;
typedef complex_wrapper<double,npy_cdouble> npy_cdouble_wrapper;
typedef complex_wrapper<long double,npy_clongdouble> npy_clongdouble_wrapper;


#endif

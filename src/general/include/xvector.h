/**
 * \file   xvector.h
 * This class adds a few features to the basic vector class:
 * *=, +, -, etc.
 * \author Jing Kong
 *
 */

#ifndef XVECTOR_H
#define XVECTOR_H
#include "libgen.h"
#include<vector>
#include <algorithm>    // std::transform
#include <functional>   // std::plus, minus

using namespace std;

template<class T> 
class XVector: public vector<T>
{
	public:
		XVector(size_t n, T v) : vector<T>(n,v) {}
		XVector(const XVector<T>& v) : vector<T>(v) {}
		XVector<T>& operator *= (const T& scale) 
		{ for (auto&& t: *this) t *= scale; return *this; }
		XVector<T>&  operator += (const XVector<T>& iv)
		{
		  transform(this->begin(), this->end(), iv.begin(), this->begin(), 
			          std::plus<T>());
			return *this;
		}
		XVector<T>& operator -= (const XVector<T>& iv)
		{
		  transform(this->begin(), this->end(), iv.begin(), this->begin(), 
			          std::minus<T>());
			return *this;
		}
};

template<class T> 
XVector<T> operator * (const T& scale, const XVector<T>& iv)
{
	XVector<T> rv(iv);
	rv *= scale;
	return rv;
}

template<class T>
XVector<T> operator + (const XVector<T>& iv1, const XVector<T>& iv2)
{
	XVector<T> rv(iv1);
	rv += iv2;
	return rv;
}

template<class T>
XVector<T> operator - (const XVector<T>& iv1, const XVector<T>& iv2)
{
	XVector<T> rv(iv1);
	rv -= iv2;
	return rv;
}

#endif


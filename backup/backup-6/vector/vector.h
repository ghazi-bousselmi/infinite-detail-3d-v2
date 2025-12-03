// Copyright(C) 2016  Ghazi Bousselmi <https://sites.google.com/site/ghazibousselmi>
//
// This program is free software : you can redistribute it and / or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <cassert>
#include <numeric>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <initializer_list>

using namespace std;

// a static vector
template<typename T, uint32_t N>
class Vec
{
public:
  T m_data[N];
  
public:
  typedef T* iterator;
  typedef const T* const_iterator;
  
public:  
  Vec()			{	fill(begin(), end(), (T)0);	} 
  template<typename T1>
  Vec(const Vec<T1, N>& v)	{	operator = (v);	}    
  template<typename T1>
  Vec(const initializer_list<T1>& l)
  {
    if (l.size() > N)
      copy(l.begin(), l.begin()+N, begin());
    else
    {
      copy(l.begin(), l.end(), begin());
      fill(begin()+l.size(), end(), (T)0);
    }
  }
  
  inline iterator begin() 		{	return m_data;		}
  inline iterator end() 		{	return m_data + N;	}
  inline const_iterator begin() const 	{	return m_data;		}
  inline const_iterator end() const 	{	return m_data + N;	}

  
  
  inline T operator ~ ()	{	return (T)sqrt( inner_product(begin(), end(), begin(), (T)0 ));	 }  
  inline T norm2 ()		{	return (T)inner_product(begin(), end(), begin(), (T)0 );	 }  
  inline Vec& normalize()	{	return operator /= (operator ~());	 }

  template<typename T1>
  inline Vec& operator = (const Vec<T1, N>& v)
  {
    copy(v.begin(), v.end(), begin());
  }
  
  template<typename T1>
  inline Vec& operator ^= (const Vec<T1, N>& v1)
  {
    assert(N==3);
    Vec v(*this); 
    m_data[0] = v.m_data[1] * v1.m_data[2] - v.m_data[2] * v1.m_data[1];    
    m_data[1] = v.m_data[2] * v1.m_data[0] - v.m_data[0] * v1.m_data[2];    
    m_data[2] = v.m_data[0] * v1.m_data[1] - v.m_data[1] * v1.m_data[0];    
    return *this;
  }

  template<typename T1>
  inline Vec& operator *=(const T1& t1)	  
  {	
    for_each(begin(), end(), [&t1](T& t){ t *= t1; }); 
    return *this;
  }  
  template<typename T1>
  inline Vec& operator /=(const T1& t1)	
  {	
    for_each(begin(), end(), [&t1](T& t){ t /= t1; }); 
    return *this;
  }  

  template<typename T1>
  inline Vec& operator += (const Vec<T1, N>& v)
  {	
    transform(begin(), end(), v.begin(), begin(), [](T t, const T1& t1)->T{ return t+t1; }); 
    return *this;
  }  
  template<typename T1>
  inline Vec& operator -= (const Vec<T1, N>& v)
  {	
    transform(begin(), end(), v.begin(), begin(), [](T t, const T1& t1)->T{ return t-t1; }); 
    return *this;
  }  

  template<typename T1>
  inline Vec& operator += (const T1& v)
  {	
    for_each(begin(), end(), [v](T& t){ t += v; }); 
    return *this;
  }  
  template<typename T1>
  inline Vec& operator -= (const T1& v)
  {	
    for_each(begin(), end(), [v](T& t){ t -= v; }); 
    return *this;
  }  
  
  
  template<typename T1>
  inline Vec& operator %= (const Vec<T1, N>& v)
  {	
    transform(begin(), end(), v.begin(), begin(), [](T t, const T1& t1)->T{ return t*t1; }); 
    return *this;
  }    

  inline Vec& operator >>= (uint32_t shift)
  {	
    for_each(begin(), end(), [shift](T& t){ t >>= shift; }); 
    return *this;
  }
  inline Vec& operator <<= (uint32_t shift)
  {	
    for_each(begin(), end(), [shift](T& t){ t <<= shift; }); 
    return *this;
  }

  template<typename T1>
  inline bool operator == (const Vec<T1, N>& v) const
  {	
    return equal(begin(), end(), v.begin());
  }  
  template<typename T1>
  inline bool operator < (const Vec<T1, N>& v) const
  {
    // FIXME: PERFORMANCE: this entails creation of an object std::pair
    return mismatch(begin(), end(), v.begin(), /*std::less<T>*/ [](const T& t, const T1& t1)->bool { return t < t1; }).first == end(); 
  }  
  template<typename T1>
  inline bool operator <= (const Vec<T1, N>& v) const
  {
    // FIXME: PERFORMANCE: this entails creation of an object std::pair
    return mismatch(begin(), end(), v.begin(), /*std::less_equal<T>*/ [](const T& t, const T1& t1)->bool { return t <= t1; }).first == end(); 
  }  
  template<typename T1>
  inline bool operator > (const Vec<T1, N>& v) const
  {
    // FIXME: PERFORMANCE: this entails creation of an object std::pair
    return mismatch(begin(), end(), v.begin(), /*std::greater<T>*/ [](const T& t, const T1& t1)->bool { return t > t1; }).first == end(); 
  }  
  template<typename T1>
  inline bool operator >= (const Vec<T1, N>& v) const
  {
    // FIXME: PERFORMANCE: this entails creation of an object std::pair
    return mismatch(begin(), end(), v.begin(), /*std::greater_equal<T>*/ [](const T& t, const T1& t1)->bool { return t >= t1; }).first == end(); 
  }  

  template<typename T1>
  inline Vec& operator &= (T1 v)
  {	
    for_each(begin(), end(), [v](T& t)->T{ return t &= v; }); 
    return *this;
  }
  
  
  // the min is 0
  template<typename T1>
  inline Vec& set_min_max(const Vec<T1, N>& v)
  {	
    transform(begin(), end(), v.begin(), begin(), [](T t, const T1& t1)->T{ return t < 0 ? 0 : t > t1 ? t1 : t; }); 
    return *this;
  }
  // the min is 0
  template<typename T1>
  inline Vec& set_min_max(T1 v)
  {	
    for_each(begin(), end(), [v](T& t)->T{ return t < 0 ? 0 : t > v ? v : t; }); 
    return *this;
  }
};


template<typename T1, typename T2, uint32_t N>
T1 operator * (const Vec<T1, N>& v1, const Vec<T2, N>& v2)	{	return 	inner_product(v1.begin(), v1.end(), v2.begin(), (T1)0);	}

template<typename T1, typename T2, uint32_t N>
Vec<T1, N> operator % (const Vec<T1, N>& v1, const Vec<T2, N>& v2)	{	return (Vec<T1, N>(v1) %= v2); }

template<typename T1, typename T2, uint32_t N>
Vec<T1, N> operator ^ (const Vec<T1, N>& v1, const Vec<T2, N>& v2)	{	return (Vec<T1, N>(v1) ^= v2); }

template<typename T1, typename T2, uint32_t N>
Vec<T1, N> operator + (const Vec<T1, N>& v1, const Vec<T2, N>& v2)	{	return (Vec<T1, N>(v1) += v2); }

template<typename T1, typename T2, uint32_t N>
Vec<T1, N> operator - (const Vec<T1, N>& v1, const Vec<T2, N>& v2)	{	return (Vec<T1, N>(v1) -= v2); }

template<typename T1, typename T2, uint32_t N>
Vec<T1, N> operator * (const Vec<T1, N>& v1, const T2& t2)	{	return (Vec<T1, N>(v1) *= t2); }
template<typename T1, typename T2, uint32_t N>
Vec<T2, N> operator * (const T1& t1, const Vec<T2, N>& v2)	{	return (Vec<T2, N>(v2) *= t1); }

template<typename T1, typename T2, uint32_t N>
Vec<T1, N> operator >> (const Vec<T1, N>& v1, const T2& t2)	{	return (Vec<T1, N>(v1) >>= t2); }

template<typename T1, typename T2, uint32_t N>
Vec<T1, N> operator / (const Vec<T1, N>& v1, const T2& t2)	{	return (Vec<T1, N>(v1) /= t2); }

template<typename T1, typename T2, uint32_t N>
Vec<T1, N> operator & (const Vec<T1, N>& v1, const T2& t2)	{	return (Vec<T1, N>(v1) &= t2); }


template<typename T1, typename T2, uint32_t N>
Vec<T1, N> operator + (const Vec<T1, N>& v1, const T2& v2)	{	return (Vec<T1, N>(v1) += v2); }

template<typename T1, typename T2, uint32_t N>
Vec<T1, N> operator - (const Vec<T1, N>& v1, const T2& v2)	{	return (Vec<T1, N>(v1) -= v2); }


template<typename T, uint32_t N>
ostream& operator << (ostream&s, const Vec<T, N>& v)
{
  s << "<";
  for_each(v.begin(), v.end()-1, [&s](const T& t){ s << t << ", ";});
  s << *(v.end()-1);
  s << ">";
  return s;
}

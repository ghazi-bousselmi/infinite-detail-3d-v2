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

#include "../vector/vector.h"

template<typename T, uint32_t N>
class Line
{
public:
  Vec<T, N> m_org; // origin point
  Vec<T, N> m_dir; // direction vector
  
public:
  Line(): m_org{0}, m_dir{1}	{}
  
  template<typename T1>
  explicit Line(const Line<T1, N>& l): m_org(l.m_org), m_dir(l.m_dir)	{}

  template<typename T1, typename T2>
  explicit Line(const Vec<T1, N>& org, const Vec<T2, N> dir): m_org(org), m_dir(dir)	{}
  
  inline Line& normalize() { m_dir.normalize(); return *this; }  
};

/*
 * V0
template<typename T, uint32_t N>
class Plane // hyper surface
{
public:
  Vec<T, N> m_org; // origin point
  Vec<T, N> m_dirs[N-1]; // direction vectors

public:
  Plane() 	{}
  Plane(const Plane& s): m_org(s.m_org)	{ copy(s.m_dirs, s.m_dirs+N-1, m_dirs); }
  
  inline Plane& normalize() { for_each(m_dirs, m_dirs+N-1, [](Vec<T, N>& dir) { dir.normalize(); }); return *this; }

  template<typename T1>
  inline T operator &(const Line<T1, N>& l)
  {
    // intersection, returns a float X such that : v.m_org + X * v.m_dir is the intersection, otherzise return -1 (or < 0)
    assert(N==3);
    Vec<T, N> n(m_dirs[0]);
    n ^= m_dirs[1];
    // n.normalize(); // assuming the m_dirs are normalized
    T dist = ((m_org - l.m_org) * n);
    if (dist == 0)
      return 0; // the origin of the line is on the surface
    
    Vec<T, N> proj = l.m_org + dist * n;
    Vec<T, N> vect = proj - l.m_org;
    T cosA_2 = (vect * l.m_dir);
    T norm2;
    if (cosA_2 <= 0)
      return -1; // either the line is parallel to the surface, or the direction is facing away (we do not want this case)
    cosA_2 *= cosA_2;
    cosA_2 /= (norm2 = vect.norm2());
    
    return sqrt(norm2 / cosA_2);
    
    //T distOnPlane = sqrt( (norm2 / cosA_2) * (1 - cosA_2) );    
    //Vec<T, N> inter = proj + distOnPlane * ( (l.m_dir * m_dirs[0]) * m_dirs[0] +  (l.m_dir * m_dirs[1]) * m_dirs[1] );
    //return inter;
  }
};
*/


// the directions are removed from the Plane classs for memory optimization
template<typename T, uint32_t N>
class PlaneDirs
{
public:
  Vec<T, N> m_dirs[N-1]; // direction vectors
  Vec<T, N> m_normal; 	 // direction vectors  
  
  PlaneDirs()	{}
  
  template<typename T1>
  explicit PlaneDirs(initializer_list<initializer_list<T1>> l)
  {
    auto i = l.begin();
    if (l.size() > 0)
      m_normal = *i;
    
    ++i;
    for (int j = 0; i != l.end() && j < N-1; ++i, ++j )
      m_dirs[j] = *i;
  }
  
  inline void calc_normal()
  {
    assert(N==3);
    m_normal = m_dirs[0] ^ m_dirs[1];
    m_normal.normalize();
  }
  
  inline PlaneDirs& normalize() { for_each(m_dirs, m_dirs+N-1, [](Vec<T, N>& dir) { dir.normalize(); }); return *this; }  
};




template<typename T, uint32_t N>
class Plane // hyper surface
{
public:
  Vec<T, N> 		m_org;  // origin point
  const PlaneDirs<T, N>*	m_dirs; // directions, just a copy, we do not own this pointer
public:
  Plane(): m_dirs(0)									{}
  template<typename T1>
  Plane(const Plane<T1, N>& s): m_org(s.m_org), m_dirs(s.m_dirs)			{}  
  template<typename T1, typename T2>
  Plane(const Vec<T1, N>& org, const PlaneDirs<T2, N>* dirs): m_org(org), m_dirs(dirs)	{}
  template<typename T1>
  Plane(const PlaneDirs<T1, N>* dirs): m_dirs(dirs)					{}

  inline Plane& normalize() { m_dirs->normalize(); return *this; }  

  template<typename T1>
  inline T operator &(const Line<T1, N>& l)
  {
    // intersection, returns a float X such that : v.m_org + X * v.m_dir is the intersection, otherzise return -1 (or < 0)
    T dist = ((m_org - l.m_org) * m_dirs->m_normal);
    if (dist == 0)
      return 0; // the origin of the line is on the surface
    
    //Vec<T, N> proj = l.m_org + dist * m_dirs.m_normal;
    //Vec<T, N> vect = proj - l.m_org;
    Vec<T, N> vect = dist * m_dirs->m_normal;
    T norm2;
    T cosA_2 = (vect * l.m_dir);
    if (cosA_2 <= 0)
      return -1; // either the line is parallel to the surface, or the direction is facing away (we do not want this case)
    cosA_2 *= cosA_2;
    cosA_2 /= (norm2 = vect.norm2());
    
    return sqrt(norm2 / cosA_2);
  }
};





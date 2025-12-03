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

template<typename TC = uint8_t>
class Color
{
public:
  union
  {
    uint32_t u32;
    struct
    {
      TC r;
      TC g;
      TC b;
      TC a;      
    } c;
    
    TC v[4];
  } __attribute__((packed));
  
  
  template<typename T1>
  inline void AddSaturate(const Color<T1>& o)
  {
    // this can be optimized through MMX
    c.r = (TC)((c.r + o.c.r > 255) ? 255 : (c.r + o.c.r));
    c.g = (TC)((c.g + o.c.g > 255) ? 255 : (c.g + o.c.g));
    c.b = (TC)((c.b + o.c.b > 255) ? 255 : (c.b + o.c.b));
    c.a = (TC)((c.a + o.c.a > 255) ? 255 : (c.a + o.c.a));
  }
  template<typename T1>
  inline void AddWeighted(const Color<T1>& o, double myFactor)
  {
    // this can be optimized through MMX
    c.r = (TC)(myFactor * c.r + (1.0 - myFactor) * o.c.r);
    c.g = (TC)(myFactor * c.g + (1.0 - myFactor) * o.c.g);
    c.b = (TC)(myFactor * c.b + (1.0 - myFactor) * o.c.b);
    c.a = (TC)(myFactor * c.a + (1.0 - myFactor) * o.c.a);
  }
  
  template<typename T1, typename T2>
  inline void SetWeighted(const Color<T1>& o1, const double& f1, const Color<T2>& o2, const double& f2)
  {
    // this can be optimized through MMX
    c.r = (TC)(f1 * o1.c.r + f2 * o2.c.r);
    c.g = (TC)(f1 * o1.c.g + f2 * o2.c.g);
    c.b = (TC)(f1 * o1.c.b + f2 * o2.c.b);
    c.a = (TC)(f1 * o1.c.a + f2 * o2.c.a);
  }
  
  template<typename T1>
  inline void AddAverage(const Color<T1>& o)
  {
    // this can be optimized through MMX
    c.r = (TC)((c.r + o.c.r) / 2);
    c.g = (TC)((c.g + o.c.g) / 2);
    c.b = (TC)((c.b + o.c.b) / 2);
    c.a = (TC)((c.a + o.c.a) / 2);
  }
  
  inline void Multiply(uint8_t factor)
  {
    c.r = (c.r * factor) / 255;
    c.g = (c.g * factor) / 255;
    c.b = (c.b * factor) / 255;
    c.a = (c.a * factor) / 255;
  }

  inline void Multiply(unsigned short factor)
  {
    c.r = (c.r * factor) / 65535;
    c.g = (c.g * factor) / 65535;
    c.b = (c.b * factor) / 65535;
    c.a = (c.a * factor) / 65535;
  }
  
  inline void Multiply(float factor)
  {
    c.r = (TC)(c.r * factor);
    c.g = (TC)(c.g * factor);
    c.b = (TC)(c.b * factor);
    c.a = (TC)(c.a * factor);
  }
  inline void Multiply(double factor)
  {
    c.r = (TC)(c.r * factor);
    c.g = (TC)(c.g * factor);
    c.b = (TC)(c.b * factor);
    c.a = (TC)(c.a * factor);
  }
  
  inline operator bool () const
  {
    return c.b || c.r || c.g || c.a;
  }
  
  template<typename T, typename T1, typename T2, uint32_t N>
  inline void add_child(const Color& _c, const Vec<T1, N>& __childPos, const Vec<T2, N>& __myPos, const T& __childEdgeSize, const T& __myEdgeSize)
  { 
    c.r = (TC)(_c.c.r);
    c.g = (TC)(_c.c.g);
    c.b = (TC)(_c.c.b);
    c.a = (TC)(_c.c.a);    
  }
  
  inline void SetEmpty()
  { 
    c.r = c.g = c.b = c.a = (TC)(0);    
  }  
}
__attribute__((packed));

template<typename TColor, typename TDensity = int32_t>
struct Colored
{
public:
  TColor	m_color;
  TDensity	m_density;
  
  static Colored m_empty;
  
public:
  inline Colored() { m_density = 0; m_color.SetEmpty(); }
  inline const TColor& GetColor() const  	{ return m_color; }
  inline const TColor& operator ()() const	{ return m_color; }
  inline void SetEmpty()
  { 
    m_color.SetEmpty();
    m_density = 0;
  }  
  
  // FIXME: the next will saturate very quickly with large values of m_density, and adding new colors will not change our value at all  
  inline void AddColorDensity(const TColor& c)  { m_density++; m_color.AddWeighted(c, (m_density - 1.0) / m_density);}
  Colored& operator += (const TColor& c)	{ AddColorDensity(c);  return *this; }
  Colored& operator  = (const TColor& c)	{ m_color = c; m_density = 1;  return *this; }
  template<typename TColor2, typename TDensity2>
  Colored& operator += (const Colored<TColor2, TDensity2>& o) 
  {  
    m_density += o.m_density; 
    if (m_density > 0)
      m_color.AddWeighted(o.m_color, (m_density - o.m_density + 0.0) / m_density);  
    return *this;     
  }

  template<typename TColor2, typename TDensity2>
  Colored& operator = (const Colored<TColor2, TDensity2>& o) {  m_density = o.m_density; m_color.AddWeighted(o.m_color, 0.0);  return *this; }
  
  inline operator bool () const {    return m_density;  }
}
__attribute__((packed));

template<typename TColor, typename TDensity>
Colored<TColor, TDensity> Colored<TColor, TDensity>::m_empty = Colored<TColor, TDensity>{};


template<typename T>
ostream& operator << (ostream&s, const Color<T>& c)
{
  s << "<" << c.c.r << "," << c.c.g << "," << c.c.b << "," << c.c.a << std::dec << ">";
  return s;
}
template<>
ostream& operator << <uint8_t> (ostream&s, const Color<uint8_t>& c)
{
  s << "<" << std::hex << (int)c.c.r << "," << (int)c.c.g << "," << (int)c.c.b << "," << (int)c.c.a << std::dec << ">";
  return s;
}
template<>
ostream& operator << <double> (ostream&s, const Color<double>& c)
{
  s << "<" << c.c.r << "," << c.c.g << "," << c.c.b << "," << c.c.a << ">";
  return s;
}


template<typename TColor, typename TDensity = int32_t>
ostream& operator << (ostream&s, const Colored<TColor, TDensity>& c)
{
  s << "<" << c.m_color << "," << c.m_density << ">";
  return s;
}


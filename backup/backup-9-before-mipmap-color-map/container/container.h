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

#include <iostream>
#include <memory>
#include <unordered_map>
#include "../line/line.h"
#include "color.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// replaced with a previous implementation in "color.h"
// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class Color
// {
// public:
//   union
//   {
//     uint32_t u32;
//     struct
//     {
//       uint8_t r, g, b, a;
//     };
//   } c;
// public:
//   explicit Color(uint32_t _c = 0) 		{ c.u32 = _c; }  
//   
//   template<typename T, typename T1, typename T2, uint32_t N>
//   inline void add_child(const Color& _c, const Vec<T1, N>& __childPos, const Vec<T2, N>& __myPos, const T& __childEdgeSize, const T& __myEdgeSize)
//   { c.u32 = _c.c.u32; }
//   Color& set_a(uint64_t __a) { c.a = __a; return *this; }
//   
//   Color& operator += (const Color& __c)
//   {
//     if (c.a == 0)
//       c = __c.c;
//     else if (__c.c.a != 0)
//     {
//       // FIXME: best guess: weigted average based on the A component, should be based on external counts of accumulated colors
//       double f1 = c.a / (0.0 + c.a + __c.c.a);
//       double f2 = __c.c.a / (0.0 + c.a + __c.c.a);
//       double f = f1 + f2;
//       int k;
//       
//       k = (int)(c.a * f1 + __c.c.a * f2) / f;
//       c.a = (uint8_t)(k > 255 ? 255 : k);
//       
//       k = (int)(c.b * f1 + __c.c.b * f2) / f;
//       c.b = (uint8_t)(k > 255 ? 255 : k);
// 
//       k = (int)(c.g * f1 + __c.c.g * f2) / f;
//       c.g = (uint8_t)(k > 255 ? 255 : k);
// 
//       k = (int)(c.r * f1 + __c.c.r * f2) / f;
//       c.r = (uint8_t)(k > 255 ? 255 : k);      
//     }
//     return *this;
//   }
// };
// 
// template<typename T>	inline bool operator == (T __c, const Color& __o)	{	return __c == __o.c.u32;	}
// template<typename T>	inline bool operator == (const Color& __o, T __c)	{	return __c == __o.c.u32;	}
// template<typename T>	inline bool operator != (T __c, const Color& __o)	{	return __c != __o.c.u32;	}
// template<typename T>	inline bool operator != (const Color& __o, T __c)	{	return __c != __o.c.u32;	}


template<typename T,
	 uint32_t N,				// space dimension
	 typename TColor			
	 >
class RenderStruct
{
public:
  Vec<T, N>	m_screenPoint;	// in world coordinates
  Vec<T, N>	m_LOS;		// line of sight, normalized
  T		m_curZ;		// where are we rendering now
  T		m_bestZ;	// where are we rendering now
  TColor	m_color;	// the color that we found with ray-tracing
  struct
  {
    T	m_a;	// at depth Z, solid angle radius = Z*m_a + m_b; 
    T	m_b;	 
  }		m_solid_angle_radius;

  static PlaneDirs<T, N> 	m_plane_directions[N]; // normalized

  // the next are used by 'render_compute_bounding_box'
  T		m_tmp_minZ;
  T		m_tmp_maxZ;
  Line<T, N>* 	m_tmp_adjustedLOS;
  uint32_t 	m_tmp_objectEdgeSize;
  Plane<T, N>	m_tmp_plane1;	// the m_org must be initialized with the ex-parameter "const Vec<T1, N>& objectCenterPos"
  Plane<T, N>	m_tmp_plane2;	// no need to init this

public:
  RenderStruct(): m_curZ(0), m_bestZ(1e+20), m_color(TColor()) { }
//   void init()
//   {
//     // owner must init m_screenPoint & m_LOS    
//     for (int i = 0; i < N; i++)
//     {
//       m_plane_directions[i].m_normal *= 0;
//       m_plane_directions[i].m_normal.m_data[i] = 1;
//       
//       int k = 0;
//       for (int j = 0 ; j < N ; j++)
//       {
// 	if (i == j)
// 	  continue;
// 	
// 	m_plane_directions[i].m_dirs[k] *= 0;
// 	m_plane_directions[i].m_dirs[k].m_data[j] = 1;
// 	
// 	k++;
//       }
//     }
//   }
};



template<typename T, typename TColor>
class RenderStruct<T, 3, TColor>
{
public:
  Vec<T, 3>	m_screenPoint;	// in world coordinates
  Vec<T, 3>	m_LOS;		// line of sight, normalized
  T		m_curZ;		// where are we rendering now
  T		m_bestZ;	// where are we rendering now
  TColor	m_color;	// the color that we found with ray-tracing
  struct
  {
    T	m_a;	// at depth Z, solid angle radius = Z*m_a + m_b; 
    T	m_b;	 
  }		m_solid_angle_radius;

  static PlaneDirs<T, 3> 	m_plane_directions[3]; // normalized

  // the next are used by 'render_compute_bounding_box'
  T		m_tmp_minZ;
  T		m_tmp_maxZ;
  Line<T, 3>* 	m_tmp_adjustedLOS;
  uint32_t 	m_tmp_objectEdgeSize;
  Plane<T, 3>	m_tmp_plane1;	// the m_org must be initialized with the ex-parameter "const Vec<T1, N>& objectCenterPos"
  Plane<T, 3>	m_tmp_plane2;	// no need to init this

public:
  RenderStruct(): m_curZ(0), m_bestZ(1e+20), m_color(TColor()) { }
};


template<typename T, typename TColor>
PlaneDirs<T, 3> RenderStruct<T, 3, TColor>::m_plane_directions[3] = 
{
  PlaneDirs<T, 3>({{1,0,0}, {0,1,0}, {0,0,1}}),
  PlaneDirs<T, 3>({{0,1,0}, {1,0,0}, {0,0,1}}),
  PlaneDirs<T, 3>({{0,0,1}, {1,0,0}, {0,1,0}})
};




template<typename T,
	 uint32_t N,				// space dimension
	 typename TColor			
	 >
class IContainer
{
public:
  
  typedef TColor TYPEColor;
  
  // function calls instead of direct access to members is expensive
  // even more so for virtual functions, but polymorphism is needed here for modularity
  virtual uint32_t 	get_Log2_ChildrenEdgeCount() const				= 0;
  virtual uint8_t	get_mipmap(uint32_t mipmap, uint32_t index) const 		= 0;
  virtual Colored<TColor>	get_mipmap_color(uint32_t mipmap, uint32_t index) const 	= 0;	// note this is not called intensively, thus it can be under-optimized for time & optimized for memory
  virtual uint32_t 	get_child_index(const Vec<uint32_t, N>& pos, int mipmap) const 	= 0;
  virtual uint32_t 	get_child_index(const Vec< int32_t, N>& pos, int mipmap) const 	= 0;
  virtual bool 		render_child(RenderStruct<T, N, TColor>& rs, uint32_t child_index, const Vec<T, N>& childCenterPos /* world units */, T childEdgeSize /* world units */) const = 0;
  
  virtual void		recalc_mipmap_color() = 0;
};


//class ColorConvert
//{
//public:
  template<typename T, uint32_t N, typename TColor>
  static inline Colored<typename IContainer<T, N, TColor>::TYPEColor> color_of(const typename IContainer<T, N, TColor>::TYPEColor& o)
  {
    Colored<typename IContainer<T, N, TColor>::TYPEColor> c;
    return (c+=o);
  }  

  template<typename T, uint32_t N, typename TColor>
  static inline Colored<TColor> color_of(const IContainer<T, N, TColor>* o)
  {
    return o->get_mipmap_color(o->get_Log2_ChildrenEdgeCount(), 0);
  }
  template<typename T, uint32_t N, typename TColor>
  static inline Colored<TColor> color_of(const IContainer<T, N, TColor>& o)
  {
    return o.get_mipmap_color(o.get_Log2_ChildrenEdgeCount(), 0);
  }
//};

template<typename T,	// basic corrds type
	 uint32_t N,	// space dimension
	 typename TColor,
	 typename TChild
	 >
class Container : public IContainer<T, N, TColor>
{
public:
  static constexpr uint32_t pow(uint32_t x, uint32_t n)
  {
    return n == 0 ? 1 : x * pow(x, n-1);
  }
  
protected:
  std::unordered_map<uint32_t, TChild> m_children;

  // mipmap[0] ==> not used, equivalent to m_children in size
  // mipmap[1] ==> edge-size divided by 2
  // mipmap[Log2_ChildrenEdgeCount-1] ==> 2 boxes per side
  // mipmap[Log2_ChildrenEdgeCount  ] ==> for convenience, the whole Container
  // mipmap(x, y, z) = m_children[x + y << m_Log2_ChildrenEdgeCount + z << (Log2_ChildrenEdgeCount * 2)]
  std::vector< std::vector<uint8_t> > m_mipmap;
  std::vector< std::vector<Colored<TColor> > >  m_mipmap_color;

  TChild 			m_emptyChild;  
  uint32_t 			m_Log2_ChildrenEdgeCount;
  Colored<TColor>	 	m_color;	// m_mipmap_color[m_Log2_ChildrenEdgeCount][0]
 
public:
  Container(uint32_t __Log2_ChildrenEdgeCount): m_emptyChild{}, m_Log2_ChildrenEdgeCount(__Log2_ChildrenEdgeCount)
  {
    init_mipmap();
  }
  ~Container()	
  {
  }
  
  const uint32_t m_mipmap_color_start_index = 2;
  
  
  inline void init_mipmap()
  {
    m_mipmap.resize(m_Log2_ChildrenEdgeCount + 1);
    for (int i = 1; i <= m_Log2_ChildrenEdgeCount ; i++)
    {
      m_mipmap[i].resize( pow(1u << (m_Log2_ChildrenEdgeCount-i), N) );
      fill(m_mipmap[i].begin(), m_mipmap[i].end(), 0);
    }
    
    m_mipmap_color.resize(m_Log2_ChildrenEdgeCount + 1);
    for (int i = m_mipmap_color_start_index; i <= m_Log2_ChildrenEdgeCount ; i++)
    {
      m_mipmap_color[i].resize( pow(1u << (m_Log2_ChildrenEdgeCount-i), N) );
      fill(m_mipmap_color[i].begin(), m_mipmap_color[i].end(), Colored<TColor>{});
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // all coords are in world units, as well as the sizes
  // the coords are actually the coords of the center (middle of the cube)    
  template<typename TChild_2, typename T1, typename T2>
  inline void add_child(const TChild_2& child, const Vec<T1, N>& __childPos, const Vec<T2, N>& __myPos, const T& __childEdgeSize, const T& __myEdgeSize)
  {
    //assert(__myEdgeSize >= (1u << m_Log2_ChildrenEdgeCount));
    double zoomFactor = ( __myEdgeSize * 1.0 / (1u << m_Log2_ChildrenEdgeCount) );
    Vec<T, N> childPos = ((__childPos - __myPos) + __myEdgeSize / 2.0) / zoomFactor;
    int32_t childEdgeSize = __childEdgeSize / zoomFactor;
    Vec<int32_t, N> deltaPos;
    fill(deltaPos.begin(), deltaPos.end(), -childEdgeSize/2);
    
    while (deltaPos.m_data[N-1] <= childEdgeSize/2)
    {
      // unsigned, not an issue, will overflow to the max
      Vec<int32_t, N> pos = childPos + deltaPos;
      
      bool coords_good_are_good = true;
      for_each(pos.begin(), pos.end(), [&coords_good_are_good, this](const uint32_t&x){ if (x < 0 || x >= (1u << m_Log2_ChildrenEdgeCount)) coords_good_are_good = false; });
      
      if (coords_good_are_good)
      {            
	//cout << "Container<" << typeid(T).name() << "," << N << "," << m_Log2_ChildrenEdgeCount << "," << typeid(TColor).name() << "," << typeid(TChild).name() << ">::add_child at " << __childPos << pos << (Vec<T,N>(__myPos) + pos * zoomFactor) + ( zoomFactor / 2.0 - __myEdgeSize / 2.0) << " __childEdgeSize=" << __childEdgeSize << " sub-Size=" << (T)zoomFactor << endl;  
	
	// set_child(get_child_index(pos, 0), child);
	m_children[get_child_index(pos, 0)].add_child(child, __childPos, (Vec<T,N>(__myPos) + pos * zoomFactor) + ( zoomFactor / 2.0 - __myEdgeSize / 2.0), __childEdgeSize, /*__myEdgeSize*/ (T)zoomFactor );
	
	// now update the mipmap
	Colored<TColor> color = color_of<T, N, TColor>(child);
	uint32_t index;
	for (int mipmap = 1; mipmap < m_Log2_ChildrenEdgeCount+1; mipmap++) 
	{
	  m_mipmap[mipmap][ index = get_child_index(pos, mipmap) ] = 1;
	  
	  // FIXME: this must be done after finishing the construction of the object, since the colors cant be accumulated this way, and we cannot offord a more accurate color type
	  //if (mipmap >= m_mipmap_color_start_index)
	  //  m_mipmap_color[mipmap][ index ] += color;	  
	}	
	
	// FIXME: take only the color that is inside this container, the following will take the
	// color of the full __child including the potion of it that is outside *this
	//cout << "00: m_color" << m_color << " + " << color << " = " << (decltype(m_color)(m_color) += color) << endl;
	//m_color += color;	
      }
      
      // advance the loop      
      deltaPos.m_data[0]++;
      for (int i = 0; i < N-1 && deltaPos.m_data[i] > childEdgeSize/2; i++)
      {
	deltaPos.m_data[i] = -childEdgeSize/2;
	deltaPos.m_data[i+1]++;
      }
    }
  }


  virtual uint32_t get_Log2_ChildrenEdgeCount() const
  {
    return m_Log2_ChildrenEdgeCount;
  }
  virtual uint8_t get_mipmap(uint32_t mipmap, uint32_t index) const
  {
    return m_mipmap[mipmap][index];
  }
  
  virtual Colored<TColor> get_mipmap_color(uint32_t mipmap, uint32_t index) const	// note this is not called intensively, thus it can be under-optimized for time & optimized for memory
  {
    if (m_Log2_ChildrenEdgeCount == mipmap)
      return m_color;
    
    if (mipmap == 0)
    {
      auto i = m_children.find(index);
      if (m_children.end() != i)
	return color_of<T, N, TColor>(i->second);
      return Colored<TColor>{};
    }
//     else if (mipmap == 1)
//     {
//       // this array is empty, lookup the color in mipmap=2
//       int32_t mask = (1u << (m_Log2_ChildrenEdgeCount-mipmap))-1;
//       int32_t X[N];
//       for (int i = 0; i < N; i++)
//       {
// 	X[i] = (index & mask) >> 1;	// move to mipmap 2
// 	index >>= (m_Log2_ChildrenEdgeCount-mipmap);
//       }
//       index = 0;
//       for (int i = N-1; i >= 0; i--)
//       {
// 	index <<= (m_Log2_ChildrenEdgeCount-mipmap-1);
// 	index |= X[i];
//       }      
//       return m_mipmap_color[2][index];
//     }
    else if (mipmap < m_mipmap_color_start_index)
    {
      // let us accumulate our children
      Colored<TColor> c;
      
      int32_t mask = (1u << (m_Log2_ChildrenEdgeCount-mipmap))-1;
      Vec<int32_t, N> X;
      for (int i = 0; i < N; i++)
      {
	X.m_data[i] = (index & mask) << mipmap; // move to mipmap 0
	index >>= (m_Log2_ChildrenEdgeCount-mipmap);
      }      

      Vec<int32_t, N> dX; // filled with 0s
      index = get_child_index(X, 0);
      while (dX.m_data[N-1] < (1 << mipmap))
      {
	auto i = m_children.find(index);
	if (m_children.end() != i)
	  c += color_of<T, N, TColor>(i->second);
		
	dX.m_data[0]++;
	index++;
	
	if (dX.m_data[0] >= (1 << mipmap))
	{
	  for (int i = 0; i < N-1 && dX.m_data[i] >= (1 << mipmap); i++)
	  {
	    dX.m_data[i] = 0;
	    dX.m_data[i+1]++;
	  }
	  index = get_child_index(X + dX, 0);
	} // if (loop boundary)
      } // while 
      
      return c;
    } // mip map <= 6
    return m_mipmap_color[mipmap][index];
  }
public:
  virtual void recalc_mipmap_color()
  {
    m_color = recalc_mipmap_color(m_Log2_ChildrenEdgeCount, 0);
  }
protected:
  virtual Colored<Color<double>, int64_t> recalc_mipmap_color(uint32_t mipmap, uint32_t __index)
  {
    // calculate the color [mipmap][index] and return it 
    Colored<Color<double>, int64_t> c;
    
    if (mipmap == 0)
    {
      auto i = m_children.find(__index);
      if (m_children.end() != i)
	c += color_of<T, N, TColor>(i->second);
      return c;
    } 
    
    int32_t mask = (1u << (m_Log2_ChildrenEdgeCount-mipmap))-1;
    Vec<int32_t, N> X;
    {
      auto __index_tmp = __index;
      for (int i = 0; i < N; i++)
      {
	X.m_data[i] = (__index_tmp & mask) << 1; // move to mipmap-1
	__index_tmp >>= (m_Log2_ChildrenEdgeCount-mipmap);
      } 
    }
    
    Vec<int32_t, N> dX; // filled with 0s
    auto index = get_child_index((X + dX) <<= (mipmap - 1), mipmap-1);
    
    while (dX.m_data[N-1] < (1 << 1))
    {
      if (mipmap==1 || get_mipmap(mipmap-1, index) /* no need to recalc if this sub-box does not have children */)
      {
	auto c1 = recalc_mipmap_color(mipmap - 1, index);      
	//cout << "00: m_mipmap_color["<<mipmap<<"]["<<index<<"]" << c << " + " << c1 << " = " << (decltype(c)(c) += c1) << endl;
	c += c1;
      }

      dX.m_data[0]++;
      index++;
      
      if (dX.m_data[0] >= (1 << 1))
      {
	for (int i = 0; i < N-1 && dX.m_data[i] >= (1 << 1); i++)
	{
	  dX.m_data[i] = 0;
	  dX.m_data[i+1]++;
	}
	index = get_child_index((X + dX) <<= (mipmap - 1), mipmap - 1);
      } // if (loop boundary)
    } // while
    
    if (mipmap >= m_mipmap_color_start_index)
    {
      m_mipmap_color[mipmap][__index] = c;    
    }
    return c;
  } // virtual Colored<Color<double>, int64_t> recalc_mipmap_color(uint32_t mipmap, uint32_t index)
  
public:
  virtual uint32_t get_child_index(const Vec<uint32_t, N>& pos, int mipmap) const
  {
    // FIXME: could be optimized for 3D
    uint32_t index = 0;
    for (int i = N-1; i >= 0; i-- )
      index = (index << (m_Log2_ChildrenEdgeCount-mipmap)) + ((pos.m_data[i]) >> mipmap);
    return index;
  }
  virtual uint32_t get_child_index(const Vec<int32_t, N>& pos, int mipmap) const
  {
    // FIXME: could be optimized for 3D
    int32_t index = 0;
    for (int i = N-1; i >= 0; i-- )
      index = (index << (m_Log2_ChildrenEdgeCount-mipmap)) + ((pos.m_data[i]) >> mipmap);
    return index;
  }  
  virtual bool 		render_child(RenderStruct<T, N, TColor>& rs, uint32_t child_index, const Vec<T, N>& childCenterPos /* world units */, T childEdgeSize /* world units */) const
  {
    //cout << "++IContainer::render_child()" << endl;
    auto i = m_children.find(child_index);
    if (m_children.end() == i)
    {
      //cout << "--IContainer::render_child() : result " << false << endl;
      return false;
    }
    bool result = render(rs, i->second, childCenterPos /* world units */, childEdgeSize /* world units */);
    //cout << "--IContainer::render_child() : result " << result << endl;
    return result;
  }
};

template<typename T,				// basic corrds type
	 uint32_t N,				// space dimension
	 typename TColor
	 >
class Bridge: public IContainer<T, N, TColor>
{
protected:
  class ChildInfo
  {
  public:
    Vec<T, N>	coordsDiff;
    T originalChildEdgeSize;
    T originalMyEdgeSize;
    // we do not own the pointer
    const IContainer<T, N, TColor>* child;
    // we could add rotation info as well
  };
  
  std::vector<ChildInfo> m_children;

  // Colored<TColor> m_color;	// m_mipmap_color[m_Log2_ChildrenEdgeCount][0]  
  
  
public:
  explicit Bridge() 		{ }  
  template<typename T1, typename T2>
  inline void add_child(const IContainer<T, N, TColor>* const& __child, const Vec<T1, N>& __childPos, const Vec<T2, N>& __myPos, const T& __childEdgeSize, const T& __myEdgeSize)  
  {
    ChildInfo c;
    c.originalChildEdgeSize = __childEdgeSize;
    c.originalMyEdgeSize = __myEdgeSize;
    c.coordsDiff = __childPos;
    c.coordsDiff -= __myPos;    
    c.child = __child;
    
    m_children.push_back(c);
    
    // FIXME: take only the color that is inside this container, the following will take the
    // color of the full __child including the potion of it that is outside *this
    // m_color += color_of(__child);
      
    //cout << "Bridge<" << typeid(T).name() << "," << N << "," << typeid(TColor).name() << ">::add_child at " << __childPos << " __myPos " << __myPos << " __childEdgeSize=" << __childEdgeSize << " __myEdgeSize=" << __myEdgeSize << endl;  
  }
  inline bool operator () () { return m_children.size() > 0; }
  
  virtual uint32_t 	get_Log2_ChildrenEdgeCount() const				{ return 0; }
  virtual uint8_t 	get_mipmap(uint32_t mipmap, uint32_t index) const 		{ return m_children.size() > 0; }
  virtual Colored<TColor>	get_mipmap_color(uint32_t mipmap, uint32_t index) const 	
  { 
    Colored<Color<double>, int32_t> m_color;
    for_each(m_children.begin(), m_children.end(), [&m_color](const ChildInfo& i)
    { 
      // cout << "11: m_color" << m_color << " + " << color_of(i.child) << " = " << (decltype(m_color)(m_color) += color_of(i.child)) << endl;
      m_color += color_of(i.child);
    });
    return (Colored<TColor>{}) += m_color;     
  }
  virtual uint32_t 	get_child_index(const Vec<uint32_t, N>& pos, int mipmap) const 	{ return 0; }
  virtual uint32_t 	get_child_index(const Vec< int32_t, N>& pos, int mipmap) const 	{ return 0; }
  virtual bool 		render_child(RenderStruct<T, N, TColor>& rs, uint32_t child_index, const Vec<T, N>& childCenterPos /* my center: world units */, T childEdgeSize /* my size: world units */) const
  {
    //cout << "Bridge::render_child(" << child_index << ", " << childCenterPos << ", " << childEdgeSize << ") ; bestZ=" << rs.m_bestZ << " curZ=" << rs.m_curZ << " m_children.size()=" << m_children.size() << endl;
    
    bool result = false;    
    
    for_each(m_children.begin(), m_children.end(), 
	[&, this](const ChildInfo& c)
	{ 
	  if (	rs.m_curZ < rs.m_bestZ	)
	  { 
	    //cout << "Bridge::render_child() ==> render( " << childCenterPos << " + ( " << c.coordsDiff << " * " << (childEdgeSize * 1.0 / c.originalMyEdgeSize) << " ) = " << childCenterPos + (c.coordsDiff * (childEdgeSize * 1.0 / c.originalMyEdgeSize) ) << ", " << c.originalChildEdgeSize * (childEdgeSize * 1.0 / c.originalMyEdgeSize) << ")"  << endl;
	    if ( render(rs, *c.child, (c.coordsDiff * (childEdgeSize * 1.0 / c.originalMyEdgeSize) ) + childCenterPos, c.originalChildEdgeSize * (childEdgeSize * 1.0 / c.originalMyEdgeSize) ) )
	      result = true;
	  }
	}      
    );
    //cout << "Bridge::render_child() : result " << result << endl;
    return result;
  }
  
  virtual void recalc_mipmap_color()
  {
    // nothing to do
  }
  
};


 
// FIXME: pass these parameters inside a struct  
template<typename T, uint32_t N, typename TColor>
void render_compute_bounding_box(RenderStruct<T, N, TColor>& rs) 
{
  // ASSUME: the lower-left-innerMost-(...) corner of the object is at the coods (0, 0, ...)
  
  // do not initialize here, the variables are supposed to be initialized
  // minZ = rs.m_curZ; 
  // maxZ = rs.m_bestZ;

  rs.m_tmp_plane2.m_org = rs.m_tmp_plane1.m_org;
  
  int nFailures = 0;
  for (int i = 0; i < N; i++)
  {
    if (0 == rs.m_tmp_adjustedLOS->m_dir * rs.m_plane_directions[i].m_normal)
    {
      // FIXME: this test should be factorized
      // do not continue, the LOS is parallel to this surface
      nFailures++;
      //cout << "2-render_compute_bounding_box : (adjustedLOS // surface) [" << i << "] " << endl;
      continue; 
    }
    
    rs.m_tmp_plane1.m_dirs = rs.m_tmp_plane2.m_dirs = &rs.m_plane_directions[i];
    
    //cout << "LOS.m_org " << rs.m_tmp_adjustedLOS->m_org << " ";
    //cout << "p.m_org " << rs.m_tmp_plane1.m_org;
    
    T z1 = rs.m_tmp_plane1 & *rs.m_tmp_adjustedLOS;    
    
    rs.m_tmp_plane2.m_org.m_data[i] += rs.m_tmp_objectEdgeSize;
    T z2 = rs.m_tmp_plane2 & *rs.m_tmp_adjustedLOS;    
    //cout << " p.m_org " << rs.m_tmp_plane2.m_org << endl;
    rs.m_tmp_plane2.m_org.m_data[i] = rs.m_tmp_plane1.m_org.m_data[i];
    
    
    if (z1 < 0 && z2 < 0)
    {
      //maxZ = 1e+20;  minZ = 2e+20;
      nFailures++;
      //cout << "2-render_compute_bounding_box : (z1 < 0 && z2 < 0) [" << i << "] " << z1 << " " << z2 << endl;
      continue;	// we do not reach the bounding surfaces for this dimension
    }
    
    if (z1 > z2)
    {
      if (rs.m_tmp_maxZ > z1) rs.m_tmp_maxZ = z1;
      if (rs.m_tmp_minZ < z2) rs.m_tmp_minZ = z2;
    }
    else
    {
      if (rs.m_tmp_maxZ > z2) rs.m_tmp_maxZ = z2;
      if (rs.m_tmp_minZ < z1) rs.m_tmp_minZ = z1;
    }
  }
  
  if (nFailures >= N)
  {
    rs.m_tmp_maxZ = 1e+20;
    rs.m_tmp_minZ = 2e+20;
  }  
}

template<typename T,
	 uint32_t N,				// space dimension
	 typename TColor,
	 typename T1
	 >
bool render(RenderStruct<T, N, TColor>& __rs, const TColor& __o, const Vec<T1, N>& __objectCenterPos /* world units */, T __objectEdgeSize /* world units */)
{
  if (! __o)
    return false;
  
  __rs.m_bestZ = __rs.m_curZ;
  __rs.m_color = __o;
  //cout << "TColor - Found!!" << endl;
  return true;
}

template<typename T,
	 uint32_t N,				// space dimension
	 typename TColor,
	 typename T1
	 >
bool render(RenderStruct<T, N, TColor>& __rs, const IContainer<T, N, TColor>& __o, const Vec<T1, N>& __objectCenterPos /* world units */, T __objectEdgeSize /* world units */)
{
  
   if (__o.get_Log2_ChildrenEdgeCount() == 0)
   {
    // get the lower-left-innerMost-(...) corner of the object to the coods [0, 0, 0, ...]
    Vec<T, N> screenPoint( __rs.m_screenPoint);
    screenPoint -= __objectCenterPos;
    screenPoint += __objectEdgeSize/2.0;    
    //cout << "\n00 rs.m_screenPoint " << __rs.m_screenPoint << "objectCenterPos " << __objectCenterPos << "screenPoint " << screenPoint << endl;
    
    Line<T, N> LOS(screenPoint, __rs.m_LOS);
    
    __rs.m_tmp_minZ = __rs.m_curZ;
    __rs.m_tmp_maxZ = __rs.m_bestZ;
    __rs.m_tmp_adjustedLOS = &LOS;
    __rs.m_tmp_objectEdgeSize = __objectEdgeSize;
    __rs.m_tmp_plane1.m_org *= 0;
    
    render_compute_bounding_box(__rs);
    
    if (__rs.m_tmp_minZ > __rs.m_tmp_maxZ)
      return false;
    
    return __o.render_child(__rs, 0, __objectCenterPos, __objectEdgeSize);
   } // if final
   else
  {
    register uint32_t LOG2_CHILDREN_EDGE_COUNT = __o.get_Log2_ChildrenEdgeCount();
    register uint32_t CHILDREN_EDGE_COUNT_minus_1 = (1u << LOG2_CHILDREN_EDGE_COUNT) - 1;
    
    //assert(__objectEdgeSize >= (1u << LOG2_CHILDREN_EDGE_COUNT));

    double zoomFactor = (0.0 + __objectEdgeSize) / (1u << LOG2_CHILDREN_EDGE_COUNT);

    // get the lower-left-innerMost-(...) corner of the object to the coods [0, 0, 0, ...]
    Vec<T, N> screenPoint( __rs.m_screenPoint);
    screenPoint -= __objectCenterPos;
    screenPoint += __objectEdgeSize/2.0;
    //cout << "\n11 rs.m_screenPoint " << __rs.m_screenPoint << "objectCenterPos " << __objectCenterPos << "screenPoint " << screenPoint << endl;
    //cout << "screenPoint " << screenPoint << " zoom=> " << screenPoint / zoomFactor << endl;
    
    screenPoint /= zoomFactor;
    Line<T, N> LOS(screenPoint, __rs.m_LOS / zoomFactor);    
    
    // this is a container object
    __rs.m_tmp_minZ = __rs.m_curZ /*/ zoomFactor*/;
    __rs.m_tmp_maxZ = __rs.m_bestZ /*/ zoomFactor*/;
    __rs.m_tmp_adjustedLOS = &LOS;
    __rs.m_tmp_objectEdgeSize = 1u << LOG2_CHILDREN_EDGE_COUNT;
    __rs.m_tmp_plane1.m_org *= 0;
    
    render_compute_bounding_box(__rs);
    
    T minZ = __rs.m_tmp_minZ;
    T maxZ = __rs.m_tmp_maxZ;
    if (minZ > maxZ)
      return false;
    if (0 == __o.get_mipmap(LOG2_CHILDREN_EDGE_COUNT, 0))
      return false;
    
    //cout << endl << "screenPoint " << screenPoint << " minZ " << minZ << " maxZ " << maxZ << endl;    
    
    register Vec<uint32_t, N> pos;
    register uint32_t mipmap = LOG2_CHILDREN_EDGE_COUNT;
    register T curZ = minZ;    
    register uint32_t index;
    register uint32_t previous_1_index;
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // the following is to optimize this call: 
    // o.render_child(__rs, index, Vec<T, N>(pos * zoomFactor + __objectCenterPos) + ( - __objectEdgeSize / 2.0 + zoomFactor / 2.0 ), zoomFactor)
    // ==>
    // o.render_child(__rs / zoomFactor, index, Vec<T, N>(pos) + ( (__objectCenterPos / zoomFactor) + ( - __objectEdgeSize / zoomFactor / 2.0 + 1.0 / 2.0 ) ), 1.0)
    // ==>
    // o.render_child(__rs / zoomFactor, index, Vec<T, N>(pos) + ( (__objectCenterPos / zoomFactor) + ( - (1u << LOG2_CHILDREN_EDGE_COUNT) / 2.0 + 1.0 / 2.0 ) ), 1.0)
    // ==>
    // o.render_child(__rs / zoomFactor - ( (__objectCenterPos / zoomFactor) + ( - (1u << LOG2_CHILDREN_EDGE_COUNT) / 2.0 + 1.0 / 2.0 ) ), index, Vec<T, N>(pos), 1.0)
    //
    // multiplying __rs ==> multiply the screenPoint AND the line of sight
    // aqdding to  __rs ==> adding to the screenPoint ONLY
    RenderStruct<T, N, TColor> rs_local;
    rs_local.m_LOS = Vec<double, N>(__rs.m_LOS) / zoomFactor;
    rs_local.m_screenPoint = (Vec<double, N>(__rs.m_screenPoint) / zoomFactor) - ( (Vec<double, N>(__objectCenterPos) / zoomFactor) + ( - (1 << LOG2_CHILDREN_EDGE_COUNT) / 2.0 + 1.0 / 2.0 ) );
    rs_local.m_bestZ = __rs.m_bestZ;				// the zoom is already included in the LOS vector
    rs_local.m_color = __rs.m_color;
    rs_local.m_tmp_adjustedLOS = &LOS;
    rs_local.m_solid_angle_radius = __rs.m_solid_angle_radius; 	// the zoom is already included in the LOS vector
    //cout << "rs_local.m_screenPoint=" << rs_local.m_screenPoint << ", rs_local.m_LOS" << rs_local.m_LOS << " : T=" << typeid(T).name() << endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    
    register uint32_t* p_mipmap_masks = new uint32_t[LOG2_CHILDREN_EDGE_COUNT+1];
    std::unique_ptr<uint32_t> p_mipmap_masks__;
    p_mipmap_masks__.reset(p_mipmap_masks);
    {
      for (int mipmap = 0; mipmap <= LOG2_CHILDREN_EDGE_COUNT; mipmap++)
	p_mipmap_masks[mipmap] = ((((1u << __o.get_Log2_ChildrenEdgeCount()) - 1) >> mipmap) << mipmap);
    }
    
    register T solid_angle_radius;
    
    while (	mipmap <= LOG2_CHILDREN_EDGE_COUNT
      &&	curZ <= maxZ
    )
    {
      pos = LOS.m_org + (curZ + 0.5*min(0.5, maxZ-curZ)) * LOS.m_dir; // the current pos
      //cout << endl << "Start pos " << pos;
      pos.set_min_max(CHILDREN_EDGE_COUNT_minus_1);
      //cout << pos << endl;
      
      solid_angle_radius = rs_local.m_solid_angle_radius.m_b + curZ * rs_local.m_solid_angle_radius.m_a;      
      
      index = __o.get_child_index(pos, mipmap);
      
      
      if (mipmap == 0)
      {
	// we zoomed up to the very last contained elements

	//cout << /*(int)(0 != __o.m_children[index]) <<*/ " = __o.m_children[index=" << index << "] " << " curZ "<< curZ << " " << pos << endl; 
	
	rs_local.m_curZ = curZ;
	if (	__o.render_child	(rs_local, index, pos, 1.0)	) 
	{
	  // we got it
	  //cout << "Found color at " << pos << " bestZ= " << rs_local.m_bestZ << " => " << Vec<double, N>(pos) * zoomFactor + __objectCenterPos - __objectEdgeSize/2  << endl;
	  __rs.m_curZ = curZ;
	  __rs.m_bestZ = rs_local.m_bestZ;
	  __rs.m_color = rs_local.m_color;
	  return true;
	}
      }
      else
      {
	if (0 == __o.get_mipmap(mipmap, index))
	{
	  //cout << (int)__o.get_mipmap(mipmap, index) << " = __o.m_mipmap[mipmap=" << mipmap << "][index=" << index << "] " << " curZ "<< curZ << " " << pos << endl; 
	  //cout << "mipmap-zoom-out " << mipmap;	  
	  
	  // skip this box and any upper level box that is empty
	  // zoom out to the topmost empty box that encloses this one
	  while (mipmap < LOG2_CHILDREN_EDGE_COUNT)
	  {
	    if (0 == __o.get_mipmap(mipmap+1, __o.get_child_index(pos, mipmap+1)))
	      mipmap++;
	    else
	      break;
	  }
	  //cout << " => " << mipmap << endl;
	}
	else
	{
	  // here we should zoom to find the smallest empty mipmap, or the final non-empty one
	  
	  if (solid_angle_radius >= (1u << mipmap))
	  {
	    // the solid angle here is bigger than the radius of the current box
	    // take the color of this box, and exit
	    //cout << "Found MIPMAP[" << mipmap << "==>" << (1u << mipmap) << "==>" << solid_angle_radius << "] color at " << pos << " bestZ= " << curZ << " => " << Vec<double, N>(pos) * zoomFactor + __objectCenterPos - __objectEdgeSize/2  << endl;
	    __rs.m_curZ  = curZ;
	    __rs.m_bestZ = curZ;
	    __rs.m_color = __o.get_mipmap_color(mipmap, index).m_color;
	    return true;
	  }	  

	  //cout << (int)__o.get_mipmap(mipmap, index) << " = __o.m_mipmap[mipmap=" << mipmap << "][index=" << index << "] " << " curZ "<< curZ << " " << pos << endl; 
  	  //cout << "mipmap-zoom-in " << mipmap;	  

	  while (mipmap > 0 && solid_angle_radius < (1u << mipmap) && 0 != __o.get_mipmap(mipmap, index))
	  {
	    index = __o.get_child_index(pos, --mipmap);
	  }
	  //cout << " => " << mipmap << endl;
	  
	  if (mipmap == 0)
	  {
	    rs_local.m_curZ = curZ;
	    if (	__o.render_child	(rs_local, index, pos, 1.0)	) 
	    { // we got it
	      //cout << "Found color at " << pos << " bestZ= " << rs_local.m_bestZ << " => " << Vec<double, N>(pos) * zoomFactor + __objectCenterPos - __objectEdgeSize/2  << endl;
	      __rs.m_curZ = curZ;
	      __rs.m_bestZ = rs_local.m_bestZ;
	      __rs.m_color = rs_local.m_color;
	      return true;
	    }
	    
	    previous_1_index = __o.get_child_index(pos, 1);
	  }
	} // if (0 == __o.get_mipmap(mipmap, index)) ... else
      } // if (mipmap == 0)
      
      {
	//if (!render_skip_box(__o, rs_local, LOS, pos & p_mipmap_masks[mipmap], curZ))
	//	return false;	  
	
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// now skip this box, at this mipmap level
	rs_local.m_tmp_minZ = curZ;
	rs_local.m_tmp_maxZ = 1e+20;
	rs_local.m_tmp_objectEdgeSize = 1u << mipmap;
	// (rs_local.m_tmp_plane1.m_org = pos) &= p_mipmap_masks[mipmap];
	rs_local.m_tmp_plane1.m_org.set_binary_and(pos, p_mipmap_masks[mipmap]);
	
	render_compute_bounding_box(rs_local /*, pos & p_mipmap_masks[mipmap]*/ );
	
	if (rs_local.m_tmp_minZ > rs_local.m_tmp_maxZ)
	{
	  //cout << "render() : ERROR : failed to compute bounding box" << endl;
	  return false;
	}	  
	// here we advance
	//cout ; T tmp_oldZ = curZ;
	if (curZ < rs_local.m_tmp_maxZ - 0.51)
	  curZ = rs_local.m_tmp_maxZ + 0.01 /* make sure that we will hit the next target */;
	else
	  curZ = rs_local.m_tmp_maxZ + 0.51 /* make sure that we will hit the next target, note that we migh miss an corner that we would have just barely touched */;
	//cout << "skipping box, curZ " << tmp_oldZ << " ==> " << curZ + 0.01 << " [] " << rs_local.m_tmp_minZ << " " << rs_local.m_tmp_maxZ << endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
      }

	
      if (0 == mipmap && previous_1_index != __o.get_child_index(pos, 1))
      {
	//cout << "mipmap-zoom-out " << mipmap << " ==> " << mipmap+1 << endl;
	mipmap++; // zoom out 1	  
      }
    } // while
    /**/
  }  
  return false;
}


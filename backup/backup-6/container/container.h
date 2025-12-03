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
#include <unordered_map>
#include "../line/line.h"

static int xxxx[20] = {0};

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
  
  PlaneDirs<T, N> 	m_plane_directions[N]; // normalized
  
public:
  RenderStruct()	{	init();		}
  void init()
  {
    m_curZ = 0;
    m_bestZ = 1e+20;
    m_color = TColor();
    
    // owner must init m_screenPoint & m_LOS
    
    for (int i = 0; i < N; i++)
    {
      m_plane_directions[i].m_normal *= 0;
      m_plane_directions[i].m_normal.m_data[i] = 1;
      
      int k = 0;
      for (int j = 0 ; j < N ; j++)
      {
	if (i == j)
	  continue;
	
	m_plane_directions[i].m_dirs[k] *= 0;
	m_plane_directions[i].m_dirs[k].m_data[j] = 1;
	
	k++;
      }
    }
  }
};



template<typename T,
	 uint32_t N,				// space dimension
	 typename TColor			
	 >
class IContainer
{
public:
  // function calls instead of direct access to members is expensive
  // even more so for virtual functions, but polymorphism is needed here for modularity
  virtual uint32_t 	get_Log2_ChildrenEdgeCount() const				= 0;
  virtual char 		get_mipmap(uint32_t mipmap, uint32_t index) const 		= 0;
  virtual uint32_t 	get_child_index(const Vec<uint32_t, N>& pos, int mipmap) const 	= 0;
  virtual uint32_t 	get_child_index(const Vec< int32_t, N>& pos, int mipmap) const 	= 0;
  virtual bool 		render_child(RenderStruct<T, N, TColor>& rs, uint32_t child_index, const Vec<T, N>& childCenterPos /* world units */, T childEdgeSize /* world units */) const = 0;
};


template<typename T,	// basic corrds type
	 uint32_t N,	// space dimension
	 uint32_t Log2_ChildrenEdgeCount, // children count = (2^Log2_ChildrenEdgeCount) ^ N
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
  
  static constexpr uint32_t TOTAL_CHILDREN_COULD = pow(1u << Log2_ChildrenEdgeCount, N);
  
protected:
  // version 1: full memory usage, optimize speed 
  // child(x, y, z) = m_children[x + y << Log2_ChildrenEdgeCount + z << (Log2_ChildrenEdgeCount * 2)]
  //Vec<Vec<TChild, 1u << Log2_ChildrenEdgeCount>*, pow(1u << Log2_ChildrenEdgeCount, N-1)> m_children;
  std::unordered_map<uint32_t, TChild> m_children;
  TChild m_emptyChild;
  
protected:
  // mipmap[0] ==> not used, equivalent to m_children in size
  // mipmap[1] ==> edge-size divided by 2
  // mipmap[Log2_ChildrenEdgeCount-1] ==> 2 boxes per side
  // mipmap[Log2_ChildrenEdgeCount  ] ==> for convenience, the whole Container
  // coords access uses the same convention as for m_children
  std::vector<uint8_t /*bool*/> m_mipmap[Log2_ChildrenEdgeCount+1];
  
public:
  // mipmap max coords
  Vec<uint32_t, N> m_max_coords[Log2_ChildrenEdgeCount+1];

public:
  Container(): m_emptyChild{}
  {
    init_mipmap();
  }
  ~Container()	
  {
  }
  
  
  inline void init_mipmap()
  {
    for (int i = 1; i <= Log2_ChildrenEdgeCount ; i++)
    {
      m_mipmap[i].resize( pow(1u << (Log2_ChildrenEdgeCount-i), N) );
      fill(m_mipmap[i].begin(), m_mipmap[i].end(), 0);
    }
    for (int i = 0; i <= Log2_ChildrenEdgeCount ; i++)
    {
      fill(m_max_coords[i].begin(), m_max_coords[i].end(), (1u << i) - 1);
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // all coords are in world units, as well as the sizes
  // the coords are actually the coords of the center (middle of the cube)    
  template<typename TChild_2, typename T1, typename T2>
  inline void add_child(const TChild_2& child, const Vec<T1, N>& __childPos, const Vec<T2, N>& __myPos, const T& __childEdgeSize, const T& __myEdgeSize)
  {
    //assert(__myEdgeSize >= (1u << Log2_ChildrenEdgeCount));
    double zoomFactor = ( __myEdgeSize * 1.0 / (1u << Log2_ChildrenEdgeCount) );
    Vec<T, N> childPos = ((__childPos - __myPos) + __myEdgeSize / 2.0) / zoomFactor;
    int32_t childEdgeSize = __childEdgeSize / zoomFactor;
    Vec<int32_t, N> deltaPos;
    fill(deltaPos.begin(), deltaPos.end(), -childEdgeSize/2);
    
    while (deltaPos.m_data[N-1] <= childEdgeSize/2)
    {
      // unsigned, not an issue, will overflow to the max
      Vec<int32_t, N> pos = childPos + deltaPos;
      
      bool coords_good_are_good = true;
      for_each(pos.begin(), pos.end(), [&coords_good_are_good](const uint32_t&x){ if (x < 0 || x >= (1u << Log2_ChildrenEdgeCount)) coords_good_are_good = false; });
      
      if (coords_good_are_good)
      {            
	// set_child(get_child_index(pos, 0), child);
	m_children[get_child_index(pos, 0)].add_child(child, __childPos, (pos * zoomFactor + __myPos) + ( zoomFactor / 2.0 - __myEdgeSize / 2.0), __childEdgeSize, /*__myEdgeSize*/ (T)zoomFactor );
	
	// now update the mipmap
	////cout << pos << endl;
	for (int mipmap = 1; mipmap < Log2_ChildrenEdgeCount+1; mipmap++) 
	{
	  if (0 == m_mipmap[mipmap][ get_child_index(pos, mipmap) ])
	  {
	    xxxx[mipmap]++;
	  }
	  m_mipmap[mipmap][ get_child_index(pos, mipmap) ] = 1;
	}
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
    return Log2_ChildrenEdgeCount;
  }
  virtual char get_mipmap(uint32_t mipmap, uint32_t index) const
  {
    return m_mipmap[mipmap][index];
  }
  virtual uint32_t get_child_index(const Vec<uint32_t, N>& pos, int mipmap) const
  {
    // FIXME: could be optimized for 3D
    uint32_t index = 0;
    for (int i = N-1; i >= 0; i-- )
      index = (index << (Log2_ChildrenEdgeCount-mipmap)) + (pos.m_data[i] >> mipmap);
    return index;
  }
  virtual uint32_t get_child_index(const Vec<int32_t, N>& pos, int mipmap) const
  {
    // FIXME: could be optimized for 3D
    int32_t index = 0;
    for (int i = N-1; i >= 0; i-- )
      index = (index << (Log2_ChildrenEdgeCount-mipmap)) + (pos.m_data[i] >> mipmap);
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


class Color
{
public:
  union
  {
    uint32_t u32;
    struct
    {
      uint8_t r, g, b, a;
    };
  } c;
public:
  explicit Color(uint32_t _c = 0) 		{ c.u32 = _c; }  
  
  template<typename T, typename T1, typename T2, uint32_t N>
  inline void add_child(const Color& _c, const Vec<T1, N>& __childPos, const Vec<T2, N>& __myPos, const T& __childEdgeSize, const T& __myEdgeSize)
  { c.u32 = _c.c.u32; }
};

template<typename T>	inline bool operator == (T c, const Color& o)	{	return c == o.c.u32;	}
template<typename T>	inline bool operator == (const Color& o, T c)	{	return c == o.c.u32;	}
template<typename T>	inline bool operator != (T c, const Color& o)	{	return c != o.c.u32;	}
template<typename T>	inline bool operator != (const Color& o, T c)	{	return c != o.c.u32;	}

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
  }
  inline bool operator () () { return m_children.size() > 0; }
  
  virtual uint32_t 	get_Log2_ChildrenEdgeCount() const				{ return 0; }
  virtual char 		get_mipmap(uint32_t mipmap, uint32_t index) const 		{ return m_children.size() > 0; }
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
	    //cout << "Bridge::render_child() ==> render( " << childCenterPos << " + ( " << c.coordsDiff << " * " << (c.originalMyEdgeSize * 1.0 / childEdgeSize) << " ) = " << childCenterPos + (c.coordsDiff * (c.originalMyEdgeSize * 1.0 / childEdgeSize) ) << ", " << c.originalChildEdgeSize * (c.originalMyEdgeSize * 1.0 / childEdgeSize) << ")"  << endl;
	    if ( render(rs, *c.child, childCenterPos + (c.coordsDiff * (c.originalMyEdgeSize * 1.0 / childEdgeSize) ), c.originalChildEdgeSize * (c.originalMyEdgeSize * 1.0 / childEdgeSize) ) )
	      result = true;
	  }
	}      
    );
    
    //cout << "Bridge::render_child() : result " << result << endl;
    return result;
  }
};


 
// FIXME: pass these parameters inside a struct  
template<typename T, uint32_t N, typename TColor, typename T1>
void render_compute_bounding_box(const RenderStruct<T, N, TColor>& rs, T& minZ, T& maxZ, const Line<T, N>& adjustedLOS, uint32_t objectEdgeSize, const Vec<T1, N>& objectCenterPos)
{
  // ASSUME: the lower-left-innerMost-(...) corner of the object is at the coods "objectCenterPos"
  
  // do not initialize here, the variables are supposed to be initialized
  // minZ = rs.m_curZ; 
  // maxZ = rs.m_bestZ;
  
  int nFailures = 0;
  for (int i = 0; i < N; i++)
  {
    if (0 == adjustedLOS.m_dir * rs.m_plane_directions[i].m_normal)
    {
      // FIXME: this test should be factorized
      // do not continue, the LOS is parallel to this surface
      nFailures++;
      //cout << "2-render_compute_bounding_box : (adjustedLOS // surface) [" << i << "] " << endl;
      continue; 
    }

    Plane<T, N> p(objectCenterPos, &rs.m_plane_directions[i]);
    
    //cout << "LOS.m_org " << adjustedLOS.m_org << " ";
    //cout << "p.m_org " << p.m_org;
    
    T z1 = p & adjustedLOS;
    
    p.m_org.m_data[i] += objectEdgeSize;
    T z2 = p & adjustedLOS;
    
    //cout << " p.m_org " << p.m_org << endl;
    
    if (z1 < 0 && z2 < 0)
    {
      //maxZ = 1e+20;  minZ = 2e+20;
      nFailures++;
      //cout << "2-render_compute_bounding_box : (z1 < 0 && z2 < 0) [" << i << "] " << z1 << " " << z2 << endl;
      continue;	// we do not reach the bounding surfaces for this dimension
    }
    
    if (z1 > z2)
    {
      if (maxZ > z1) maxZ = z1;
      if (minZ < z2) minZ = z2;
    }
    else
    {
      if (maxZ > z2) maxZ = z2;
      if (minZ < z1) minZ = z1;
    }
  }
  
  if (nFailures >= N)
  {
    maxZ = 1e+20;
    minZ = 2e+20;
  }  
}

template<typename T, uint32_t N, typename TColor>
void render_compute_bounding_box(const RenderStruct<T, N, TColor>& rs, T& minZ, T& maxZ, const Line<T, N>& adjustedLOS, uint32_t objectEdgeSize)
{
  // ASSUME: the lower-left-innerMost-(...) corner of the object is at the coods [0, 0, 0, ...]
  render_compute_bounding_box(rs, minZ, maxZ, adjustedLOS, objectEdgeSize, Vec<uint32_t, N>());
}



// FIXME: pass these parameters inside a structure
template<typename T,
	 uint32_t N,				// space dimension
	 class TObject,				// the object to render
	 typename TColor			
	 >
inline bool render_skip_box(const TObject& o, RenderStruct<T, N, TColor>& rs, const Line<T, N>& LOS, uint32_t mipmap, register Vec<uint32_t, N> pos, T& curZ)
{
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // now skip this box, at this mipmap level
  T __minZ = curZ, __maxZ = 1e+20;
  render_compute_bounding_box(rs, __minZ, __maxZ, LOS, 1u << mipmap, /*(pos >> mipmap) << mipmap*/ pos & ((((1u << o.get_Log2_ChildrenEdgeCount()) - 1) >> mipmap) << mipmap) );
  if (__minZ > __maxZ)
  {
    //cout << "render() : ERROR : failed to compute bounding box" << endl;
    return false;
  }	  
  // here we advance
  T oldZ = curZ;
  if (curZ < __maxZ - 0.51)
    curZ = __maxZ + 0.01 /* make sure that we will hit the next target */;
  else
    curZ = __maxZ + 0.51 /* make sure that we will hit the next target, note that we migh miss an corner that we would have just barely touched */;
  //cout << "skipping box, curZ " << oldZ << " ==> " << curZ + 0.01 << " [] " << __minZ << " " << __maxZ << endl;
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  return true;
}

template<typename T,
	 uint32_t N,				// space dimension
	 typename TColor			
	 >
bool render(RenderStruct<T, N, TColor>& rs, const TColor& o, const Vec<T, N>& objectCenterPos /* world units */, T objectEdgeSize /* world units */)
{
  if (0 == o)
    return false;
  
  rs.m_bestZ = rs.m_curZ;
  rs.m_color = o;
  //cout << "TColor - Found!!" << endl;
  return true;
}




template<typename T,
	 uint32_t N,				// space dimension
	 typename TColor			
	 >
bool render(RenderStruct<T, N, TColor>& rs, const IContainer<T, N, TColor>& o, const Vec<T, N>& objectCenterPos /* world units */, T objectEdgeSize /* world units */)
{
  
   if (o.get_Log2_ChildrenEdgeCount() == 0)
   {
     return o.render_child(rs, 0, objectCenterPos, objectEdgeSize);
   } // if final
   else
  {
    register uint32_t LOG2_CHILDREN_EDGE_COUNT = o.get_Log2_ChildrenEdgeCount();
    register uint32_t CHILDREN_EDGE_COUNT_minus_1 = (1u << LOG2_CHILDREN_EDGE_COUNT) - 1;
    
    assert(objectEdgeSize >= (1u << LOG2_CHILDREN_EDGE_COUNT));

    double zoomFactor = (0.0 + objectEdgeSize) / (1u << LOG2_CHILDREN_EDGE_COUNT);

    // get the lower-left-innerMost-(...) corner of the object to the coods [0, 0, 0, ...]
    Vec<T, N> screenPoint = rs.m_screenPoint - objectCenterPos + objectEdgeSize/2.0;
    //cout << "rs.m_screenPoint " << rs.m_screenPoint << "objectCenterPos " << objectCenterPos << "screenPoint " << screenPoint << endl;
    //cout << "screenPoint " << screenPoint << " zoom=> " << screenPoint / zoomFactor << endl;
    
    screenPoint /= zoomFactor;
    Line<T, N> LOS(screenPoint, rs.m_LOS / zoomFactor);
    
    
    /**/
    // this is a container object
    T minZ = rs.m_curZ / zoomFactor;
    T maxZ = rs.m_bestZ / zoomFactor;
    render_compute_bounding_box(rs, minZ, maxZ, LOS, 1u << LOG2_CHILDREN_EDGE_COUNT);
    if (minZ > maxZ)
      return false;
    if (0 == o.get_mipmap(LOG2_CHILDREN_EDGE_COUNT, 0))
      return false;
    
    register Vec<uint32_t, N> pos;

    //cout << endl << "screenPoint " << screenPoint << " minZ " << minZ << " maxZ " << maxZ << endl;

    uint32_t mipmap = LOG2_CHILDREN_EDGE_COUNT;
    register T curZ = minZ;
    
    register uint32_t index;
    register uint32_t previous_1_index;
    while (	mipmap <= LOG2_CHILDREN_EDGE_COUNT
      &&	curZ <= maxZ
    )
    {
      pos = LOS.m_org + (curZ + 0.5*min(0.5, maxZ-curZ)) * LOS.m_dir; // the current pos
      //cout << endl << "Start pos " << pos;
      pos.set_min_max(CHILDREN_EDGE_COUNT_minus_1);
      //cout << pos << endl;
      
      
      index = o.get_child_index(pos, mipmap);

      if (mipmap == 0)
      {
	//cout << /*(int)(0 != o.m_children[index]) <<*/ " = o.m_children[index=" << index << "] " << " curZ "<< curZ << " " << pos << endl; 
	
	// we zoomed up to the very last contained elements
	rs.m_curZ = curZ * zoomFactor;	    
	if (	o.render_child	(rs, index, Vec<T, N>(pos*zoomFactor + objectCenterPos) + ( - objectEdgeSize / 2.0 + zoomFactor / 2.0 ), zoomFactor)	) 
	{
	  // we got it
	  //cout << "Found color at " << pos << " => " << Vec<double, N>(pos) * zoomFactor + objectCenterPos - objectEdgeSize/2  << endl;
	  return true;
	}
      }
      else
      {
	if (0 == o.get_mipmap(mipmap, index))
	{
	  //cout << (int)o.get_mipmap(mipmap, index) << " = o.m_mipmap[mipmap=" << mipmap << "][index=" << index << "] " << " curZ "<< curZ << " " << pos << endl; 
	  //cout << "mipmap-zoom-out " << mipmap;	  
	  // skip this box and any upper level box that is empty
	  // zoom out to the topmost empty box that encloses this one
	  while (mipmap < LOG2_CHILDREN_EDGE_COUNT)
	  {
	    if (0 == o.get_mipmap(mipmap+1, o.get_child_index(pos, mipmap+1)) )
	      mipmap++;
	    else
	      break;	    
	  }
	  //cout << " => " << mipmap << endl;
	}
	else
	{
	  //cout << (int)o.get_mipmap(mipmap, index) << " = o.m_mipmap[mipmap=" << mipmap << "][index=" << index << "] " << " curZ "<< curZ << " " << pos << endl; 
  	  //cout << "mipmap-zoom-in " << mipmap;	  

	  // here we should zoom to find the smallest non empty mipmap
	  while (mipmap > 0 && 0 != o.get_mipmap(mipmap, index))
	  {
	    index = o.get_child_index(pos, --mipmap);
	  }
	  //cout << " => " << mipmap << endl;
	  
	  if (mipmap == 0)
	  {
	    rs.m_curZ = curZ * zoomFactor;	    
	    if (	o.render_child	(rs, index, Vec<T, N>(pos*zoomFactor + objectCenterPos) + ( - objectEdgeSize / 2.0 + zoomFactor / 2.0 ), zoomFactor)		) 
	    { // we got it
	      //cout << "Found color at " << pos << " => " << Vec<double, N>(pos) * zoomFactor + objectCenterPos - objectEdgeSize/2  << endl;
	      return true;
	    }	    
	    previous_1_index = o.get_child_index(pos, 1);
	  }
	} // if (0 == o.get_mipmap(mipmap, index)) ... else
      } // if (mipmap == 0)
      
      if (!render_skip_box(o, rs, LOS, mipmap, pos, curZ))
	return false;	  
	
      if (0 == mipmap && previous_1_index != o.get_child_index(pos, 1))
      {
	//cout << "mipmap-zoom-out " << mipmap << " ==> " << mipmap+1 << endl;
	mipmap++; // zoom out 1	  
      }
    } // while
    /**/
  }  
  return false;
}


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
#include "../line/line.h"

#ifndef __in
#	define	__in
#endif
#ifndef __out
#	define	__out
#endif
#ifndef __inout
#	define	__inout
#endif
#ifndef __tmp
#	define	__tmp
#endif
#ifndef __const
#	define	__const
#endif


template<typename T,	// basic corrds type
	 uint32_t N,	// space dimension
	 uint32_t Log2_ChildrenEdgeCount, // children count = (2^Log2_ChildrenEdgeCount) ^ N
	 typename TChild
	 >
class Container
{
public:
  static constexpr uint32_t pow(uint32_t x, uint32_t n)
  {
    return n == 0 ? 1 : x * pow(x, n-1);
  }
  
  static constexpr uint32_t TOTAL_CHILDREN_COULD = pow(1u << Log2_ChildrenEdgeCount, N);
  enum
  {
    LOG2_CHILDREN_EDGE_COUNT = Log2_ChildrenEdgeCount
  };
  
public:
  // version 1: full memory usage, optimize speed 
  // child(x, y, z) = m_children[x + y << Log2_ChildrenEdgeCount + z << (Log2_ChildrenEdgeCount * 2)]
  TChild m_children[TOTAL_CHILDREN_COULD];
  
  // version 2: optimize memory usage  
  // N-1 dimensions are static
  // std::vector<TChild>	m_children[pow(1u << Log2_ChildrenEdgeCount, N-1)];
  
  /*bool*/ 
  // mipmap[0] ==> not used, equivalent to m_children in size
  // mipmap[1] ==> edge-size divided by 2
  // mipmap[Log2_ChildrenEdgeCount-1] ==> 2 boxes per side
  // mipmap[Log2_ChildrenEdgeCount  ] ==> for convenience, the whole Container
  // coords access uses the same convention as for m_children
protected:
  std::vector<uint8_t> m_mipmap[Log2_ChildrenEdgeCount+1];

public:  
  // mipmap max coords
  Vec<uint32_t, N> m_max_coords[Log2_ChildrenEdgeCount+1];
  
  inline char get_mipmap(uint32_t mipmap, uint32_t index) const
  {
    if (m_mipmap[mipmap].size()==0)
      return 0;
    return m_mipmap[mipmap][index >> 3] & (1 << (index & 7));
  }
  inline void set_mipmap(uint32_t mipmap, uint32_t index)
  {
    // [0] is not allowed
    // !!!!! if (mipmap == 0) return; !!!!!
    
    if (m_mipmap[mipmap].size()==0)
    {
      m_mipmap[mipmap].resize( 1 + pow(1u << (Log2_ChildrenEdgeCount-mipmap), N) / 8 );
      fill(m_mipmap[mipmap].begin(), m_mipmap[mipmap].end(), 0);
    }    
    m_mipmap[mipmap][index >> 3] |= (1 << (index & 7));
  }

public:
  Container()	
  {
    init_mipmap();
  }
  
  
  inline void init_mipmap()
  {
//     for (int i = 1; i <= Log2_ChildrenEdgeCount ; i++)
//     {
//       m_mipmap[i].resize( pow(1u << (Log2_ChildrenEdgeCount-i), N) );
//       fill(m_mipmap[i].begin(), m_mipmap[i].end(), 0);
//     }
    for (int i = 0; i <= Log2_ChildrenEdgeCount ; i++)
    {
      fill(m_max_coords[i].begin(), m_max_coords[i].end(), (1u << i) - 1);
    }
  }
  
  inline uint32_t get_child_index(const Vec<uint32_t, N>& pos, int mipmap) const
  {
    // FIXME: could be optimized for 3D
    uint32_t index = 0;
    for (int i = N-1; i >= 0; i-- )
      index = (index << (Log2_ChildrenEdgeCount-mipmap)) + (pos.m_data[i] >> mipmap);
    return index;
  }
  inline uint32_t get_child_index(const Vec<int32_t, N>& pos, int mipmap) const
  {
    // FIXME: could be optimized for 3D
    uint32_t index = 0;
    for (int i = N-1; i >= 0; i-- )
      index = (index << (Log2_ChildrenEdgeCount-mipmap)) + (pos.m_data[i] >> mipmap);
    return index;
  }
  
  // the coords are in my coords system, i.e. from 0 to [(1u << Log2_ChildrenEdgeCount) - 1]
  inline void add_child(const Vec<uint32_t, N>& pos, const TChild& child)
  {
    // FIXME: could be optimized for 3D
    for_each(pos.begin(), pos.end(), [](const uint32_t&x){ assert(x >= 0 && x < (1u << Log2_ChildrenEdgeCount)); });
    
    m_children[ get_child_index(pos, 0) ] = child;
    
    // now update the mipmap
    for (int mipmap = 1; mipmap < Log2_ChildrenEdgeCount+1; mipmap++) 
      set_mipmap(mipmap, get_child_index(pos, mipmap));
  }
  
};


class Color
{
public:
  enum
  {
    LOG2_CHILDREN_EDGE_COUNT = 0
  };
  
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
  explicit Color(uint32_t _c = 0) { c.u32 = _c; }
};

template<typename T>	inline bool operator == (T c, const Color& o)	{	return c == o.c.u32;	}
template<typename T>	inline bool operator == (const Color& o, T c)	{	return c == o.c.u32;	}
template<typename T>	inline bool operator != (T c, const Color& o)	{	return c != o.c.u32;	}
template<typename T>	inline bool operator != (const Color& o, T c)	{	return c != o.c.u32;	}



template<typename T, uint32_t N, typename TColor>
class RenderStruct
{
public:
  Vec<T, N>	m_screenPoint;	// in world coordinates
  Vec<T, N>	m_LOS;		// line of sight, normalized
  T		m_maxZ;		// absolute maximum Z
  T		m_curZ;		// where are we rendering now
  T		m_bestZ;	// where are we rendering now
  TColor	m_color;	// the color that we found with ray-tracing
  
  PlaneDirs<T, N> 	m_plane_directions[N]; // normalized
  
public:
  RenderStruct()	{	init();		}
  void init()
  {
    m_curZ = 0;
    m_maxZ = m_bestZ = 1e+20;
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

template<typename T, uint32_t N, typename TColor, typename T1 = uint32_t>
class BoundingBoxParams
{
public:
  __in const RenderStruct<T, N, TColor>* 	m_rs;
  __in const Line<T, N>* 			m_adjustedLOS;
  
  __inout T		m_curZ;			// the current Z within the rendering loop    
  __in    uint32_t 	m_mipmap;		// input for render_skip_box
  
  __tmp   T		m_minZ;			// temporary for render_skip_box, output for render_compute_bounding_box
  __tmp   T		m_maxZ;			// temporary for render_skip_box, output for render_compute_bounding_box
  
  __tmp uint32_t 	m_objectEdgeSize;	// temporary for render_skip_box, input for render_compute_bounding_box
  __tmp Vec<T1, N> 	m_objectCenterPos;	// temporary for render_skip_box, input for render_compute_bounding_box
  
  __const uint32_t 	m_children_count_mask; 	// ((1u << TObject::LOG2_CHILDREN_EDGE_COUNT) - 1)
}; // currently in rendering loop


// FIXME: pass these parameters inside a struct  
template<typename T, uint32_t N, typename TColor, typename T1>
void render_compute_bounding_box(BoundingBoxParams<T, N, TColor, T1>& bbp)
{
  // ASSUME: the lower-left-innerMost-(...) corner of the object is at the coods "objectCenterPos"
  
  // do not initialize here, the variables are supposed to be initialized
  // bbp.m_minZ = rs.m_curZ; 
  // bbp.m_maxZ = rs.m_maxZ;
  
  int nFailures = 0;
  for (int i = 0; i < N; i++)
  {
    if (0 == bbp.m_adjustedLOS->m_dir * bbp.m_rs->m_plane_directions[i].m_normal)
    {
      // FIXME: this test should be factorized
      // do not continue, the LOS is parallel to this surface
      nFailures++;
      //cout << "2-render_compute_bounding_box : (adjustedLOS // surface) [" << i << "] " << endl;
      continue; 
    }

    Plane<T, N> p(bbp.m_objectCenterPos, &bbp.m_rs->m_plane_directions[i]);
    
    //cout << "LOS.m_org " << adjustedLOS.m_org << " ";
    //cout << "p.m_org " << p.m_org;
    
    T z1 = p & *bbp.m_adjustedLOS;
    
    p.m_org.m_data[i] += bbp.m_objectEdgeSize;
    T z2 = p & *bbp.m_adjustedLOS;
    
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
      if (bbp.m_maxZ > z1) bbp.m_maxZ = z1;
      if (bbp.m_minZ < z2) bbp.m_minZ = z2;
    }
    else
    {
      if (bbp.m_maxZ > z2) bbp.m_maxZ = z2;
      if (bbp.m_minZ < z1) bbp.m_minZ = z1;
    }
  }
  
  if (nFailures >= N)
  {
    bbp.m_maxZ = 1e+20;
    bbp.m_minZ = 2e+20;
  }  
}

// FIXME: pass these parameters inside a struct  
template<typename T, uint32_t N, typename TColor>
void render_compute_bounding_box(const RenderStruct<T, N, TColor>& rs, T& minZ, T& maxZ, const Line<T, N>& adjustedLOS, uint32_t objectEdgeSize)
{
  // ASSUME: the lower-left-innerMost-(...) corner of the object is at the coods [0, 0, 0, ...]
  BoundingBoxParams<T, N, TColor> bbp;
  bbp.m_rs = &rs;
  bbp.m_adjustedLOS = &adjustedLOS;
  bbp.m_objectEdgeSize = objectEdgeSize;
  bbp.m_objectCenterPos = Vec<uint32_t, N>();  
  bbp.m_minZ = minZ;
  bbp.m_maxZ = maxZ;
  render_compute_bounding_box(bbp);
  
  minZ = bbp.m_minZ;
  maxZ = bbp.m_maxZ;
}




// FIXME: pass these parameters inside a structure
template<typename T,
	 uint32_t N,				// space dimension
	 typename TColor			// must define an enumeration { LOG2_CHILDREN_EDGE_COUNT }
	 >
inline bool render_skip_box(BoundingBoxParams<T, N, TColor>& bbp, const Vec<uint32_t, N>& pos)
{
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // now skip this box, at this mipmap level
  bbp.m_minZ = bbp.m_curZ;
  bbp.m_maxZ = 1e+20;
  bbp.m_objectEdgeSize = 1u << bbp.m_mipmap;
  bbp.m_objectCenterPos = pos & ((bbp.m_children_count_mask >> bbp.m_mipmap) << bbp.m_mipmap);
  
  render_compute_bounding_box(bbp);
  if (bbp.m_minZ > bbp.m_maxZ)
  {
    //cout << "render() : ERROR : failed to compute bounding box" << endl;
    return false;
  }	  
  // here we advance
  T oldZ = bbp.m_curZ;
  if (bbp.m_curZ < bbp.m_maxZ - 0.51)
    bbp.m_curZ = bbp.m_maxZ + 0.01 /* make sure that we will hit the next target */;
  else
    bbp.m_curZ = bbp.m_maxZ + 0.51 /* make sure that we will hit the next target, note that we migh miss an corner that we would have just barely touched */;
  //cout << "skipping box, curZ " << oldZ << " ==> " << curZ + 0.01 << " [] " << __minZ << " " << __maxZ << endl;
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  return true;
}


template<typename T,
	 uint32_t N,				// space dimension
	 class TObject,				// the object to render
	 typename TColor			// must define an enumeration { LOG2_CHILDREN_EDGE_COUNT }
	 >
void render(RenderStruct<T, N, TColor>& rs, const TObject& o, const Vec<T, N>& objectCenterPos /* world units */, uint32_t objectEdgeSize /* world units */)
{
  double zoomFactor = (0.0 + objectEdgeSize) / (1u << TObject::LOG2_CHILDREN_EDGE_COUNT);

  // get the lower-left-innerMost-(...) corner of the object to the coods [0, 0, 0, ...]
  Vec<T, N> screenPoint = rs.m_screenPoint - objectCenterPos + objectEdgeSize/2.0;
  //cout << "rs.m_screenPoint " << rs.m_screenPoint << "objectCenterPos " << objectCenterPos << "screenPoint " << screenPoint << endl;
  //cout << "screenPoint " << screenPoint << " zoom=> " << screenPoint / zoomFactor << endl;
  
  screenPoint /= zoomFactor;
  Line<T, N> LOS(screenPoint, rs.m_LOS / zoomFactor);
  
  if (TObject::LOG2_CHILDREN_EDGE_COUNT == 0)
  {
    /*
    // this is a final object, a color point
    // let us check if the line of sight hits it
    T minZ = rs.m_curZ;
    T maxZ = rs.m_maxZ;
    render_compute_bounding_box(rs, minZ, maxZ, LOS, objectEdgeSize);
    if (minZ <= maxZ)
    {
      if (minZ <= rs.m_bestZ)
      {
	rs.m_bestZ = minZ;
	rs.m_color = o;
      }
    }
    else
    {
      // the LOS does not hit this object
    }
    */
  } // if final
  else
  {
    assert(objectEdgeSize >= (1u << TObject::LOG2_CHILDREN_EDGE_COUNT));
    
    
    /**/
    // this is a container object
    T minZ = rs.m_curZ / zoomFactor;
    T maxZ = rs.m_maxZ / zoomFactor;
    render_compute_bounding_box(rs, minZ, maxZ, LOS, 1u << TObject::LOG2_CHILDREN_EDGE_COUNT);
    if (minZ > maxZ)
      return;
    if (0 == o.get_mipmap(TObject::LOG2_CHILDREN_EDGE_COUNT, 0) /*o.m_mipmap[TObject::LOG2_CHILDREN_EDGE_COUNT][0]*/)
      return;
    
    // for render_skip_box
    BoundingBoxParams<T, N, TColor, uint32_t> bbp;
    bbp.m_rs = &rs;
    bbp.m_adjustedLOS = &LOS;
    bbp.m_children_count_mask = ((1u << TObject::LOG2_CHILDREN_EDGE_COUNT) - 1);
    
    register Vec<uint32_t, N> pos;

    //cout << endl << "screenPoint " << screenPoint << " minZ " << minZ << " maxZ " << maxZ << endl;

    register uint32_t& mipmap = bbp.m_mipmap;
    register T& curZ = bbp.m_curZ;
    
     mipmap = TObject::LOG2_CHILDREN_EDGE_COUNT;
     curZ =  minZ;
    
    register uint32_t index;
    register uint32_t previous_1_index;
    while (	mipmap <= TObject::LOG2_CHILDREN_EDGE_COUNT
      &&	curZ <= maxZ
    )
    {
      pos = LOS.m_org + (curZ + 0.5*min(0.5, maxZ-curZ)) * LOS.m_dir; // the current pos
      //cout << endl << "Start pos " << pos;
      pos.set_min_max(1u << TObject::LOG2_CHILDREN_EDGE_COUNT - 1);
      //cout << pos << endl;
      
      
      index = o.get_child_index(pos, mipmap);

      if (mipmap == 0)
      {
	//cout << (int)(0 != o.m_children[index]) << " = o.m_children[index=" << index << "] " << " curZ "<< curZ << " " << pos << endl; 
	
	// we zoomed up to the very last contained elements
	if (0 != o.m_children[index])
	{
	  // we got it
	  rs.m_color = o.m_children[index];
	  rs.m_bestZ = curZ * zoomFactor;
	  //cout << "Found color at " << pos << " => " << pos * zoomFactor + objectCenterPos - objectEdgeSize/2  << endl;
	  return;
	}
	
	// skip this box, do not zoom out just yet 
	if (!render_skip_box(bbp, pos))
	  return;	  
	
	if (previous_1_index != o.get_child_index(pos, 1))
	{
	  //cout << "mipmap-zoom-out " << mipmap << " ==> " << mipmap+1 << endl;
	  mipmap++; // zoom out 1	  
	}
      }
      else
      {
	if (0 == o.get_mipmap(mipmap, index) /*o.m_mipmap[mipmap][index]*/)
	{
	  // skip this box
	  // and may be go up in the mipmap ?
	  //cout << (int)o.m_mipmap[mipmap][index] << " = o.m_mipmap[mipmap=" << mipmap << "][index=" << index << "] " << " curZ "<< curZ << " " << pos << endl; 
	  //cout << "mipmap-zoom-out " << mipmap;	  
	  // zoom out to the topmost empty box that encloses this one
	  while (mipmap < TObject::LOG2_CHILDREN_EDGE_COUNT)
	  {
	    if (0 == o.get_mipmap(mipmap+1, o.get_child_index(pos, mipmap+1)) /* o.m_mipmap[mipmap+1][o.get_child_index(pos, mipmap+1)] */ )
	      mipmap++;
	    else
	      break;	    
	  }
	  //cout << " => " << mipmap << endl;
	  
	  if (!render_skip_box(bbp, pos))
	    return;
	}
	else
	{
	  //cout << (int)o.m_mipmap[mipmap][index] << " = o.m_mipmap[mipmap=" << mipmap << "][index=" << index << "] " << " curZ "<< curZ << " " << pos << endl; 
  	  //cout << "mipmap-zoom-in " << mipmap;	  

	  // here we should zoom to find the smallest non empty mipmap
	  while (mipmap > 0 && 0 != o.get_mipmap(mipmap, index) /*o.m_mipmap[mipmap][index]*/)
	  {
	    index = o.get_child_index(pos, --mipmap);
	  }
	  //cout << " => " << mipmap << endl;
	  
	  if (mipmap == 0)
	  {
	    if (0 != o.m_children[index])
	    {
	      // we got it
	      rs.m_color = o.m_children[index];
	      rs.m_bestZ = curZ * zoomFactor;
	      //cout << "Found color at " << pos << " => " << (pos * zoomFactor) + objectCenterPos - objectEdgeSize/2  << endl;
	      return;
	    }
	    
	    previous_1_index = o.get_child_index(pos, 1);
	  }	  
	} // if ([mipmap] == 0) ... else
      } // 
    } // while
    /**/
  }
  
}


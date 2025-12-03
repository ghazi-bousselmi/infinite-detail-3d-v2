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

#include "../container/color.h"
#include "../vector/vector.h"


#ifndef __in 
#	define __in
#endif
#ifndef __out
#	define __out
#endif
#ifndef __inout
#	define __inout
#endif
#ifndef ASSERT
#   if CUBOID_VERBOSE
#	define CUBOID_ASSERT(x) { if (!(x)) printf("Warning false == #x\n"); }
#   else
#	define CUBOID_ASSERT(x) 
#   endif
#endif

#ifndef max
#	define max(x, y) ((x) > (y) ? (x) : (y))
#endif
#ifndef min
#	define min(x, y) ((x) < (y) ? (x) : (y))
#endif




/*
 * Abstract camera object, for space dimensions >= 2
 * 
 * the dimension Ds of the screen is 1 less than the dimension D of the space
 * the screen is the place where the image of the space will be projected
 * the screen is a "hyper-rectangle" in the hyper-plane perpendicular to the last dimension
 * the screen has (Ds = D-1) resolutions for each of the first (Ds = D-1) dimensions of space
 * the screen coordinates range from (0, ..., 0) ... (R[0]-1, ..., R[Ds-1]-1) ; in a 3d space, (0, 0) to (res_x-1, res_y-1)
 * 
 * the camera is the representation of the light detector that would be present inside the space
 * the detector is a hyper-surface of dimension Ds
 * the detector surface has indicative Ds edge sizes, for each of the first Ds dimensions of space
 * the center of the camera detector is at the center of the space (for flat detector, perpendicular to the last dimension of space)
 * the last dimension of space is the default line of sight, the direction is from negative to positive coordinates along
 * 
 * to each point on the screen corresponds a point on the surface of the camera detector
 * 
 * on a camera with a flat detector matrix:
 * the screen coordinates correspond to the space coordinates (-edge_size[0]/2, ..., -edge_size[Ds-1]/2) ... (edge_size[0] - edge_size[0]/2, ..., edge_size[Ds-1] - edge_size[Ds-1]/2)
 * 
 * for a camera with a hyper-sphere detector matrix:
 * the surface is a part of hyper-sphere, having a center placed along the last dismension of space, with a coordinate "center_ordinate" (ranging from -infinity to +infinity)
 * the flat detector camera is a limit case of hyper-sphere detector where the "center_ordinate = -infinity"
 * the detector_h_center = DHC = (0, ..., 0, center_ordinate)  (in D-space)
 * the radius of the sphere is 
 * 	Radius  = length( (edge_size[0]/2, ..., edge_size[Ds-1]/2, 0) - DHC )
 * 		= length( (edge_size[0]/2, ..., edge_size[Ds-1]/2, -center_ordinate) )
 * 		= sqrt( (edge_size[0]/2)^2 + ... + (edge_size[Ds-1]/2)^2 + center_ordinate*center_ordinate  )
 * 
 * 
 * 
 * 
 * the eye point (a real point inside space) is at the space coordinate focal_ordinate on the last dimension
 * "focal_ordinate" is a value in the coordinate system of the space
 * note that "focal_ordinate" is a relative number and can be negative
 * 
 * 
 * 
 * 
 * 
 * */

template<typename Element, typename TFloat, int N /* N = SpaceDimentionality, must be >= 2 */, typename LOSElement = Element, bool bNormalizeLOS = true>
class CAbstractCamera
{
public:
  explicit CAbstractCamera() {}
  virtual ~CAbstractCamera() {}
  
  typedef Vec<Element,	N>		TSpace;
  typedef Vec<TFloat,	N>		TSpaceF;
  typedef Vec<LOSElement,N>		TSpaceLOS;
  typedef Vec<Element,	N-1 >		TEdgeSizes;
  typedef Vec<int,N-1 >			TResolution;

  typedef struct
  {
    TSpaceLOS	delta0;	// delta = delta0 + z * A ; where z is the usual z depth, in world units
    TSpaceLOS	A;	
  }
  STLineOfSight_delta;

  typedef struct
  {
    __out TSpace	spacePoint;
    __out TSpaceF	cameraSurfaceNormal;
    __out TSpaceF	cameraSurfaceNormal_NotNull;	// every 0 component of cameraSurfaceNormal replaced with 1.0
    __out TSpaceF	cameraSurfaceVectors[N-1];
    __out TSpaceLOS	lineOfSight;
    __out LOSElement	lineOfSightLength;    
    __out TFloat	solidAngleRadiusA;
    __out TFloat	solidAngleRadiusB; 
    
    __out STLineOfSight_delta	lineOfSight_delta_DConst[N-1][N];		// [m] is the delta with the next adjacent at the resolution dimension [m]
										// [m][l] is the delta with space dimension [l] constant
  } 
  STSpaceSurfacePoint;	// a point on the camera surface corresponding to a screen point
  
  typedef struct
  {
    __in  TResolution		screenPoint;
    __out STSpaceSurfacePoint   spacePoint;  
  } 
  STScreenToSpace;  
  
  virtual bool SetEyeOrdinate	(const Element		&eyeOrdinate);
  virtual bool SetEdgeSizes	(const TEdgeSizes	&edgeSizes);
  virtual bool SetResolution	(const TResolution	&resolution);
  
  inline Element 	       	GetEyeOrdinate() { return m_eyePoint[N - 1];	}
  inline const TSpace&		GetEyePoint()	 { return m_eyePoint;		}
  inline const TEdgeSizes&  	GetEdgeSizes() 	 { return m_edgeSizes;		}
  inline const TResolution& 	GetResolution()  { return m_resolution; 	}
  

  // retrieves the space point on the camera detector surface corresponding to a screen point
  virtual void GetSpacePoint(
    __in const 		TResolution& screenPoint, 
    __out TSpace&	spacePoint, 
    __out TSpaceF&	cameraSurfaceNormal, 
    __out TSpaceF*	cameraSurfaceVectors,
    __out TSpaceLOS&	lineOfSight,
    __out LOSElement&	lineOfSightLength,
    __out TFloat&	solidAngleRadiusA,
    __out TFloat&	solidAngleRadiusB  
    ) = 0;
  
  virtual void GetSpacePoint(__in const TResolution& screenPoint, __out STSpaceSurfacePoint& s);
  virtual void GetSpacePoint(STScreenToSpace& s2s);
  
  inline void GetSpacePoint (__in const TResolution& screenPoint, __out STSpaceSurfacePoint*& s)
  {
    register int nIndex = screenPoint[N-2];
    for (int i = N-3; i >= 0; i--)
      nIndex = nIndex * m_resolution[i] + screenPoint[i];
    s = &m_arrayCameraPoints[nIndex];
  }
  inline void GetSpacePoint(int nIndex, __out STSpaceSurfacePoint*& s)
  {
    s = &m_arrayCameraPoints[nIndex];
  }

protected:
  inline void RecalcCameraPoints()
  {
    unsigned int nPoints = 1;
    for (unsigned int i = 0; i < N-1 ; i++)
      nPoints = nPoints * m_resolution[i];
  
    m_arrayCameraPoints.resize(nPoints);
    
    TResolution r;
    for (unsigned int i = 0; i < nPoints ; i++)
    {
      unsigned int n = i;
      for (unsigned int j = 0; j < N-1 ; j++)
      {
	r[j] = n % m_resolution[j];
	n = n / m_resolution[j];
      }
      
      GetSpacePoint(r, m_arrayCameraPoints[i]);
    }
    
    if (nPoints <= 0)
      return;
    
    return;
    
    TResolution rweigts;
    for (unsigned int i = 0 ; i < N - 1 ; i++)
    {
      r[i] = 0;
      if (i == 0)
	rweigts[i] = 1;
      else
	rweigts[i] = rweigts[i-1] * m_resolution[i-1];
    }

    if (false)
    {
      printf("resolution ( ");
      for (unsigned int i = 0 ; i < N - 1 ; i++)
	printf("%i ", (int)m_resolution[i]);
      printf(")\n");

      printf("weights ( ");
      for (unsigned int i = 0 ; i < N - 1 ; i++)
	printf("%i ", rweigts[i]);
      printf(")\n");
    }
    
    printf("calculating delta line of sight\n");
    for (;;)
    {
      unsigned int index = 0;
      for (unsigned int i = 0 ; i < N - 1 ; i++)
	index += r[i] * rweigts[i];
      STSpaceSurfacePoint* P0 = &m_arrayCameraPoints[index];

      if (false)
      {
	printf("( ");
	for (unsigned int i = 0 ; i < N - 1 ; i++)
	  printf("%i ", r[i]);
	printf(") = %i : SpacePoint( ", index);
	for (unsigned int i = 0 ; i < N ; i++)
	  printf("%7.5f ", (float)P0->spacePoint[i]);
	printf(") LOS( ");
	for (unsigned int i = 0 ; i < N ; i++)
	  printf("%7.5f ", (float)P0->lineOfSight[i]);
	printf(")\n");
      }
      
      
      for (int nScreenDimension = 0 ; nScreenDimension < N-1 ; nScreenDimension++)
      {
	STSpaceSurfacePoint* P1 = 0;
	TSpace deltaSpacePoint;
	if ( r[nScreenDimension] < m_resolution[nScreenDimension] )
	{
	  P1 = &m_arrayCameraPoints[index + rweigts[ nScreenDimension] ];
	  deltaSpacePoint = P1->spacePoint - P0->spacePoint;
	}
	  
	for (int nConstSpaceDimension  = 0 ; nConstSpaceDimension < N    ; nConstSpaceDimension++ )
	{
	  // calculate the delta, when the dimension "nConstSpaceDimension" is constant
	  if (P1 && P0->lineOfSight[ nConstSpaceDimension ] != 0.0f && P1->lineOfSight[ nConstSpaceDimension ] != 0.0f)
	  {
	    // here we can have the same value for dimension [ nConstSpaceDimension ] in both lineOfSight
	    P0->lineOfSight_delta_DConst[ nScreenDimension ][ nConstSpaceDimension ].delta0 = deltaSpacePoint + ((0.0f - deltaSpacePoint[ nConstSpaceDimension ]) / P1->lineOfSight[ nConstSpaceDimension ]) * P1->lineOfSight ;
	    P0->lineOfSight_delta_DConst[ nScreenDimension ][ nConstSpaceDimension ].A      = ((0.0f + P0->lineOfSight[ nConstSpaceDimension ]) / P1->lineOfSight[ nConstSpaceDimension ]) * P1->lineOfSight - P0->lineOfSight ;

	    P0->lineOfSight_delta_DConst[ nScreenDimension ][ nConstSpaceDimension ].delta0[ nConstSpaceDimension ] = 0;
	    P0->lineOfSight_delta_DConst[ nScreenDimension ][ nConstSpaceDimension ].A	   [ nConstSpaceDimension ] = 0;
	  }
	  else
	  {
	    // we could have 0 or 1 common point, or the full 2 lines fulfilling the conditions.
	    // for all intends and purposes, not useful
	    // just set this delta so that it is huge, so that it overflows any optimization
	    for (int i = 0 ; i < N ; i++ )
	    {
	      P0->lineOfSight_delta_DConst[ nScreenDimension ][ nConstSpaceDimension ].delta0[ i ]	= 1e20f;
	      P0->lineOfSight_delta_DConst[ nScreenDimension ][ nConstSpaceDimension ].A[ i ]		= 1e20f;
	    }
	  }
	
	  if (false)
	  {
	    printf("  S:%i ; %i=( ", nScreenDimension, nConstSpaceDimension);
	    for (int i = 0 ; i < N ; i++ )
	      printf("%s%7.5f ", P0->lineOfSight_delta_DConst[ nScreenDimension ][ nConstSpaceDimension ].delta0[ i ] >= 0.0 ? " " : "", (float)P0->lineOfSight_delta_DConst[ nScreenDimension ][ nConstSpaceDimension ].delta0[ i ]);
	    printf(") A( ");
	    for (int i = 0 ; i < N ; i++ )
	      printf("%s%5.3e ", P0->lineOfSight_delta_DConst[ nScreenDimension ][ nConstSpaceDimension ].A     [ i ] >= 0.0 ? " " : "", (float)P0->lineOfSight_delta_DConst[ nScreenDimension ][ nConstSpaceDimension ].A     [ i ]);
	    printf(")\n");
	  }
	  
	}	
      }
      
      // NEXT index
      for (int i = 0 ; i < N-1 ; i++)
      {
	r[i]++;
	if (i < N-2 && r[i] >= m_resolution[i])
	  r[i] = 0;
	else
	  break;
      }
      
      // EXIT condition
      if (r[N-2] >= m_resolution[N-2])
	break;
    }
  }
  
  virtual void OnParametersChange();
  
  
protected:
  TSpace		m_eyePoint;
  TEdgeSizes		m_edgeSizes;
  TResolution		m_resolution;
public:
  std::vector<STSpaceSurfacePoint> m_arrayCameraPoints;
};





template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
bool CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS>::SetEyeOrdinate	(const Element		&eyeOrdinate)
{
  for (int i = 0; i < N ; i++)
    m_eyePoint[i] = 0;
  m_eyePoint[N - 1] = eyeOrdinate;
  
  OnParametersChange();
  return true;
}

template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
bool CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS>::SetEdgeSizes	(const TEdgeSizes	&edgeSizes)
{
  for (int i = 0; i < N - 1 ; i++)
    if (edgeSizes[i] <= 0)
      return false;
  m_edgeSizes = edgeSizes;
  
  OnParametersChange();
  return true;    
}

template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
bool CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS>::SetResolution	(const TResolution	&resolution)
{
  for (int i = 0; i < N - 1 ; i++)
    if (resolution[i] < 4)
      return false;
  m_resolution = resolution;
  
  OnParametersChange();
  return true;
}

template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
void CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS>::GetSpacePoint(__in const TResolution& screenPoint, __out STSpaceSurfacePoint& s)
{
  GetSpacePoint(screenPoint,s.spacePoint,s.cameraSurfaceNormal,s.cameraSurfaceVectors,s.lineOfSight,s.lineOfSightLength,s.solidAngleRadiusA,s.solidAngleRadiusB);
}

template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
void CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS>::GetSpacePoint(STScreenToSpace& s2s)
{
  GetSpacePoint(s2s.screenPoint,
		s2s.spacePoint.spacePoint,s2s.spacePoint.cameraSurfaceNormal,s2s.spacePoint.cameraSurfaceVectors,
		s2s.spacePoint.lineOfSight,s2s.spacePoint.lineOfSightLength,
		s2s.spacePoint.solidAngleRadiusA,
		s2s.spacePoint.solidAngleRadiusB);
}

template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
void CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS>::OnParametersChange()
{
  RecalcCameraPoints();
}







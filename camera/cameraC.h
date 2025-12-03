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

#include "cameraB.h"
#include "math.h"

/////////////////////////////////////////////////////////////////////////////
//
// hybrid flat-Spherical detector camera
// the surface of the detector on each axis is either flat or circular
// the flat camera is the limit of the hybrid where the surface is flat on all axes
// the spherical camera is the limit of the hybrid where the surface is circular on all axes
// 

template<typename Element, typename TFloat, int N /* N = SpaceDimentionality, must be >= 2 */, typename LOSElement = Element, bool bNormalizeLOS = true>
class CCameraC: public CCameraB<Element, TFloat, N, LOSElement, bNormalizeLOS>
{
public:
  explicit CCameraC() {}
  virtual ~CCameraC() {}

  typedef CCameraB<Element, TFloat, N, LOSElement, bNormalizeLOS> TParent;
  typedef Vec<bool, N>  TCircular;

  virtual bool SetCircular		(int nDim, const bool&bCircular){ m_circular[nDim] = bCircular;	OnParametersChange(); return true; }
  virtual bool SetCircular		(const TCircular &circular)	{ m_circular = circular; 	OnParametersChange(); return true; }
  inline  bool GetCircular		(int nDim)			{ return m_circular[nDim]; }
  inline  const TCircular& GetCircular	(const TCircular &circular)	{ return m_circular; }
  
  
  virtual void GetSpacePoint(
    __in const typename TParent::TResolution& 	screenPoint, 
    __out typename TParent::TSpace&		spacePoint, 
    __out typename TParent::TSpaceF& 		cameraSurfaceNormal, 
    __out typename TParent::TSpaceF*		cameraSurfaceVectors,
    __out typename TParent::TSpaceLOS& 		lineOfSight,
    __out LOSElement&				lineOfSightLength,
    __out TFloat&  				solidAngleRadiusA,	// a x + b
    __out TFloat&  				solidAngleRadiusB  
    );
  
protected:
  virtual void OnParametersChange();
  
protected:
  TCircular	m_circular;
};









template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
void CCameraC<Element, TFloat, N, LOSElement, bNormalizeLOS>::GetSpacePoint(
  __in const typename TParent::TResolution& 	screenPoint, 
  __out typename TParent::TSpace&		spacePoint, 
  __out typename TParent::TSpaceF& 		cameraSurfaceNormal, 
  __out typename TParent::TSpaceF*		cameraSurfaceVectors,
  __out typename TParent::TSpaceLOS&		lineOfSight,
  __out LOSElement&				lineOfSightLength,
  __out TFloat&  				solidAngleRadiusA,	// a x + b
  __out TFloat&  				solidAngleRadiusB  
  )
{
  solidAngleRadiusA = 0;
  solidAngleRadiusB = 0;
  memset(&cameraSurfaceNormal[0], 0, N * sizeof(cameraSurfaceNormal[0]));
  memset(&spacePoint[0],          0, N * sizeof(spacePoint[0]));
  cameraSurfaceNormal[N-1] = 1;

#ifndef PI
#	define PI 3.14159265359
#endif
  
  for (int i = 0; i < N-1; i++)
  {
    if (m_circular[i])
    {
      // spherical part
      TFloat A = this->m_apertures[i] * (-1.0 + screenPoint[i] * 2.0 / this->m_resolution[i]);
      cameraSurfaceNormal    *= cos( A );
      cameraSurfaceNormal[i]  = sin( A );
      solidAngleRadiusB = solidAngleRadiusB 
	+   0.25 * (this->m_edgeSizes[i] * 2*this->m_apertures[i] / (2*PI) / (TFloat)this->m_resolution[i])
	         * (this->m_edgeSizes[i] * 2*this->m_apertures[i] / (2*PI) / (TFloat)this->m_resolution[i]);  
    }
    else
    {
      // flat part
      spacePoint[i] = -this->m_edgeSizes[i] / 2 + (this->m_edgeSizes[i] * screenPoint[i]) / this->m_resolution[i];
      cameraSurfaceNormal[i] = 0;
      solidAngleRadiusB = solidAngleRadiusB 
	+ 0.25 * (this->m_edgeSizes[i] * 1.0 / (TFloat)this->m_resolution[i]) 
	*        (this->m_edgeSizes[i] * 1.0 / (TFloat)this->m_resolution[i]);
    }
  }
  spacePoint = spacePoint + (this->m_CenterPoint + cameraSurfaceNormal * this->m_Radius);
   
  lineOfSight = spacePoint - this->m_eyePoint;
  lineOfSightLength = lineOfSight.length();
  
  solidAngleRadiusB = sqrt(solidAngleRadiusB);
  solidAngleRadiusB *= (lineOfSight * cameraSurfaceNormal) / lineOfSightLength; // cos( lineOfSight . surfaceNormal )
  solidAngleRadiusA = solidAngleRadiusB / lineOfSightLength;
  
  if (bNormalizeLOS)
  {
    lineOfSight.normalize();
  }  
}  

template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
void CCameraC<Element, TFloat, N, LOSElement, bNormalizeLOS>::OnParametersChange()
{
  /////////////////////////////////////////////////////////////
  // first calculate the circular axes
  typename TParent::TSpaceF P0;
  memset(&P0[0], 0, N*sizeof(P0[0]));
  int nFirstCircularIndex = N;
  for (int i = 0; i < N-1; i++)
  {
    if (m_circular[i])
    {
      P0[i] = this->m_edgeSizes[i] / 2;
      nFirstCircularIndex = i;
    }
  }
  P0 = P0 - this->m_CenterPoint;
  
  this->m_Radius = P0.length();
  P0.normalize();
  
  for (int i = N-2; i >= 0; i--)
  {
    if (m_circular[i])
    {
      this->m_apertures[i] = asin( P0[i] );
      if (i > nFirstCircularIndex)
      {
	P0[i] = 0;
	P0 = P0 / cos( this->m_apertures[i] );
      }
    }
  }
  
  if (nFirstCircularIndex < N)
  {
    if (this->m_apertures[nFirstCircularIndex] < 0.0)
      this->m_apertures[nFirstCircularIndex] += 2 * asin( 1.0 ); // += PI
  }
  
  // we curcicuited CCameraB::OnParametersChange();
  CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS>::OnParametersChange(); 
}


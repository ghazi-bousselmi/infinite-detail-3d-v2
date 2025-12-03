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

#include "camera.h"
#include "math.h"

/////////////////////////////////////////////////////////////////////////////
//
// Spherical detector camera
//
//

template<typename Element, typename TFloat, int N /* N = SpaceDimentionality, must be >= 2 */, typename LOSElement = Element, bool bNormalizeLOS = true>
class CCameraB: public CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS>
{
public:
  explicit CCameraB() { m_Radius = 1.0; }
  virtual ~CCameraB() {}

  typedef CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS> TParent;

  virtual bool SetEyeOrdinate	(const Element		&eyeOrdinate);
  virtual bool SetCenterOrdinate(const Element		&centerOrdinate);
  
  inline Element 				GetCenterOrdinate()	{ return m_CenterPoint[N - 1]; }
  inline const typename TParent::TSpace&	GetCenterPoint()	{ return m_CenterPoint;        }

  
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
  typename TParent::TSpace	m_CenterPoint;
  Vec<TFloat, N>		m_apertures;	// angles, half aperture angles
  double m_Radius;
};





template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
bool CCameraB<Element, TFloat, N, LOSElement, bNormalizeLOS>::SetEyeOrdinate	(const Element		&eyeOrdinate)
{
  //if (eyeOrdinate >= 0)
  //  return false;
  return TParent::SetEyeOrdinate(eyeOrdinate);
}

template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
bool CCameraB<Element, TFloat, N, LOSElement, bNormalizeLOS>::SetCenterOrdinate(const Element		&centerOrdinate)
{
  for (int i = 0; i < N - 1; i++)
    m_CenterPoint[i] = 0;
  m_CenterPoint[N - 1] = centerOrdinate;
  
  OnParametersChange();
  return true;
}

template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
void CCameraB<Element, TFloat, N, LOSElement, bNormalizeLOS>::GetSpacePoint(
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
  cameraSurfaceNormal[N-1] = 1;

#ifndef PI
#	define PI 3.14159265359
#endif
  
  for (int i = 0; i < N-1; i++)
  {
    TFloat A = m_apertures[i] * (-1.0 + screenPoint[i] * 2.0 / this->m_resolution[i]);
    cameraSurfaceNormal    *= cos( A );
    cameraSurfaceNormal[i]  = sin( A );
    solidAngleRadiusB = solidAngleRadiusB 
      + 0.25 * (this->m_edgeSizes[i] * 2*m_apertures[i] / (2*PI) / (TFloat)this->m_resolution[i])
	     * (this->m_edgeSizes[i] * 2*m_apertures[i] / (2*PI) / (TFloat)this->m_resolution[i]);
  }
  spacePoint = m_CenterPoint + cameraSurfaceNormal * m_Radius;
   
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
void CCameraB<Element, TFloat, N, LOSElement, bNormalizeLOS>::OnParametersChange()
{
  typename TParent::TSpaceF P0;
  P0[N-1] = 0;
  for (int i = 0; i < N-1; i++)
    P0[i] = this->m_edgeSizes[i] / 2;
  
  P0 = P0 - m_CenterPoint;
  
  m_Radius = P0.length();
  P0.normalize();
  
  for (int i = N-2; i >= 0; i--)
  {
    m_apertures[i] = asin( P0[i] );
    if (i > 0)
    {
      P0[i] = 0;
      P0 = P0 / cos( m_apertures[i] );
    }
  }
  
  if (m_apertures[0] < 0.0)
    m_apertures[0] += 2 * asin( 1.0 ); // += PI
  
  TParent::OnParametersChange();  
}


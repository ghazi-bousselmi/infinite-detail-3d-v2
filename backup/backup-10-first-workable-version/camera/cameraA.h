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
#include <math.h>

/////////////////////////////////////////////////////////////////////////////
//
// flat detector camera
//
//

template<typename Element, typename TFloat, int N /* N = SpaceDimentionality, must be >= 2 */, typename LOSElement = Element, bool bNormalizeLOS = true>
class CCameraA: public CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS>
{
public:
  explicit CCameraA() {}
  virtual ~CCameraA() {}
  
  typedef CAbstractCamera<Element, TFloat, N, LOSElement, bNormalizeLOS> TParent;

  virtual bool SetEyeOrdinate	(const Element		&eyeOrdinate);
  
  virtual void GetSpacePoint(
    __in const typename TParent::TResolution& 	screenPoint, 
    __out typename TParent::TSpace&		spacePoint, 
    __out typename TParent::TSpaceF&		cameraSurfaceNormal, 
    __out typename TParent::TSpaceF*		cameraSurfaceVectors,
    __out typename TParent::TSpaceLOS&		lineOfSight,
    __out LOSElement&				lineOfSightLength,
    __out TFloat&  				solidAngleRadiusA,	// a x + b
    __out TFloat&  				solidAngleRadiusB  
    );
};



template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
bool CCameraA<Element, TFloat, N, LOSElement, bNormalizeLOS>::SetEyeOrdinate	(const Element		&eyeOrdinate)
{
  if (eyeOrdinate >= 0)
    return false;
  return TParent::SetEyeOrdinate(eyeOrdinate);
}

template<typename Element, typename TFloat, int N, typename LOSElement, bool bNormalizeLOS>
void CCameraA<Element, TFloat, N, LOSElement, bNormalizeLOS>::GetSpacePoint(
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
  for (int i = 0; i < N-1; i++)
  {
    spacePoint[i] = -this->m_edgeSizes[i] / 2 + (this->m_edgeSizes[i] * screenPoint[i]) / this->m_resolution[i];
    cameraSurfaceNormal[i] = 0;
    solidAngleRadiusB = solidAngleRadiusB 
      + 0.25 * (this->m_edgeSizes[i] * 1.0 / (TFloat)this->m_resolution[i]) 
             * (this->m_edgeSizes[i] * 1.0 / (TFloat)this->m_resolution[i]);

    memset(& ((cameraSurfaceVectors[i])[0]), 0, N * sizeof( (cameraSurfaceVectors[i])[0] )    );
    (cameraSurfaceVectors[i])[i] = 1;
  }
  spacePoint[N-1] = 0;
  cameraSurfaceNormal[N-1] = 1;
  
  //printf("Camera (%i %i) (%1.6f %1.6f %1.6f)\n", screenPoint[0], screenPoint[1], spacePoint[0], spacePoint[1], spacePoint[2]);
  
  lineOfSight = spacePoint - this->m_eyePoint;
  lineOfSightLength = lineOfSight.length();
  
  solidAngleRadiusB = sqrt(solidAngleRadiusB);
  //solidAngleRadiusB *= (lineOfSight * cameraSurfaceNormal) / lineOfSightLength; // cos( lineOfSight . surfaceNormal(=Z) )
  solidAngleRadiusA = solidAngleRadiusB / lineOfSightLength;
  
  //printf("A %f, B %f\n", solidAngleRadiusA, solidAngleRadiusB);
  
  if (bNormalizeLOS)
  {
    lineOfSight.normalize();
  }  
}  

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

#include <chrono>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <sys/time.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <ctime>
#include <ratio>
#include <chrono>
#include <mutex>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <atomic>



#include "vector/vector.h"
#include "line/line.h"
#include "container/container.h"
#include "camera/camera.h"
#include "camera/cameraA.h"
#include "camera/cameraB.h"
#include "camera/cameraC.h"

#include "bitmap/bitmap.h"

using namespace std;
typedef double TFloat;
// for the camera
typedef int32_t BaseType;
typedef double LosType;
const bool bNormalizeLOS = true;

const float granularity = 1;

#define MEMORY_MAX 5*1024*1024

template <class Rep, class Period = std::ratio<1> >
class CDurationAccu
{
public:
  CDurationAccu()   { m_count = 0; }
  inline void start()
  {
    //m_mutex.lock();
    m_starts.push_back(std::chrono::high_resolution_clock::now());
    //m_mutex.unlock();
  }
  inline void stop()
  {
    //m_mutex.lock();
    if (m_starts.begin() != m_starts.end())
    {
      m_count ++;
      m_accu += std::chrono::duration_cast<std::chrono::duration<Rep, Period> > ( std::chrono::high_resolution_clock::now() - m_starts.front() );
      m_starts.pop_front();
    }
    //m_mutex.unlock();
  }
  inline int getRec()
  {
    return m_starts.size();
  }
  inline void stop_all()
  {
    //m_mutex.lock();
    while (m_starts.begin() != m_starts.end())
    {
      m_count ++;
      m_accu += std::chrono::duration_cast<std::chrono::duration<Rep, Period> > ( std::chrono::high_resolution_clock::now() - m_starts.front() );
      m_starts.pop_front();
    }
    //m_mutex.unlock();
  }
  inline std::chrono::duration<Rep, Period> get()
  {
    return m_accu;
  }
  inline int getCount()
  {
    return m_count;
  }
  inline void clear()
  {
    m_count = 0;
    m_accu = std::chrono::duration<Rep, Period>();
    m_starts.clear();
  }
  
protected:
  int m_count;
  std::recursive_mutex	m_mutex;
  std::chrono::duration<Rep, Period>	m_accu;
  std::list<std::chrono::time_point<std::chrono::high_resolution_clock, std::chrono::duration<Rep, Period> > >	m_starts;
};

template<typename TContainer, typename TColor, typename T1>
inline void CreateCube(TContainer* container, double containerEdgeSize,
		       const Vec<T1, 3>& center, double Radius0, 
		       const TColor& c1, const TColor& c2, const TColor& c3, const TColor& c4, const TColor& c5, const TColor& c6 , 
		       int nInterlace)
{
  Vec<T1, 3> coords;
  int nc = 0;
  for (double y = -Radius0; y <= Radius0; y += nInterlace )
  {

    for (double x = -Radius0; x <= Radius0; x += nInterlace , nc += 6)
    {
      coords[0] = center[0] + x;
      coords[1] = center[1] + y;
      coords[2] = center[2] - Radius0;
      container->add_child(c1, coords, Vec<T1, 3>{}, 1, containerEdgeSize);

      coords[0] = center[0] + x;
      coords[1] = center[1] + y;
      coords[2] = center[2] + Radius0;
      container->add_child(c2, coords, Vec<T1, 3>{}, 1, containerEdgeSize);

      coords[0] = center[0] + x;
      coords[1] = center[1] - Radius0;
      coords[2] = center[2] + y;
      container->add_child(c3, coords, Vec<T1, 3>{}, 1, containerEdgeSize);

      coords[0] = center[0] + x;
      coords[1] = center[1] + Radius0;
      coords[2] = center[2] + y;
      container->add_child(c4, coords, Vec<T1, 3>{}, 1, containerEdgeSize);

      
      coords[0] = center[0] - Radius0;
      coords[1] = center[1] + x;
      coords[2] = center[2] + y;
      container->add_child(c5, coords, Vec<T1, 3>{}, 1, containerEdgeSize);

      coords[0] = center[0] + Radius0;
      coords[1] = center[1] + x;
      coords[2] = center[2] + y;
      container->add_child(c6, coords, Vec<T1, 3>{}, 1, containerEdgeSize);
    }   
    
    struct rusage _ru;
    getrusage(RUSAGE_SELF, &_ru);
    
    if ( _ru.ru_maxrss >= MEMORY_MAX )
    {
      printf("Current memory usage = %i MB, max is %i MB, stopping\n",  _ru.ru_maxrss, MEMORY_MAX);
      break;
    }
  }
  
  printf("Cube : points = %i, granularity=%f, Radius=%f\n", nc, granularity);  
}

template<typename TContainer, typename TColor, typename T1>
inline void CreateCircle(TContainer* container, double containerEdgeSize,  
			 const Vec<T1, 3>& center, 
			 double RadiusX, double RadiusY,
			 const Vec<T1, 3>& baseX,
			 const Vec<T1, 3>& baseY,
			 const TColor& c1,
			 const TColor& c2,
			 int nInterlace,
			 int * pnc = 0,
			 bool bVerbose = false
			)
{
  Vec<T1, 3> coords;
  int nc = 0;
  T1 R = max(RadiusX, RadiusY);
  double dA = R > granularity ? granularity / R : 3.1415 / 5;
  TColor rgba;
  for (double a = 0; a <= (2*3.1416); a += dA, nc++)
  {
    double s = sin(a);
    double c = cos(a);
    coords = center + (RadiusX * c) * baseX + (RadiusY * s) * baseY;
    
    rgba.SetWeighted(c1, s*s,      c2, c*c);
    
    container->add_child(rgba, coords, Vec<T1, 3>{}, 1, containerEdgeSize);
  }   

  if (bVerbose)
    printf("circle points = %i, granularity=%f\n", nc, granularity);
  if (pnc)
    *pnc = nc;
}

template<typename TContainer, typename TColor, typename T1>
inline void CreateSphere(TContainer* container, double containerEdgeSize,  
			 const Vec<T1, 3>& center, double Radius0, const TColor& c1, const TColor& c2, int nInterlace)
{
  Vec<T1, 3> coords;
  int nc = 0;
  
  double Radius2 = Radius0 * Radius0;
  double in2 = nInterlace * nInterlace;
  double stmp, ctmp, s, c;
  TColor rgba;

  for (double a1 = 0; a1 <= 3.1416; )
  {
    c  = cos(a1);
    s  = sin(a1);
    double z  = Radius0 * c;
    double R2 = Radius0 * s;
    double dA = R2 > granularity ? granularity / R2 : 3.1415 / 5;
    
    for (double a = 0; a <= (2*3.1416); a += dA, nc++)
    {
      double s__ = sin(a);
      double c__ = cos(a);
      coords[0] = center[0] + R2 * c__;
      coords[1] = center[1] + R2 * s__;
      coords[2] = center[2] + z;
      
      rgba.SetWeighted(c1, s__ * s__,		c2, c__ * c__);
      
      container->add_child(rgba, coords, Vec<T1, 3>{}, 1, containerEdgeSize);
    }   
    
    // now move the angle a1
    double dA1 = 3.1416 / 5;
    do
      dA1 *= 0.9;
    while(Radius2 * (c - (ctmp=cos(a1+dA1)) )*(c - ctmp) + (s - (stmp=sin(a1+dA1)) ) *(s - stmp) > in2);
    a1 += dA1;
  }
  
  printf("Sphere points = %i, granularity=%f\n", nc, granularity);
}

template<typename TContainer, typename TColor, typename T1>
inline void CreateSphere(TContainer* container, double containerEdgeSize,  
			 const Vec<T1, 3>& center, 
			 double RadiusX, double RadiusY, double RadiusZ, 
			 const Vec<T1, 3>& baseX,
			 const Vec<T1, 3>& baseY,
			 const Vec<T1, 3>& baseZ,
			 const TColor& c1,
			 const TColor& c2,
			 const TColor& c3,
			 int nInterlace)
{
  Vec<T1, 3> coords;
  int nc = 0, nnnn;  
  
  if (RadiusZ >= RadiusX && RadiusZ >= RadiusY && RadiusX >= RadiusY)
  {
    double Rx2 = RadiusX * RadiusX;
    double Rz2 = RadiusZ * RadiusZ;
    double in2 = nInterlace * nInterlace;
    double stmp, ctmp, s, c;
    
    for (double a1 = 0; a1 <= 3.1416; )
    {
      c = cos(a1);
      s = sin(a1);

      CreateCircle(container, containerEdgeSize, center + c * RadiusZ * baseZ, s * RadiusX, s * RadiusY, baseX, baseY, c1, c2, nInterlace, &nnnn);
      nc += nnnn;
      
      // now move the angle a1
      double dA1 = 3.1416 / 5;
      do
	dA1 *= 0.9;
      while( Rz2 * (c - (ctmp=cos(a1+dA1)) )*(c - ctmp) + Rx2 * (s - (stmp=sin(a1+dA1)) ) * (s - stmp) > in2);
      a1 += dA1;
    }
  }
  else if (RadiusY >= RadiusX && RadiusY >= RadiusZ)
  {
    if (RadiusX >= RadiusZ)
      CreateSphere(container, containerEdgeSize, center, RadiusX, RadiusZ, RadiusY, baseX, baseZ, baseY, c1, c3, c2, nInterlace);
    else
      CreateSphere(container, containerEdgeSize, center, RadiusZ, RadiusX, RadiusY, baseZ, baseX, baseY, c3, c1, c2, nInterlace);
    return;
  }
  else // if (RadiusX is the max)
  {
    if (RadiusY >= RadiusZ)
      CreateSphere(container, containerEdgeSize, center, RadiusY, RadiusZ, RadiusX, baseY, baseZ, baseX, c2, c3, c1, nInterlace);
    else
      CreateSphere(container, containerEdgeSize, center, RadiusZ, RadiusY, RadiusX, baseZ, baseY, baseX, c3, c2, c1, nInterlace);
    return;
  }  
  
  printf("Sphere points = %i, granularity=%f\n", nc, granularity);
}


template<typename TContainer, typename TColor, typename T1>
inline void CreateTaurus(TContainer* container, double containerEdgeSize,  
			 const Vec<T1, 3>& center, 
			 double RadiusX, double RadiusY, double RadiusZ,
			 const Vec<T1, 3>& baseX,
			 const Vec<T1, 3>& baseY,
			 const TColor& c1,
			 const TColor& c2,
			 int nInterlace)
{
  Vec<T1, 3> coords;
  int nc = 0;
  double R = min(RadiusX, RadiusY);
  double dA = R > granularity ? granularity / R : 3.1415 / 5;
  //Vec<T1, 3> baseZ = baseX ^ baseY;
  Vec<T1, 3> baseZ;
  // cross product
  baseZ[0] = baseX[1] * baseY[2] - baseX[2] * baseY[1];
  baseZ[1] = baseX[2] * baseY[0] - baseX[0] * baseY[2];
  baseZ[2] = baseX[0] * baseY[1] - baseX[1] * baseY[0];
  
  baseZ.normalize();
  
  for (double a = 0; a <= (2*3.1416); a += dA, nc++)
  {
    double s = sin(a);
    double c = cos(a);
    Vec<T1, 3> tmp = (RadiusX * c) * baseX + (RadiusY * s) * baseY;    
    coords = center + tmp;
    
    tmp.normalize();
    
    int nnnn;
    CreateCircle(container, containerEdgeSize, coords, RadiusZ, RadiusZ, tmp, baseZ, c1, c2, nInterlace, &nnnn, false);
    nc += nnnn;
  }   

  printf("taurus points = %i, granularity=%f\n", nc, granularity);
}

template<typename TContainer, typename TColor, typename T1>
inline void CreateLine(TContainer* container, double containerEdgeSize,  
				  const Vec<T1, 3>& center, 
				  double Radius,
				  const Vec<T1, 3>& vector,
				  const TColor& c1,
				  const TColor& c2,
				  int nInterlace,
				  int * pnc = 0,
				  bool bVerbose = false
 				)
{
  Vec<T1, 3> coords;
  int nc = 0;
  TColor rgba;
  
  for (double x = -Radius; x <= Radius ; x += nInterlace, nc++)
  {
    coords = center + x * vector;
    double xx = (x + Radius) / Radius / 2;
    rgba.SetWeighted(
      c1,	xx * xx	,
      c2,	1.0 - xx * xx
    );
    
    container->add_child(rgba, coords, Vec<T1, 3>{}, 1, containerEdgeSize);
  }
  if (bVerbose)
    printf("line points = %i, granularity=%f\n", nc, granularity);
  if (pnc)
    *pnc = nc;
}

template<typename TContainer, typename TColor, typename T1>
inline void CreateParallelogramme(TContainer* container, double containerEdgeSize,  
				  const Vec<T1, 3>& center, 
				  double RadiusX, double RadiusY,
				  const Vec<T1, 3>& baseX,
				  const Vec<T1, 3>& baseY,
				  const TColor& c1,
				  const TColor& c2,
				  bool bFilled,
				  int nInterlace,
				  int *pnc = 0,
				  bool bVerbose = false
				  )
{
  Vec<T1, 3> coords;
  int nc = 0;  
  int nnnn;
  
  if (!bFilled)
  {  
    CreateLine(container, containerEdgeSize, center + RadiusX * baseX, RadiusY, baseY, c1, c2, nInterlace, &nnnn);	nc += nnnn;
    CreateLine(container, containerEdgeSize, center - RadiusX * baseX, RadiusY, baseY, c1, c2, nInterlace, &nnnn);	nc += nnnn;

    CreateLine(container, containerEdgeSize, center + RadiusY * baseY, RadiusX, baseX, c2, c1, nInterlace, &nnnn);	nc += nnnn;
    CreateLine(container, containerEdgeSize, center - RadiusY * baseY, RadiusX, baseX, c2, c1, nInterlace, &nnnn);	nc += nnnn;
  }
  else
  {
    if (RadiusX < RadiusY)
    {
      for (double x = -RadiusX ; x <= RadiusX; x += nInterlace)
      {
	CreateLine(container, containerEdgeSize, center + x * baseX, RadiusY, baseY, c1, c2, nInterlace, &nnnn, false);	
	nc += nnnn;
      }
    }
    else
    {
      for (double y = -RadiusY ; y <= RadiusY; y += nInterlace)
      {
	CreateLine(container, containerEdgeSize, center + y * baseY, RadiusX, baseX, c2, c1, nInterlace, &nnnn, false);	
	nc += nnnn;
      }
    }    
  }
  
  if (bVerbose)
    printf("parallelogramme%s points = %i, granularity=%f\n", bFilled ? "-filled" : "", nc, granularity);
  if (pnc)
    *pnc = nc;
}

template<typename TContainer, typename TColor, typename T1>
inline void CreateParallipipede(TContainer* container, double containerEdgeSize,  
				const Vec<T1, 3>& center, 
				double RadiusX, double RadiusY, double RadiusZ,
				const Vec<T1, 3>& baseX,
				const Vec<T1, 3>& baseY,
				const Vec<T1, 3>& baseZ,
				const TColor& c1,
				const TColor& c2,
				const TColor& c3,
				bool bFilled,
				int nInterlace,
				int *pnc = 0,
				bool bVerbose = false
				)
{
  Vec<T1, 3> coords;
  int nc = 0;  
  int nnnn;
  
  CreateParallelogramme(container, containerEdgeSize, center - RadiusX * baseX, RadiusY, RadiusZ, baseY, baseZ, c2, c3, bFilled, nInterlace, &nnnn); nc += nnnn;
  CreateParallelogramme(container, containerEdgeSize, center + RadiusX * baseX, RadiusY, RadiusZ, baseY, baseZ, c2, c3, bFilled, nInterlace, &nnnn); nc += nnnn;
  
  CreateParallelogramme(container, containerEdgeSize, center - RadiusY * baseY, RadiusX, RadiusZ, baseX, baseZ, c1, c3, bFilled, nInterlace, &nnnn); nc += nnnn;
  CreateParallelogramme(container, containerEdgeSize, center + RadiusY * baseY, RadiusX, RadiusZ, baseX, baseZ, c1, c3, bFilled, nInterlace, &nnnn); nc += nnnn;

  CreateParallelogramme(container, containerEdgeSize, center - RadiusZ * baseZ, RadiusX, RadiusY, baseX, baseY, c1, c2, bFilled, nInterlace, &nnnn); nc += nnnn;
  CreateParallelogramme(container, containerEdgeSize, center + RadiusZ * baseZ, RadiusX, RadiusY, baseX, baseY, c1, c2, bFilled, nInterlace, &nnnn); nc += nnnn;
  
  if (bVerbose)
    printf("parallepipede%s points = %i, granularity=%f\n", bFilled ? "-filled" : "", nc, granularity);
  if (pnc)
    *pnc = nc;
}




// int main(int argc, char** argv)
// {
//   Vec<TFloat, 3> v1 {1, 2, 3}, v2 {4, 5, 6};
//   
//   cout << v1 << endl;
//   cout << v2 << endl;
//   cout << v1*v2 << endl;
//   cout << v1%v2 << endl;
//   cout << (v1^v2) << endl;
//   
//   Line<TFloat, 3> l;
//   l.m_org = {0, 0, 0};
//   l.m_dir = {1, 1, 0};
//   l.m_dir.normalize();
// 
//   PlaneDirs<TFloat, 3> d;
//   d.m_dirs[0] = {0, 1, 0};
//   d.m_dirs[1] = {0, 0, 1};
//   d.calc_normal();
//   
//   Plane<TFloat, 3> p(Vec<TFloat, 3>{4, 0, 0}, &d);
//   
//   auto x = p & l;
//   auto pnt = l.m_org + x * l.m_dir;
//   cout << pnt << endl;
//   
//   RenderStruct<TFloat, 3, Color> rs;
//   rs.m_screenPoint = {0, 0, 0};
//   rs.m_LOS = {1, 0, 0};
//   rs.m_solid_angle_radius.m_b = 1;
//   rs.m_solid_angle_radius.m_a = 1 / 10.0; // the eye point is at 10 units behing the screenPoint
//   
//   Color c;
//   c.u32 = 0xFF0000FF;
// 
//   auto o = new Container<TFloat, 3, Color/*TColor*/, Color/*TChild*/>(4);
//   o->add_child(c, Vec<TFloat, 3>{0.0, 0.0, 0.0}, Vec<TFloat, 3>{0, 0, 0}, 4, 1.0 * (1 << 4));
//   
//   constexpr uint32_t Log2 = 10;
//   
//   auto container = new Container<TFloat, 3, Color/*TColor*/, Bridge<TFloat, 3, Color>/*TChild*/>(Log2);
//   //for (int y = 0; y < (1u << Log2); y++)
//   //for (int z = 0; z < (1u << Log2); z++)
//   //  container->add_child(o, Vec<TFloat, 3>{0.0, y*10.0 - (1u << Log2)/2*10, z*10.0 - (1u << Log2)/2*10}, Vec<TFloat, 3>{0, 0, 0}, 50, 10.0 * (1u << Log2));
//   for (int y = 0; y <= 0; y++)
//   for (int z = 0; z <= 0; z++)
//     container->add_child(o, Vec<TFloat, 3>{0.0, y*10.0, z*10.0}, Vec<TFloat, 3>{0, 0, 0}, 50, 10.0 * (1u << Log2));
//   
//   cout << "rendering started" << endl;
//  
//   using namespace std::chrono;
//   high_resolution_clock::time_point t1 = high_resolution_clock::now();
// 
//   for (int i = 0; i < 1/*5*1000*1000 /**/; i++)
//     render(rs, *container, Vec<TFloat, 3>{500, 0, 0}, 1.0 * (1u << Log2));
//   
//   cout << "rendering done" << endl;
//   
//   high_resolution_clock::time_point t2 = high_resolution_clock::now();
// 
//   duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
// 
//   cout << "It took me " << time_span.count() << " seconds." << endl;
//   
//   cout << "found color " << rs.m_color.u32 << endl;
//   
//   int i;
//   cin >> i;
//   
//   return 0;
// }


template<typename T>
void ____exit(const T t)
{
  cout << t << endl;
  exit(-1);
}


int main(int argc, char ** argv)
{
  struct rusage _ru;    
  
  if (argc > 1 || !argv[0])
  {
    printf("No arguments accepted (for now). Exiting.\n");
    exit(-1);
  }
  
  const BaseType tcameraSizeFactor = 2;
  
  auto camera = new CCameraA<BaseType, TFloat, 3, LosType, bNormalizeLOS>();
    
  CAbstractCamera<BaseType, TFloat, 3, LosType, bNormalizeLOS>::TEdgeSizes sizes;
  CAbstractCamera<BaseType, TFloat, 3, LosType, bNormalizeLOS>::TResolution resolution;

  if (!camera->SetEyeOrdinate(-800* tcameraSizeFactor ))
    ____exit("error in cam.SetEyeOrdinate()"); 
  
  resolution[0] = 1200;
  resolution[1] = 600;
  sizes[0] = resolution[0] * tcameraSizeFactor;
  sizes[1] = resolution[1] * tcameraSizeFactor;
  if (!camera->SetEdgeSizes(sizes))
    ____exit("error in cam.SetEdgeSizes()");
  

  if (!camera->SetResolution(resolution))
    ____exit("error in cam.SetResolution()");
  
  printf("Debug 2\n"); fflush(stdout);

  auto bitmap = new Color<uint8_t>[camera->GetResolution()[0] * camera->GetResolution()[1]];
  
  
  //int then, now;
  CDurationAccu<double> dur;
 
  constexpr uint32_t Log2 = 10;
  auto container = new Container<TFloat, 3, Color<uint8_t>/*TColor*/, Bridge<TFloat, 3, Color<uint8_t>>/*TChild*/>(Log2);
  getrusage(RUSAGE_SELF, &_ru); cout << "Current memory usage = " <<  _ru.ru_maxrss / 1024.0 << " MB" <<  endl;
  
  {	
    int k; 
    printf("paused, press a key"); scanf("%c", &k);  	    
  }
  
  //then = time(0);
  dur.clear(); dur.start();

  TFloat object_edgeSizeOverall = 1024;
  TFloat container_EdgeSizeOverall = 1024*1024;
  TFloat world_edgeSizeOverall = container_EdgeSizeOverall;
  TFloat Radius0 = object_edgeSizeOverall / 2;
      
  Vec<TFloat, 3> zero_center{};
  
  
  Color<uint8_t> rgba; 
  Color<uint8_t> rgba__;
  Color<uint8_t> rgba_green;
  Color<uint8_t> rgba_red;  
  Color<uint8_t> rgba_blue; 

  Color<uint8_t> rgba_rb;   
  Color<uint8_t> rgba_rg;   
  Color<uint8_t> rgba_bg;   

  rgba	    .u32 = 0xf0F0f000;
  rgba__    .u32 = 0xf0F0f000;
  rgba_green.u32 = 0xf000F000;
  rgba_red  .u32 = 0xf0F00000;
  rgba_blue .u32 = 0xf00000F0;

  rgba_rb  .u32	= 0xf0F000F0;
  rgba_rg  .u32	= 0xf0F0F000;
  rgba_bg  .u32	= 0xf000F0F0;

  printf("original color x0%08X, Radius0=%1.2f\n", rgba, Radius0);
  
  
  /**/
  auto cube 		= new Container<TFloat, 3, Color<uint8_t>/*TColor*/, Color<uint8_t>/*TChild*/>(10);
  auto sphere 		= new Container<TFloat, 3, Color<uint8_t>/*TColor*/, Color<uint8_t>/*TChild*/>(10);
  auto taurus 		= new Container<TFloat, 3, Color<uint8_t>/*TColor*/, Color<uint8_t>/*TChild*/>(10);
  auto parallelipipede 	= new Container<TFloat, 3, Color<uint8_t>/*TColor*/, Color<uint8_t>/*TChild*/>(10);
  auto sphere_2 	= new Container<TFloat, 3, Color<uint8_t>/*TColor*/, Color<uint8_t>/*TChild*/>(10);
  printf("Just created empty objects : ");
  getrusage(RUSAGE_SELF, &_ru); cout << "Current memory usage = " <<  _ru.ru_maxrss / 1024.0 << " MB" <<  endl;
  
  {
    Vec<TFloat, 3> baseX({ 1.0f, 0.0f, 0.0f });
    Vec<TFloat, 3> baseY_Bias({ 0.0, sqrt(2.0)/2.0, sqrt(2.0)/2.0 });
    Vec<TFloat, 3> baseY({ 0.0f, 1.0f, 0.0f });
    Vec<TFloat, 3> baseZ({ 0.0f, 0.0f, 1.0f});

    CreateParallipipede(parallelipipede, object_edgeSizeOverall, zero_center, Radius0 * 0.5, Radius0 * 0.5, Radius0 * 0.5, baseX, baseY, baseZ, rgba_blue, rgba_green, rgba_red, true, granularity, 0, true);
    getrusage(RUSAGE_SELF, &_ru); cout << "color" << color_of(parallelipipede	) << " ; Current memory usage = " <<  _ru.ru_maxrss / 1024.0 << " MB" <<  endl;
    
    CreateSphere(sphere, object_edgeSizeOverall, zero_center, Radius0 * 0.5, Radius0 * 0.5, Radius0 * 0.5, baseX, baseY, baseZ, rgba_blue, rgba_green, rgba_green, granularity);
    getrusage(RUSAGE_SELF, &_ru); cout << "color" << color_of(sphere		) << " ; Current memory usage = " <<  _ru.ru_maxrss / 1024.0 << " MB" <<  endl;
    
    CreateTaurus(taurus, object_edgeSizeOverall, zero_center, Radius0 * 0.8, Radius0 * 0.8, Radius0 / 6, baseX, baseY, rgba_green, rgba_red, granularity);
    getrusage(RUSAGE_SELF, &_ru); cout << "color" << color_of(taurus		) << " ; Current memory usage = " <<  _ru.ru_maxrss / 1024.0 << " MB" <<  endl;
    
    //CreateCube(cube, zero_center, Radius0, granularity);
    CreateParallipipede(cube, object_edgeSizeOverall, zero_center, Radius0*0.8, Radius0  * 0.5, Radius0  * 0.5, baseX, baseY_Bias, baseZ, rgba_red, rgba_green, rgba_blue, true, granularity, 0, true);
    getrusage(RUSAGE_SELF, &_ru); cout << "color" << color_of(cube			) << " ; Current memory usage = " <<  _ru.ru_maxrss / 1024.0 << " MB" <<  endl;
    
    CreateSphere(sphere_2, object_edgeSizeOverall, zero_center, Radius0 * 0.33, Radius0 * 0.75, Radius0 * 0.1, baseX, baseY_Bias, baseZ, rgba_blue, rgba_blue, rgba_red, granularity);
    getrusage(RUSAGE_SELF, &_ru); cout << "color" << color_of(sphere_2		) << " ; Current memory usage = " <<  _ru.ru_maxrss / 1024.0 << " MB" <<  endl;

    
    {	
      int k; 
      printf("paused, press a key"); scanf("%c", &k);  	    
    }    
  }
  /**/
  
  
  cube 			->recalc_mipmap_color();
  sphere 		->recalc_mipmap_color();
  taurus 		->recalc_mipmap_color();
  parallelipipede 	->recalc_mipmap_color();
  sphere_2 		->recalc_mipmap_color();
  cout << "color" << color_of(cube	) <<  endl;
  cout << "color" << color_of(sphere	) <<  endl;
  cout << "color" << color_of(taurus	) <<  endl;
  cout << "color" << color_of(parallelipipede	) <<  endl;
  cout << "color" << color_of(sphere_2	) <<  endl;
  
  
  dur.stop();
  printf("DT = %lf\n", dur.get().count());

  //then = time(0);
  dur.clear(); dur.start();
  int nO = 0;
  int x=0, y=0, z=0;
  for (z = -10; z <= 10; z++)
  {
    printf("z=%i, ", z);
    for (y = -10; y <= 10; y++)
    {
      for (x = -10; x <= 10; x++, nO++)
      {
	const double dObjCoordsFactor = 4;
	Vec<TFloat, 3> C;
	C[0] = x * object_edgeSizeOverall * 2;
	C[1] = y * object_edgeSizeOverall * 2;
	C[2] = z * object_edgeSizeOverall * 2;
	
	auto p = cube;
	switch( (x + 100 + (y + 100) * 201 + (z+50) * 201 * 201) % 5   )
	{
	  case 0:   p = cube;
	    break;
	  case 1:   p = sphere;
	    break;
	  case 2:   p = sphere_2;
	    break;
	  case 3:   p = parallelipipede;
	    break;
	  case 4:   p = taurus;
	    break;
	}
	 
	container->add_child(p, C, zero_center, object_edgeSizeOverall, container_EdgeSizeOverall);
      }
    }
  }
  
  container->recalc_mipmap_color();
  cout << "container-color" << color_of(sphere_2	) <<  endl;

  
  {	
    int k; 
    printf("paused, press a key"); scanf("%c", &k);  	    
  }    
  
  
  dur.stop();
  printf("\n");
  printf("DT = %lf, %i Objects added\n", dur.get().count(), nO);
  getrusage(RUSAGE_SELF, &_ru); cout << "Current memory usage = " <<  _ru.ru_maxrss / 1024.0 << " MB" <<  endl;
  
  
/**/
  printf("Debug 2.1\n"); fflush(stdout);  
  
  // CMTRenderer<BaseType, TFloat, 3, LosType, bNormalizeLOS> mtRenderer(camera, bitmap);  
  
  printf("Debug 2.2\n"); fflush(stdout);  
  //mtRenderer.Init(8);  
  
  dur.clear(); dur.start();
  for (int i = 24; i < 50; i++)
  {
    printf("rendering image %i\n", i);
    
    Vec<TFloat, 3> cameraCoords;
    cameraCoords[0] = 0;
    cameraCoords[1] = 0;
    cameraCoords[2] = -500 * object_edgeSizeOverall * 3 + (i * 0.1) * object_edgeSizeOverall;
    
    // mtRenderer.Render(camera, bitmap);
    std::vector<std::unique_ptr<std::thread>> array_threads;
    std::atomic<uint32_t> scan_line;
    scan_line.store(0);
    
    array_threads.resize(8);
    for (int i = 0; i < array_threads.size(); i++)
    {
      array_threads[i].reset(new std::thread(
	[camera, bitmap, &scan_line, container, &cameraCoords, &zero_center, container_EdgeSizeOverall]
	{
	  auto resX = camera->GetResolution()[0];
	  auto resY = camera->GetResolution()[1];	    
	  auto __bitmap = bitmap;
	  auto Y = scan_line.fetch_add(16);
	  
	  for (; Y < resY; Y = scan_line.fetch_add(16))
	  {
	    for (int i = 0; i < 16 && Y + i < resY; i++)
	    {
	      register uint32_t offs = (Y + i) * resX;
	      register Color<uint8_t>* bitmap = __bitmap + offs;
	      register uint32_t offsMax = offs + resX;    
	      auto ssp = &camera->m_arrayCameraPoints[offs];	    
	      
	      //cout << "Y= " << Y + i << " ; ";

	      for ( ; offs < offsMax; offs++, bitmap++, ssp++ )
	      {
		RenderStruct<TFloat, 3, Color<uint8_t>> rs;
		rs.m_screenPoint = cameraCoords + ssp->spacePoint;
		rs.m_LOS = ssp->lineOfSight;
		//double z = ~(ssp->lineOfSight);
		//cout << "~(ssp->lineOfSight)=" << z << endl; 
		rs.m_LOS.normalize();
		rs.m_solid_angle_radius.m_b = ssp->solidAngleRadiusB /*/ z*/;
		rs.m_solid_angle_radius.m_a = ssp->solidAngleRadiusA /*/ z*/; // the eye point is at 10 units behing the screenPoint
		
		if (render(rs, *container, zero_center, container_EdgeSizeOverall))
		  *bitmap = rs.m_color;
	      }
	    } // for i	    
	    //cout << endl;
 	  }	  
	}
      ));
    }
    
    for (int i = 0; i < array_threads.size(); i++)
      array_threads[i]->join();
    
    //if (false)
    {
      BITMAPINFO bitmapinfo;
      memset(&bitmapinfo, 0, sizeof(bitmapinfo));
      bitmapinfo.bmiHeader.biSize   = sizeof(bitmapinfo);
      bitmapinfo.bmiHeader.biWidth  = camera->GetResolution()[0];
      bitmapinfo.bmiHeader.biHeight = camera->GetResolution()[1];
      bitmapinfo.bmiHeader.biPlanes = 1;
      bitmapinfo.bmiHeader.biBitCount = sizeof(Color<uint8_t>)*8;
      bitmapinfo.bmiHeader.biSizeImage = bitmapinfo.bmiHeader.biWidth * bitmapinfo.bmiHeader.biHeight * (bitmapinfo.bmiHeader.biBitCount / 8);
     
      char fn[128];
      sprintf(fn, "result-%03i.bmp", i);
      SaveDIBitmap(fn, &bitmapinfo, (char*)bitmap);
    }
  }
  //mtRenderer.Uninit();
  //now = time(0);
  dur.stop();
  printf("DT = %lf\n", dur.get().count());
  
  getrusage(RUSAGE_SELF, &_ru);
  printf("Current memory usage = %lf MB\n",  _ru.ru_maxrss / 1024.0);
  
  delete bitmap;
  delete camera;
  delete container;
  delete cube 		;
  delete sphere 	;
  delete taurus 	;
  delete parallelipipede;
  delete sphere_2 	;
  
  {	int k; scanf("%i", &k);  	}
  return 0;
}



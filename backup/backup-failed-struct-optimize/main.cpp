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

#include "vector/vector.h"
#include "line/line.h"
#include "container/container.h"
#include <chrono>

using namespace std;
typedef double TFlot;

int main(int argc, char** argv)
{
  Vec<TFlot, 3> v1 {1, 2, 3}, v2 {4, 5, 6};
  
  cout << v1 << endl;
  cout << v2 << endl;
  cout << v1*v2 << endl;
  cout << v1%v2 << endl;
  cout << (v1^v2) << endl;
  
  Line<TFlot, 3> l;
  l.m_org = {0, 0, 0};
  l.m_dir = {1, 1, 0};
  l.m_dir.normalize();

  PlaneDirs<TFlot, 3> d;
  d.m_dirs[0] = {0, 1, 0};
  d.m_dirs[1] = {0, 0, 1};
  d.calc_normal();
  
  Plane<TFlot, 3> p(Vec<TFlot, 3>{4, 0, 0}, &d);
  
  auto x = p & l;
  auto pnt = l.m_org + x * l.m_dir;
  cout << pnt << endl;
  
  RenderStruct<TFlot, 3, Color> rs;
  rs.m_screenPoint = {0, 0, 0};
  rs.m_LOS = {1, 0, 0};
  
  Color c;
  c.c.u32 = 0x000000FF;
  
  constexpr uint32_t Log2 = 10;
  
  auto container = new Container<uint32_t, 3, Log2, Color>;
  for (int y = 0; y < (1u << Log2); y++)
  for (int z = 0; z < (1u << Log2); z++)
    container->add_child(Vec<uint32_t, 3>{(1 << Log2)/2, y, z}, c);
 
  using namespace std::chrono;
  high_resolution_clock::time_point t1 = high_resolution_clock::now();

  for (int i = 0; i < 1000000; i++)
    render(rs, *container, Vec<TFlot, 3>{2500, 250, 250}, 10.0 * (1u << Log2));
  
  cout << "rendering done" << endl;
  
  high_resolution_clock::time_point t2 = high_resolution_clock::now();

  duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

  cout << "It took me " << time_span.count() << " seconds." << endl;
  
  cout << "found color " << rs.m_color.c.u32 << endl;
  
  int i;
  cin >> i;
  
  return 0;
}
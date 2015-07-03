#ifndef SAFE_MAP_READ_HPP
#define SAFE_MAP_READ_HPP 1

#include <utility>
#include <map>

using std::map;

template<class S, class T> const T& safe_map_read(const map<S,T>& m,
						  const S& s)
{
  assert(m.count(s)==1);
  return m.find(s)->second;
}

#endif // SAFE_MAP_READ_HPP

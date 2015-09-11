template<class S, class T> const T& safe_retrieve
(const boost::container::flat_map<S,T>& m,
 const S& s)
{
  assert(m.find(s)!=m.end());
  return m.find(s)->second;
}

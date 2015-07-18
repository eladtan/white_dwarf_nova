#ifndef MULTIPLE_MANIPULATION_HPP
#define MULTIPLE_MANIPULATION_HPP 1

#include "source/newtonian/test_2d/main_loop_2d.hpp"

class MultipleManipulation: public Manipulate
{
public:

  MultipleManipulation(const vector<Manipulate*>& manip_list);

  void operator()(hdsim& sim);

  ~MultipleManipulation(void);

private:
  const vector<Manipulate*> manip_list_;
};

#endif // MULTIPLE_MANIPULATION_HPP

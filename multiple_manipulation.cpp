#include "multiple_manipulation.hpp"

MultipleManipulation::MultipleManipulation
(const vector<Manipulate*>& manip_list):
  manip_list_(manip_list) {}

void MultipleManipulation::operator()(hdsim& sim)
{
  for(size_t i=0;i<manip_list_.size();++i)
    (*manip_list_[i])(sim);
}

MultipleManipulation::~MultipleManipulation(void)\
{
  for(size_t i=0;i<manip_list_.size();++i)
    delete manip_list_[i];
}

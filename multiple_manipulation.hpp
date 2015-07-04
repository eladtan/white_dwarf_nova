#ifndef MULTIPLE_MANIPULATION_HPP
#define MULTIPLE_MANIPULATION_HPP 1

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

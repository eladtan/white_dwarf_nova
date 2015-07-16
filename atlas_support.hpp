#ifndef ATLAS_SUPPORT
#define ATLAS_SUPPORT 1

#include "source/newtonian/test_2d/main_loop_2d.hpp"

class AtlasSupport: public Manipulate
{
public:

  AtlasSupport(void);

  void operator()(hdsim& sim);
};

#endif // ATLAS_SUPPORT

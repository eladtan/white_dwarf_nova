#include "bracketed.hpp"

bool bracketed(int low, int arg, int high)
{
  return arg>=low && high>arg;
}

#ifndef COSMOCALC_H_SEEN
#define COSMOCALC_H_SEEN

#include "luc.h"
#include "interpcosmoc.h"
#include "constcosmo.h"

inline double ZErr2CoDistErr(SimpleUniverse su, double zerr)
{ 
  return zerr * (1+su.ZE()) * ( (SPEED_OF_LIGHT_MS)/(su.HZE()*1e3) );
}

#endif

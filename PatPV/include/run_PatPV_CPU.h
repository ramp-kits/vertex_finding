#pragma once

#ifndef PATPV_RUN_PATPV_CPU_H
#define PATPV_RUN_PATPV_CPU_H

#include "AdaptivePV3DFitter.h"
#include "patPV_Definitions.cuh"
#include <algorithm>
#include <string>

extern "C" {

uint32_t reconstructMultiPVFromTracks(VeloState tracks2use[], VertexBase outvtxvec[], 
                                      bool tracks2disable[], XYZPoint seeds[], 
                                      uint32_t number_of_tracks) ;

}

#endif

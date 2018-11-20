
#ifndef PATPV_ADAPTIVEPV3DFITTER_H
#define PATPV_ADAPTIVEPV3DFITTER_H


#include "patPV_Definitions.cuh"






extern "C" {
// Fitting
bool fitVertex( XYZPoint& seedPoint,
            VeloState host_velo_states[],
           Vertex& vtx, uint32_t number_of_tracks, bool tracks2disable[], bool tracks2remove[]) ;

// Get Tukey's weight
double getTukeyWeight(double trchi2, int iter) ;
}

#endif // ADAPTIVE_H

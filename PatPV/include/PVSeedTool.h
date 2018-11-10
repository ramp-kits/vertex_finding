#pragma once

#ifndef PATPV_PVSEEDTOOL_H
#define PATPV_PVSEEDTOOL_H 



#include "patPV_Definitions.cuh"


/** @class PVSeedTool PVSeedTool.h tmp/PVSeedTool.h
 *
 *
 *  @author Mariusz Witek
 *  @date   2005-11-19
 */


extern "C" {

uint32_t getSeeds( VeloState inputTracks[], const XYZPoint& beamspot, 
                   XYZPoint  seeds[], bool tracks2disable[],
                   uint32_t number_of_tracks)  ;

int findClusters(vtxCluster * vclus, double * zclusters, int number_of_clusters) ;
void errorForPVSeedFinding(double tx, double ty, double &sigzaq) ;

double zCloseBeam(  VeloState track, const XYZPoint& beamspot) ;

}

#endif //PATPV_PVSEEDTOOL_H 

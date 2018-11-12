#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <typeinfo>
#include "readData.h"
#include "PVChecker.h"
#include "PVSeedTool.h"
#include "AdaptivePV3DFitter.h"
#include "run_PatPV_CPU.h"
#include "patPV_Definitions.cuh"

using namespace PatPV;

int main (int argc, char *argv[]) {
  signed char c;

  char* filename = "../data/RapidVPData_6719289_74676.json";
  while ((c = getopt(argc, argv, "f:")) != -1) {
    switch (c) {
    case 'f':
      filename = (char*)(optarg);
      break;
    default:
      return -1;
    }
  }

  //allocate arrays
  VeloState host_velo_states[max_tracks];
  
  //seeds
  XYZPoint  seeds[max_seeds];
  //output vertices
  VertexBase outvtxvec[max_vertices];
  //vector containing all mc vertices of all event-> one entry in vector is vector of mc vertices of that event
  MCVertex mcvertices_all[max_vertices];

  uint32_t number_of_tracks = readTracks(host_velo_states, 
                                         filename);
  uint32_t num_mc_vertices  = readMCVertices(mcvertices_all, 
                                            filename);

  std::cout << "---- begin fitter --------- " << std::endl;

  bool  tracks2disable[max_tracks];

  //stride because we saved tracks from straight line fit and Kalman fit in same vector
  for(int i = 0; i < number_of_tracks; i++) tracks2disable[i] = false;
  
  //reconstruct PVs 
  uint32_t num_rec_vertices = reconstructMultiPVFromTracks(host_velo_states, outvtxvec, 
                                                           tracks2disable, seeds, number_of_tracks);
  //check PVs
  checkPVs(mcvertices_all,
           num_mc_vertices,
           outvtxvec,
           num_rec_vertices);

  return 0;
}


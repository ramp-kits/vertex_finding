#include "run_PatPV_CPU.h"
#include "PVSeedTool.h"
#include <cmath>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;
//configuration

XYZPoint beamspot = {0.,0.,0.};


uint32_t reconstructMultiPVFromTracks(VeloState tracks2use[], VertexBase outvtxvec[], 
                                  bool tracks2disable[], XYZPoint seeds[],
                                  uint32_t number_of_tracks) 
{

  pt::ptree configtree;
  pt::read_json("config.json",configtree);

  //separation to accept reconstructed PV
  const auto  m_pvsChi2Separation = configtree.get("pvsChi2Separation",25.);
  const auto  m_pvsChi2SeparationLowMult = configtree.get("pvsChi2SeparationLowMult",91);


  VeloState * rtracks = tracks2use; 

  // reconstruct vertices

  uint32_t nvtx_after  =  0;
  bool continue_fitting = true;

  while(continue_fitting) {

    uint32_t number_of_seeds = getSeeds( rtracks, beamspot, seeds, tracks2disable, number_of_tracks);
    uint32_t before_fit = nvtx_after;
    for(int i=0; i < number_of_seeds; i++) {
      XYZPoint seed = seeds[i]; 
      Vertex recvtx;

      // fitting

      bool tracks2remove[number_of_tracks];

      for(int i = 0; i < number_of_tracks; i++) tracks2remove[i] = false;


      bool scvfit = fitVertex( seed, rtracks, recvtx, number_of_tracks, tracks2disable, tracks2remove);
      if (!scvfit) continue;


      //only accept vertex if it is not too close too already found one

   
      double chi2min = 1e10;
      for(int i = 0; i < nvtx_after; i++) {
        int index = i;
        double z1 = outvtxvec[index].z;
        double z2 = recvtx.vb.z;
        double sigma2z1 = outvtxvec[index].c22;
        double sigma2z2 = recvtx.vb.c22;
        double chi2 = pow((z1 - z2), 2) / (sigma2z1 + sigma2z2);
        if (chi2 < chi2min) chi2min = chi2;
      }

      bool vsepar = true;

      if ( chi2min < m_pvsChi2Separation ) vsepar = false;
  // protect secondary vertices of B signal
      //if ( chi2min < m_pvsChi2SeparationLowMult && recvtx.ndof < 7 ) vsepar = false;
      if ( chi2min < m_pvsChi2SeparationLowMult && 0.5*(recvtx.vb.ndof+3) < 7 ) vsepar = false;

      if(!vsepar) continue;

      //remove tracks used

      for(int i = 0; i < number_of_tracks; i++) {
            if(tracks2remove[i]) tracks2disable[i] = true;
      }

      //add reconstructed vertex to output array
      outvtxvec[nvtx_after] = recvtx.vb;
      nvtx_after++;

    } //iterate on seeds

    if(before_fit == nvtx_after) continue_fitting = false;
    
    
  }

  return nvtx_after;
  

}

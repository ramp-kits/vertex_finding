#ifndef PATPV_PVCHECKER_H
#define PATPV_PVCHECKER_H

#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>
#include "patPV_Definitions.cuh"


extern "C" {

//configuration for PV checker -> check values
int m_nTracksToBeRecble = 4;
double m_dzIsolated = 10; //mm
bool m_matchByTracks = false;

void checkPVs( MCVertex mcvertices_all[],
               uint32_t num_mc_vertices, 
               VertexBase rec_vertex[],
               uint32_t num_rec_vertices);

void match_mc_vertex_by_distance(int ipv, std::vector<RecPVInfo>& rinfo, std::vector<MCPVInfo>& mcpvvec) {

  double mindist = 999999.;
  int indexmc = -1;

  for(int imc = 0; imc < (int) mcpvvec.size(); imc++) {
    if ( mcpvvec[imc].indexRecPVInfo  > -1) continue;
    double dist = fabs(mcpvvec[imc].pMCPV->z -
                       rinfo[ipv].z);
    if(dist < mindist) {
      mindist = dist;
      indexmc = imc;
    }
  }
  if ( indexmc > -1 ) {
    if(mindist < 5.0 * rinfo[ipv].positionSigma.z) {
      rinfo[ipv].indexMCPVInfo = indexmc;
      mcpvvec[indexmc].indexRecPVInfo = ipv;
    }
  }

}

void printRat(std::string mes, int a, int b) {

  double rat = 0.;
  if(b>0) rat = 1.0*a/b;

  // reformat message
  unsigned int len = 20;
  std::string pmes = mes;
  while(pmes.length() < len) {
    pmes+=" ";
  }
  pmes+= " : ";

  std::cout << pmes << " " << rat << "( " << a << " / " << b << " )" << std::endl;

}

std::vector<MCPVInfo>::iterator closestMCPV(std::vector<MCPVInfo>& rblemcpv,
                                            std::vector<MCPVInfo>::iterator& itmc) {

  std::vector<MCPVInfo>::iterator itret = rblemcpv.end();
  double mindist = 999999.;
  if(rblemcpv.size() < 2) return itret;
  std::vector<MCPVInfo>::iterator it;
  for (it = rblemcpv.begin(); it != rblemcpv.end(); it++) {
    if(it->pMCPV != itmc->pMCPV ) {
      double diff_x = it->pMCPV->x - itmc->pMCPV->x;
      double diff_y = it->pMCPV->y - itmc->pMCPV->y;
      double diff_z = it->pMCPV->z - itmc->pMCPV->z;
      double dist = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);
      
      if(dist < mindist) {
        mindist = dist;
        itret = it;
      }
    }
  }
  return itret;
}

}

#endif // 

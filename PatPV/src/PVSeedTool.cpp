// Include files

// STL
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

// local
#include "PVSeedTool.h"


bool  vtxcomp( vtxCluster *first, vtxCluster *second ) {
    return first->z < second->z;
}
bool  multcomp( vtxCluster *first, vtxCluster *second ) {
    return first->ntracks > second->ntracks;
}


// steering parameters for merging procedure
double m_maxChi2Merge = 25.;
double m_factorToIncreaseErrors = 15.;

//try parameters from RecoUpgradeTracking.py
int    m_minClusterMult = 4;
int    m_minCloseTracksInCluster = 3;


// steering parameters for final cluster selection
// int    m_minClusterMult = 3;
double m_dzCloseTracksInCluster = 5.; // unit: mm
// int    m_minCloseTracksInCluster = 3;
int    m_highMult = 10;
double m_ratioSig2HighMult = 1.0;
double m_ratioSig2LowMult = 0.9;

int m_max_clusters = 200; // maximmum number of clusters

double m_x0MS = 0.01;// X0 (tunable) of MS to add for extrapolation of
                                                       // track parameters to PV

//don't forget to actually calculate this!!
//double  m_scatCons = 0;     // calculated from m_x0MS
double X0 = m_x0MS;
double m_scatCons = (13.6*sqrt(X0)*(1.+0.038*log(X0)));






//=============================================================================
// getSeeds
//=============================================================================
uint32_t getSeeds( VeloState inputTracks[], const XYZPoint& beamspot,
                   XYZPoint seeds[], bool tracks2disable[],
                   uint32_t number_of_tracks)  {
  
  
  

  vtxCluster  vclusters[1200];

  int number_of_clusters = 0;
  for (int i = 0; i < number_of_tracks; i++) {

    //don't reuse used tracks
    if (tracks2disable[i]) continue;
    double sigsq;
    double zclu;
    auto trk = inputTracks[i];
   
    zclu = zCloseBeam(trk,beamspot);
    errorForPVSeedFinding(trk.tx, trk.ty,sigsq);

    if ( fabs(zclu)>2000.) continue;
    vtxCluster clu;
    clu.z = zclu;
    clu.sigsq = sigsq;
    clu.sigsqmin = clu.sigsq;
    clu.ntracks = 1;
    vclusters[number_of_clusters] = clu;
    
    number_of_clusters++;

  }

  double  zseeds[number_of_clusters];

  int number_final_clusters = findClusters(vclusters, zseeds, number_of_clusters);

  for(int i = 0; i < number_final_clusters; i++) seeds[i] = XYZPoint{ beamspot.x, beamspot.y, zseeds[i]};
  

  return number_final_clusters;

}


int findClusters(vtxCluster * vclus, double * zclusters, int number_of_clusters)  {



  

  
  for(int i = 0; i < number_of_clusters; i++) {
    vclus[i].sigsq *= m_factorToIncreaseErrors*m_factorToIncreaseErrors; // blow up errors
    vclus[i].sigsqmin = vclus[i].sigsq;
  }


  //maybe sort in z before merging? -> does not seem to help

  bool no_merges = false;
  while(!no_merges) {
    //reset merged flags
    for (int j = 0; j < number_of_clusters; j++) vclus[j].merged = false;

    no_merges = true;
    for(int index_cluster = 0; index_cluster < number_of_clusters - 1; index_cluster++) {


      //skip cluster which have already been merged
      if(vclus[index_cluster].ntracks == 0) continue;

   //sorting by chi2dist seems to increase efficiency in nominal code


      
      for(int index_second_cluster = 0; index_second_cluster < number_of_clusters ; index_second_cluster++) {
        if(vclus[index_second_cluster].merged || vclus[index_cluster].merged) continue;
        //skip cluster which have already been merged
        if(vclus[index_second_cluster].ntracks == 0) continue;
        if(index_cluster == index_second_cluster) continue;
        double z1 = vclus[index_cluster].z;
        double z2 = vclus[index_second_cluster].z;
        double s1 = vclus[index_cluster].sigsq;
        double s2 = vclus[index_second_cluster].sigsq;
        double s1min = vclus[index_cluster].sigsqmin;
        double s2min = vclus[index_second_cluster].sigsqmin;
        double sigsqmin = s1min;
        if(s2min<s1min) sigsqmin = s2min;


        double zdist = z1 - z2;
        double chi2dist = zdist*zdist/(s1+s2);
        //merge if chi2dist is smaller than max
        if (chi2dist<m_maxChi2Merge ) {
          no_merges = false;
          double w_inv = (s1*s2/(s1+s2));
          double zmerge = w_inv*(z1/s1+z2/s2);

          vclus[index_cluster].z        = zmerge;
          vclus[index_cluster].sigsq    = w_inv;
          vclus[index_cluster].sigsqmin = sigsqmin;
          vclus[index_cluster].ntracks += vclus[index_second_cluster].ntracks;
          vclus[index_second_cluster].ntracks = 0;  // mark second cluster as used
          vclus[index_cluster].merged = true;
          vclus[index_second_cluster].merged = true;

          //break;
        } 
      }
    }
  }
 
  

  int return_number_of_clusters = 0;
  //count final number of clusters
  vtxCluster pvclus[2000];
  for(int i = 0; i < number_of_clusters; i++) {
    if(vclus[i].ntracks != 0)     {pvclus[return_number_of_clusters] = vclus[i]; return_number_of_clusters++;}
  } 


  //clean up clusters, do we gain much from this?

  // Select good clusters.

  int number_good_clusters = 0;

  for(int index = 0; index < return_number_of_clusters; index++) {

    int n_tracks_close = 0;
    for(int i = 0; i < number_of_clusters; i++) if(fabs(vclus[i].z - pvclus[index].z ) < m_dzCloseTracksInCluster ) n_tracks_close++;
  

    double dist_to_closest = 1000000.;
    if(return_number_of_clusters > 1) {
      for(int index2 = 0; index2 < return_number_of_clusters; index2++) {
        if( index!=index2 && ( fabs( pvclus[index2].z - pvclus[index].z) < dist_to_closest) )  dist_to_closest = fabs( pvclus[index2].z - pvclus[index].z);
      }
    }

    // ratio to remove clusters made of one low error track and many large error ones
    double rat = pvclus[index].sigsq/pvclus[index].sigsqmin;
    bool igood = false;
    int ntracks = pvclus[index].ntracks;
    if( ntracks >= m_minClusterMult ) {
      if( dist_to_closest>10. && rat<0.95) igood=true;
      if( ntracks >= m_highMult && rat < m_ratioSig2HighMult)  igood=true;
      if( ntracks <  m_highMult && rat < m_ratioSig2LowMult )  igood=true;
    }
    // veto
    if( n_tracks_close < m_minCloseTracksInCluster ) igood = false;
    if(igood) {zclusters[number_good_clusters] = pvclus[index].z; number_good_clusters++;}


  }



  //return return_number_of_clusters;
  return number_good_clusters;

}

void errorForPVSeedFinding(double tx, double ty, double &sigz2)  {

    // the seeding results depend weakly on this eror parametrization

    double pMean = 3000.; // unit: MeV

    double tanTheta2 =  tx * tx + ty * ty;
    double sinTheta2 =  tanTheta2 / ( 1. + tanTheta2 );

    // assume that first hit in VD at 8 mm
    double distr        = 8.; // unit: mm
    double dist2        = distr*distr/sinTheta2;
    double sigma_ms2    = m_scatCons * m_scatCons * dist2 / (pMean*pMean);
    double fslope2      = 0.0005*0.0005;
    double sigma_slope2 = fslope2*dist2;

    sigz2 = (sigma_ms2 + sigma_slope2) / sinTheta2;
    if(sigz2 == 0) sigz2 = 100.;

}



double zCloseBeam( VeloState track, const XYZPoint& beamspot) {

  XYZPoint tpoint;
  tpoint.x = track.x;
  tpoint.y = track.y;
  tpoint.z = track.z;
  XYZPoint tdir;
  tdir.x = track.tx;
  tdir.y = track.ty;
  tdir.z = 1.;

  double wx = ( 1. + tdir.x * tdir.x ) / track.c00;
  double wy = ( 1. + tdir.y * tdir.y ) / track.c11;

  double x0 = tpoint.x - tpoint.z * tdir.x - beamspot.x;
  double y0 = tpoint.y - tpoint.z * tdir.y - beamspot.y;
  double den = wx * tdir.x * tdir.x + wy * tdir.y * tdir.y;
  double zAtBeam = - ( wx * x0 * tdir.x + wy * y0 * tdir.y ) / den ;

  double xb = tpoint.x + tdir.x * ( zAtBeam - tpoint.z ) - beamspot.x;
  double yb = tpoint.y + tdir.y * ( zAtBeam - tpoint.z ) - beamspot.y;
  double r2AtBeam = xb*xb + yb*yb ;

  return r2AtBeam < 0.5*0.5 ? zAtBeam : 10e8;
}

//=============================================================================
 

#include <vector>
#include <cmath>
#include <iostream>

#include "AdaptivePV3DFitter.h"


//configuration
size_t m_minTr = 4;
int    m_Iterations = 20;
int    m_minIter = 5;
double m_maxDeltaZ = 0.0005; // unit:: mm
double m_minTrackWeight = 0.00000001;
double m_TrackErrorScaleFactor = 1.0;
double m_maxChi2 = 400.0;
double m_trackMaxChi2 = 12.;
double m_trackChi = std::sqrt(m_trackMaxChi2);     // sqrt of trackMaxChi2
double m_trackMaxChi2Remove = 25.;
double m_maxDeltaZCache = 1.; //unit: mm



//=============================================================================
// Least square adaptive fitting method
//=============================================================================
bool fitVertex( XYZPoint& seedPoint, VeloState host_velo_states[], Vertex& vtx,
                uint32_t number_of_tracks, bool tracks2disable[], bool tracks2remove[]) 
{




  double tr_state_x[number_of_tracks] ;
  double tr_state_y[number_of_tracks] ;
  double tr_state_z[number_of_tracks] ;

  double tr_state_tx[number_of_tracks];
  double tr_state_ty[number_of_tracks] ;

  double tr_state_c00[number_of_tracks] ;
  double tr_state_c11[number_of_tracks] ;
  double tr_state_c20[number_of_tracks] ;
  double tr_state_c22[number_of_tracks] ;
  double tr_state_c31[number_of_tracks] ;
  double tr_state_c33[number_of_tracks] ;

  // position at which derivatives are evaluated

  XYZPoint vtxpos = seedPoint ;

  // prepare tracks
 

  int pvTrack_counter = 0;

  for(int index = 0; index < number_of_tracks; index++) {  
    //don't use disabled tracks
    if(tracks2disable[index]) continue;

    double new_z = vtxpos.z;

    double m_state_x = host_velo_states[index].x;
    double m_state_y = host_velo_states[index].y;
    double m_state_z = host_velo_states[index].z;

    double m_state_tx = host_velo_states[index].tx;
    double m_state_ty = host_velo_states[index].ty;

    double m_state_c00 = host_velo_states[index].c00;
    double m_state_c11 = host_velo_states[index].c11;
    double m_state_c20 = host_velo_states[index].c20;
    double m_state_c22 = host_velo_states[index].c22;
    double m_state_c31 = host_velo_states[index].c31;
    double m_state_c33 = host_velo_states[index].c33;

    const double dz = new_z - m_state_z ;
    const double dz2 = dz*dz ;

    m_state_x += dz * m_state_tx ;
    m_state_y += dz * m_state_ty ;
    m_state_z = new_z;
    m_state_c00 += dz2 * m_state_c22 + 2*dz* m_state_c20 ;
    m_state_c20 += dz* m_state_c22 ;
    m_state_c11 += dz2* m_state_c33 + 2* dz*m_state_c31 ;
    m_state_c31 += dz* m_state_c33 ;

    Vector2 res{ vtxpos.x - m_state_x, vtxpos.y - m_state_y };


    double  tr_chi2 = res.x*res.x / m_state_c00 +res.y*res.y / m_state_c11;
      


    if(tr_chi2 < m_maxChi2) {
 // have to use updated values!!
      tr_state_x[pvTrack_counter] = m_state_x;
      tr_state_y[pvTrack_counter] = m_state_y;
      tr_state_z[pvTrack_counter] = m_state_z;

      tr_state_tx[pvTrack_counter] = m_state_tx;
      tr_state_ty[pvTrack_counter] = m_state_ty;

      tr_state_c00[pvTrack_counter] = m_state_c00;
      tr_state_c11[pvTrack_counter] = m_state_c11;
      tr_state_c20[pvTrack_counter] = m_state_c20;
      tr_state_c22[pvTrack_counter] = m_state_c22;
      tr_state_c31[pvTrack_counter] = m_state_c31;
      tr_state_c33[pvTrack_counter] = m_state_c33;

      pvTrack_counter++;
    }

      
  }
    

  if( pvTrack_counter < m_minTr ) {
    std::cout << "Too few tracks to fit PV" << std::endl;
    return false;
  }


  double vtxcov[6] ;
  bool converged = false;
  double maxdz = m_maxDeltaZ;
  int nbIter = 0;
  int tracks_in_vertex = 0;
  while( (nbIter < m_minIter) || (!converged && nbIter < m_Iterations) ) {
    ++nbIter;

    double halfD2Chi2DX2_00 = 0.;
    double halfD2Chi2DX2_10 = 0.;
    double halfD2Chi2DX2_11 = 0.;
    double halfD2Chi2DX2_20 = 0.;
    double halfD2Chi2DX2_21 = 0.;
    double halfD2Chi2DX2_22 = 0.;
    XYZPoint halfDChi2DX;
    halfDChi2DX.x = 0.;
    halfDChi2DX.y = 0.;
    halfDChi2DX.z = 0.;


    
    

    // add contribution from all tracks
    double chi2(0) ;
    size_t ntrin(0) ;
    for( int index = 0; index < pvTrack_counter; index++) {
      //update cache

      double new_z = vtxpos.z;

      

      double m_state_x = tr_state_x[index];
      double m_state_y = tr_state_y[index];
      double m_state_z = tr_state_z[index];

      double m_state_tx = tr_state_tx[index];
      double m_state_ty = tr_state_ty[index];

      double m_state_c00 = tr_state_c00[index];
      double m_state_c11 = tr_state_c11[index];
      double m_state_c20 = tr_state_c20[index];
      double m_state_c22 = tr_state_c22[index];
      double m_state_c31 = tr_state_c31[index];
      double m_state_c33 = tr_state_c33[index];

      const double dz = new_z - m_state_z ;
      const double dz2 = dz*dz ;

      m_state_x += dz * m_state_tx ;
      m_state_y += dz * m_state_ty ;
      m_state_z = new_z;
      m_state_c00 += dz2 * m_state_c22 + 2*dz* m_state_c20 ;
      m_state_c20 += dz* m_state_c22 ;
      m_state_c11 += dz2* m_state_c33 + 2* dz*m_state_c31 ;
      m_state_c31 += dz* m_state_c33 ;

      Vector2 res{ vtxpos.x - m_state_x, vtxpos.y - m_state_y };

      double tr_halfD2Chi2DX2_00 = 1. / m_state_c00;
      double tr_halfD2Chi2DX2_10 = 0.;
      double tr_halfD2Chi2DX2_11 = 1. / m_state_c11;
      double tr_halfD2Chi2DX2_20 = - m_state_tx / m_state_c00;
      double tr_halfD2Chi2DX2_21 = - m_state_ty / m_state_c11;
      double tr_halfD2Chi2DX2_22 = m_state_tx * m_state_tx / m_state_c00 + m_state_ty * m_state_ty / m_state_c11;

      double tr_halfDChi2DX_x = res.x / m_state_c00;
      double tr_halfDChi2DX_y = res.y / m_state_c11;
      double tr_halfDChi2DX_z = -m_state_tx*res.x / m_state_c00-m_state_ty*res.y / m_state_c11;
      double tr_chi2          = res.x*res.x / m_state_c00 +res.y*res.y / m_state_c11;


      double weight = getTukeyWeight(tr_chi2, nbIter) ;

      // add the track
      if ( weight > m_minTrackWeight ) {
        ++ntrin;
        halfD2Chi2DX2_00 += weight * tr_halfD2Chi2DX2_00 ;
        halfD2Chi2DX2_10 += weight * tr_halfD2Chi2DX2_10 ;
        halfD2Chi2DX2_11 += weight * tr_halfD2Chi2DX2_11 ;
        halfD2Chi2DX2_20 += weight * tr_halfD2Chi2DX2_20 ;
        halfD2Chi2DX2_21 += weight * tr_halfD2Chi2DX2_21 ;
        halfD2Chi2DX2_22 += weight * tr_halfD2Chi2DX2_22 ;



        halfDChi2DX.x   += weight * tr_halfDChi2DX_x ;
        halfDChi2DX.y   += weight * tr_halfDChi2DX_y ;
        halfDChi2DX.z   += weight * tr_halfDChi2DX_z ;

        chi2 += weight * tr_chi2 ;


      }
    }
    
    // check nr of tracks that entered the fit
    if(ntrin < m_minTr) {
      std::cout << "Too few tracks after PV fit" << std::endl;
      return false;
    }

    // compute the new vertex covariance


    //repalce Cholesky inverter by analytical solution
    double a00 = halfD2Chi2DX2_00;
    double a10 = halfD2Chi2DX2_10;
    double a11 = halfD2Chi2DX2_11;
    double a20 = halfD2Chi2DX2_20;
    double a21 = halfD2Chi2DX2_21;
    double a22 = halfD2Chi2DX2_22;

    double det = a00 * (a22 * a11 - a21 * a21) - a10 * (a22 * a10 - a21 * a20) + a20 * (a21*a10 - a11*a20);
    if (det == 0) return false;


   vtxcov[0] = (a22*a11 - a21*a21) / det;
   vtxcov[1] = -(a22*a10-a20*a21) / det;
   vtxcov[2] = (a22*a00-a20*a20) / det;
   vtxcov[3] = (a21*a10-a20*a11) / det;
   vtxcov[4] = -(a21*a00-a20*a10) / det;
   vtxcov[5] = (a11*a00-a10*a10) / det;

    

    // compute the delta w.r.t. the reference
    XYZPoint delta{0.,0.,0.};
    delta.x = -1.0 * (vtxcov[0] * halfDChi2DX.x + vtxcov[1] * halfDChi2DX.y + vtxcov[3] * halfDChi2DX.z );
    delta.y = -1.0 * (vtxcov[1] * halfDChi2DX.x + vtxcov[2] * halfDChi2DX.y + vtxcov[4] * halfDChi2DX.z );
    delta.z = -1.0 * (vtxcov[3] * halfDChi2DX.x + vtxcov[4] * halfDChi2DX.y + vtxcov[5] * halfDChi2DX.z );

    // note: this is only correct if chi2 was chi2 of reference!
    chi2  += delta.x * halfDChi2DX.x + delta.y * halfDChi2DX.y + delta.z * halfDChi2DX.z;

    // deltaz needed for convergence
    const double deltaz = vtxpos.z + delta.z - vtxpos.z ;

    // update the position
    vtxpos.x = ( vtxpos.x + delta.x ) ;
    vtxpos.y = ( vtxpos.y + delta.y ) ;
    vtxpos.z = ( vtxpos.z + delta.z ) ;

    vtx.setChi2AndDoF( chi2, 2*ntrin-3 ) ;

    // loose convergence criteria if close to end of iterations
    if ( 1.*nbIter > 0.8*m_Iterations ) maxdz = 10.*m_maxDeltaZ;
    converged = std::abs(deltaz) < maxdz ;
    tracks_in_vertex = ntrin;
  } // end iteration loop

  if(!converged) return false;

  
  // set position and covariance
  vtx.setPosition( vtxpos ) ;
  vtx.setCovMatrix( vtxcov ) ;
  // Set tracks. Compute final chi2.
  vtx.clearTracks();




/*
  //radial cut against fakes?
      double m_beamSpotRCut  = 0.6;
    double m_beamSpotRCutHMC = 0.4;
    int m_beamSpotRMT =  10;
    std::cout << "tracks_in_vertex: " << tracks_in_vertex << std::endl;
    double r2 = std::pow(vtxpos.x- 0.,2) + std::pow( vtxpos.y- 0.,2);
       double r  = (  tracks_in_vertex <  m_beamSpotRMT ? m_beamSpotRCut : m_beamSpotRCutHMC );
        if ( r2 > r*r ) return false;
        */



  //disable tracks added to this vertex
  for(int index = 0; index < number_of_tracks; index++) {
    //don't use disabled tracks
    if(tracks2disable[index]) continue;

    double new_z = vtxpos.z;


    double m_state_x = host_velo_states[index].x;
    double m_state_y = host_velo_states[index].y;
    double m_state_z = host_velo_states[index].z;

    double m_state_tx = host_velo_states[index].tx;
    double m_state_ty = host_velo_states[index].ty;

    double m_state_c00 = host_velo_states[index].c00;
    double m_state_c11 = host_velo_states[index].c11;
    double m_state_c20 = host_velo_states[index].c20;
    double m_state_c22 = host_velo_states[index].c22;
    double m_state_c31 = host_velo_states[index].c31;
    double m_state_c33 = host_velo_states[index].c33;

    const double dz = new_z - m_state_z ;
    const double dz2 = dz*dz ;

    m_state_x += dz * m_state_tx ;
    m_state_y += dz * m_state_ty ;
    m_state_z = new_z;
    m_state_c00 += dz2 * m_state_c22 + 2*dz* m_state_c20 ;
    m_state_c20 += dz* m_state_c22 ;
    m_state_c11 += dz2* m_state_c33 + 2* dz*m_state_c31 ;
    m_state_c31 += dz* m_state_c33 ;

    Vector2 res{ vtxpos.x - m_state_x, vtxpos.y - m_state_y };

    double tr_chi2          = res.x*res.x / m_state_c00 + res.y*res.y / m_state_c11;

    if( tr_chi2 > m_minTrackWeight) {
      VeloState state = host_velo_states[index];
      vtx.addToTracks( state, tr_chi2 ) ;
    }

    //if( tr_chi2 < m_trackMaxChi2Remove) tracks2disable[index] = true;
    if( tr_chi2 < m_trackMaxChi2Remove) tracks2remove[index] = true;
   
  }
  

  
  return true;
}

//=============================================================================
// Get Tukey's weight
//=============================================================================

double getTukeyWeight(double trchi2, int iter) 
{
  if (iter<1 ) return 1.;
  double ctrv = m_trackChi * std::max(m_minIter -  iter,1);
  double cT2 = trchi2 / std::pow(ctrv*m_TrackErrorScaleFactor,2);
  return cT2 < 1. ? std::pow(1.-cT2,2) : 0. ;
}




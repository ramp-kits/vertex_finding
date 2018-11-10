#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <vector>
#include <stdint.h>
extern "C" {

namespace PatPV {
  static constexpr uint32_t max_tracks   = 1200;
  static constexpr uint32_t max_seeds    = 100;
  static constexpr uint32_t max_vertices = 30;
}

struct VeloState { // 48 B
  float x, y, z, tx, ty;
  float c00, c20, c22, c11, c31, c33;
  float chi2;
  bool backward;
};

struct vtxCluster final {

  double  z = 0;            // z of the cluster
  double  sigsq = 0;        // sigma**2 of the cluster
  double  sigsqmin = 0;     // minimum sigma**2 of the tracks forming cluster
  int     ntracks = 1;      // number of tracks in the cluster
  bool     merged = false;   // flag for iterative merging

  vtxCluster() = default;

};

struct XYZPoint {
  double x,y,z;
};

struct Vector2 {
  double x;
  double y;

  Vector2(double m_x, double m_y) : x(m_x), y(m_y){}
};

//Split Vertex into a simple VertexBase, to be accessed
//from python, and the more complex full Vertex, which
//we can leave for a later use as/when needed.
struct VertexBase{
  double x,y,z;
  double c00, c10, c11, c20, c21, c22;
  double chi2;
  int ndof;
};


class Vertex {
  public:
    VertexBase vb;    
    Vertex() {vb.x=0;vb.y=0;vb.z=0;
              vb.chi2=0;vb.ndof=0;
              vb.c00=0;vb.c10=0;vb.c11=0;
              vb.c20=0;vb.c21=0;vb.c22=0;}

    std::vector<VeloState> tracks;
    std::vector<double> weights;
    void setChi2AndDoF(double m_chi2, int m_ndof) {
      vb.chi2 = m_chi2;
      vb.ndof = m_ndof;
    }
    void setPosition(XYZPoint& point) {
      vb.x = point.x;
      vb.y = point.y;
      vb.z = point.z;
    }
    void setCovMatrix(double * m_cov) {
      vb.c00 = m_cov[0];
      vb.c10 = m_cov[1];
      vb.c11 = m_cov[2];
      vb.c20 = m_cov[3];
      vb.c21 = m_cov[4];
      vb.c22 = m_cov[5];
    }

    void clearTracks() {
      tracks.clear();
      weights.clear();
    };
    void addToTracks(VeloState track, double weight) {
      tracks.push_back(track);
      weights.push_back(weight);
    };
};

typedef struct {
 public:
  int nTracks;                    // number of tracks
  int nVeloTracks;                // number of velo tracks in a vertex       
  int nLongTracks;
  double minTrackRD;              //                                                                                                                                        
  double maxTrackRD;              //                                                                             
  double chi2;
  double nDoF;
  double d0;
  double d0nTr;
  double chi2nTr;
  double mind0;
  double maxd0;
  int mother;
  //XYZPoint position;       // position
  double x;
  double y;
  double z;
  XYZPoint positionSigma;  // position sigmas
  int indexMCPVInfo;              // index to MCPVInfo
  Vertex* pRECPV;        // pointer to REC PV
} RecPVInfo;

struct MCVertex {
  double x;
  double y;
  double z;
  int numberTracks;
};

typedef struct {  
  MCVertex* pMCPV;     // pointer to MC PV 
  int nRecTracks;            // number of reconstructed tracks from this MCPV  
  int nRecBackTracks;        // number of reconstructed backward tracks
  int indexRecPVInfo;        // index to reconstructed PVInfo (-1 if not reco)
  int nCorrectTracks;        // correct tracks belonging to reconstructed PV
  int multClosestMCPV;       // multiplicity of closest reconstructable MCPV
  double distToClosestMCPV;  // distance to closest reconstructible MCPV
  int decayCharm;             // type of mother particle                                                                                                                                           
  int decayBeauty;
  //std::vector<LHCb::MCParticle*> m_mcPartInMCPV;
  //std::vector<LHCb::Track*> m_recTracksInMCPV;
} MCPVInfo;

}
#endif // DEFINITIONS_H

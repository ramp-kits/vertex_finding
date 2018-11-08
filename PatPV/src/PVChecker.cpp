#include "PVChecker.h"

void checkPVs( MCVertex mcvertices_all[], 
               uint32_t num_mc_vertices, 
               VertexBase rec_vertex[],
               uint32_t num_rec_vertices)
{

  std::cout << "Checking PVs: " << std::endl;

  //counters for efficiency
  int number_reconstructible_vertices = 0; 
  int number_reconstructed_vertices = 0;
  int number_fake_vertices = 0;
    

  //counters for efficiencies/fake rate
  int m_nMCPV = 0;
  int m_nRecMCPV  = 0;
  int m_nMCPV_isol = 0;
  int m_nRecMCPV_isol = 0;
  int m_nMCPV_close = 0;
  int m_nRecMCPV_close = 0;
  int m_nFalsePV = 0;
  int m_nFalsePV_real = 0;

  //vectors to collect the pulls and erros
  std::vector<double> vec_diff_x;
  std::vector<double> vec_diff_y;
  std::vector<double> vec_diff_z;

  std::vector<double> vec_err_x;
  std::vector<double> vec_err_y;
  std::vector<double> vec_err_z;
  
  std::vector<Vertex*> vecOfVertices;
  //first fill vector with reconstructed vertices
  for(uint32_t i = 0; i < num_rec_vertices; i++) {
    Vertex* fullrecvtx = new Vertex();
    fullrecvtx->vb = rec_vertex[i];
    vecOfVertices.push_back(fullrecvtx);
  }
  // Fill reconstucted PV info
  std::vector<RecPVInfo> recpvvec;
  std::vector<Vertex*>::iterator itRecV;
  for(itRecV = vecOfVertices.begin(); vecOfVertices.end() != itRecV;
               itRecV++) {
    Vertex* pv;
    pv = *itRecV;
    RecPVInfo recinfo;
    recinfo.pRECPV= pv;
    recinfo.x = pv->vb.x;
    recinfo.y = pv->vb.y;
    recinfo.z = pv->vb.z;
    
    double sigx = sqrt(pv->vb.c00);
    double sigy = sqrt(pv->vb.c11);
    double sigz = sqrt(pv->vb.c22);
    XYZPoint a3d;
    a3d.x = sigx;
    a3d.y = sigy;
    a3d.z = sigz;
    recinfo.positionSigma = a3d;
    recinfo.nTracks = pv->tracks.size();
    double minRD = 99999.;
    double maxRD = -99999.;
    double chi2 = pv->vb.chi2;
    double nDoF = pv->vb.ndof;  

    int mother = 0;
    int velo = 0;
    int lg = 0;
    double d0 = 0;
    double mind0 = 99999.0;
    double maxd0 = -99999.0;
    double trackChi2 = 0.0;
    int tr = 0;

    recinfo.minTrackRD = minRD;
    recinfo.maxTrackRD = maxRD;
    recinfo.mother = mother;
    recinfo.chi2 = chi2;
    recinfo.nDoF = nDoF;
    recinfo.d0 = d0;
    recinfo.d0nTr = (double)d0/(double)tr;
    recinfo.chi2nTr = (double)trackChi2/(double)tr; 
    recinfo.mind0 = mind0;
    recinfo.maxd0 = maxd0;
    recinfo.nVeloTracks = velo; 
    recinfo.nLongTracks = lg; 
    recinfo.indexMCPVInfo = -1;
    recpvvec.push_back(recinfo);
    
  }
  
  //vector with MCPVinfo
  std::vector<MCPVInfo> mcpvvec;
                   
  for(uint32_t mc_vtx_offset = 0; mc_vtx_offset < num_mc_vertices; mc_vtx_offset++ ) {
    
        MCPVInfo mcprimvert;
        mcprimvert.pMCPV = &(mcvertices_all[mc_vtx_offset]);
        //mcprimvert.nRecTracks = 0;
        mcprimvert.nRecTracks = mcvertices_all[mc_vtx_offset].numberTracks;
        //mcprimvert.nRecTracks = 99;
        mcprimvert.nRecBackTracks = 0;
        mcprimvert.indexRecPVInfo = -1;
        mcprimvert.nCorrectTracks = 0;
        mcprimvert.multClosestMCPV = 0;
        mcprimvert.distToClosestMCPV = 999999.;
        mcprimvert.decayBeauty = 0;
        mcprimvert.decayCharm  = 0;
        
        mcpvvec.push_back(mcprimvert);
      
    
  } 

  std::vector<MCPVInfo> rblemcpv;
  std::vector<MCPVInfo> not_rble_but_visible;
  std::vector<MCPVInfo> not_rble;
  int nmrc = 0;

  //count not reconstructible MC PVs
  std::vector<MCPVInfo>::iterator itmc;
  for (itmc = mcpvvec.begin(); mcpvvec.end() != itmc; itmc++) {
    rblemcpv.push_back(*itmc);

    if (itmc->nRecTracks < m_nTracksToBeRecble)
      {
        nmrc++;
      }
    if(itmc->nRecTracks < m_nTracksToBeRecble && itmc->nRecTracks > 1)
      {
        not_rble_but_visible.push_back(*itmc);
      }
    if(itmc->nRecTracks < m_nTracksToBeRecble && itmc->nRecTracks < 2)
      {
        not_rble.push_back(*itmc);  
      }

  }

  //match by distance
  for(int ipv = 0; ipv < (int) recpvvec.size(); ipv++) {
      match_mc_vertex_by_distance(ipv, recpvvec, rblemcpv);
  };
  
  // find nr of false PV
  int nFalsePV = 0;
  int nFalsePV_real = 0;
  for(int ipv = 0; ipv < (int) recpvvec.size(); ipv++) {
    int fake = 0;
    double x = recpvvec[ipv].x;
    double y = recpvvec[ipv].y;
    double z = recpvvec[ipv].z;
    double r = std::sqrt(x*x + y*y);
    double errx = recpvvec[ipv].positionSigma.x;
    double erry = recpvvec[ipv].positionSigma.y;
    double errz = recpvvec[ipv].positionSigma.z;
    double errr = std::sqrt(((x*errx)*(x*errx)+(y*erry)*(y*erry))/(x*x+y*y));
    double minRDTrack = recpvvec[ipv].minTrackRD;
    double maxRDTrack = recpvvec[ipv].maxTrackRD;
    int mother = recpvvec[ipv].mother;
    double velo = recpvvec[ipv].nVeloTracks;
    double lg = recpvvec[ipv].nLongTracks;
    double d0 = recpvvec[ipv].d0;
    double d0nTr = recpvvec[ipv].d0nTr;
    double chi2nTr = recpvvec[ipv].chi2nTr;
    double mind0 = recpvvec[ipv].mind0;
    double maxd0 = recpvvec[ipv].maxd0;
    double chi2 = recpvvec[ipv].chi2;
    double nDoF = recpvvec[ipv].nDoF;


    if (recpvvec[ipv].indexMCPVInfo < 0 ) {
      nFalsePV++; 
      fake = 1;
      bool vis_found = false;
      for(unsigned int imc = 0; imc < not_rble_but_visible.size() ; imc++) {
        if ( not_rble_but_visible[imc].indexRecPVInfo > -1 ) continue;
        double dist = fabs(mcpvvec[imc].pMCPV->z - recpvvec[ipv].z);
        if (  dist < 5.0 * recpvvec[ipv].positionSigma.z ) {
          vis_found = true;
          not_rble_but_visible[imc].indexRecPVInfo = 10;
          break;
        }
      } // imc
      if ( !vis_found ) nFalsePV_real++;
    }
  }

  // Fill distance to closest recble MC PV and its multiplicity
  std::vector<MCPVInfo>::iterator itmcl;
  for(itmcl = rblemcpv.begin(); rblemcpv.end() != itmcl; itmcl++) {
    std::vector<MCPVInfo>::iterator cmc = closestMCPV(rblemcpv,itmcl);
    double dist = 999999.;
    int mult = 0;
    if(cmc != rblemcpv.end()) {
      
      double diff_x = cmc->pMCPV->x - itmcl->pMCPV->x;
      double diff_y = cmc->pMCPV->y - itmcl->pMCPV->y;
      double diff_z = cmc->pMCPV->z - itmcl->pMCPV->z;
      double dist = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);
      mult = cmc->nRecTracks;
    }
    itmcl->distToClosestMCPV = dist;
    itmcl->multClosestMCPV = mult;
  }

   // Counters
  int nMCPV                 = rblemcpv.size()-nmrc;
  int nRecMCPV              = 0;
  int nMCPV_isol            = 0;
  int nRecMCPV_isol         = 0;
  int nMCPV_close           = 0;
  int nRecMCPV_close        = 0;



  for(itmc = rblemcpv.begin(); rblemcpv.end() != itmc; itmc++) {
    if(itmc->distToClosestMCPV > m_dzIsolated) nMCPV_isol++;
    if(itmc->distToClosestMCPV < m_dzIsolated) nMCPV_close++;
    if(itmc->indexRecPVInfo > -1) {
      nRecMCPV++;
      if(itmc->distToClosestMCPV > m_dzIsolated) nRecMCPV_isol++;
      if(itmc->distToClosestMCPV < m_dzIsolated) nRecMCPV_close++;
    }
  }

  m_nMCPV                 +=  nMCPV;
  m_nRecMCPV              +=  nRecMCPV;
  m_nMCPV_isol            +=  nMCPV_isol;
  m_nRecMCPV_isol         +=  nRecMCPV_isol;
  m_nMCPV_close           +=  nMCPV_close;
  m_nRecMCPV_close        +=  nRecMCPV_close;
  m_nFalsePV              +=  nFalsePV;
  m_nFalsePV_real         +=  nFalsePV_real;


  //loop over matched MC PVs and get pull and errors
  for(auto mc_vertex_info : rblemcpv) {
    int rec_index = mc_vertex_info.indexRecPVInfo;
    if(rec_index < 0) continue;
    MCVertex* mc_vertex = mc_vertex_info.pMCPV; 
    double diff_x = recpvvec[rec_index].x - mc_vertex->x;
    double diff_y = recpvvec[rec_index].y - mc_vertex->y;
    double diff_z = recpvvec[rec_index].z - mc_vertex->z;
    std::cout << "true PV z " << mc_vertex->z << " " <<  mc_vertex->numberTracks <<  std::endl;
    vec_diff_x.push_back(diff_x);
    vec_diff_y.push_back(diff_y);
    vec_diff_z.push_back(diff_z);

    double err_x = recpvvec[rec_index].positionSigma.x;
    double err_y = recpvvec[rec_index].positionSigma.y;
    double err_z = recpvvec[rec_index].positionSigma.z;

    vec_err_x.push_back(err_x);
    vec_err_y.push_back(err_y);
    vec_err_z.push_back(err_z);
    
  }
  std::cout << "-----" << std::endl;

  //loop over reconstructed vertices and print associated tracks
  std::cout << "reconstructed vertices:" << std::endl;
  for(auto p_rec_vertex : vecOfVertices ) {
    auto rec_vertex = * p_rec_vertex;
    std::cout << "x: " << rec_vertex.vb.x << std::endl;
    std::cout << "y: " << rec_vertex.vb.y << std::endl;
    std::cout << "z: " << rec_vertex.vb.z << std::endl;
    //for(auto track : rec_vertex.tracks) std::cout << track.x << std::endl;

  std::cout << "-----" << std::endl;
  }

  std::cout.precision(4);
  std::cout << " ============================================" << std::endl;
  std::cout << " Efficiencies for reconstructible MC vertices: " << std::endl;
  std::cout << " ============================================" << std::endl;
  std::cout << " " << std::endl;

  std::cout << " MC PV is reconstructible if at least " << m_nTracksToBeRecble
         << "  tracks are reconstructed" << std::endl;
  std::cout << " MC PV is isolated if dz to closest reconstructible MC PV >  "
         << m_dzIsolated << " mm" << std::endl;
  std::string ff = "by counting tracks";
  if ( !m_matchByTracks ) ff = "by dz distance";
  std::cout << " REC and MC vertices matched:  "
         <<  ff << std::endl;


  std::cout << " " << std::endl;


            printRat("All",       m_nRecMCPV ,       m_nMCPV );
            printRat("Isolated",  m_nRecMCPV_isol,   m_nMCPV_isol );
            printRat("Close",     m_nRecMCPV_close,  m_nMCPV_close );
            printRat("False rate",m_nFalsePV ,       m_nRecMCPV+m_nFalsePV );


          printRat("Real false rate", m_nFalsePV_real , m_nRecMCPV+m_nFalsePV_real);

   std::cout <<   "new found: " <<     m_nRecMCPV << " / " << m_nMCPV  << std::endl;
   std::cout << "new fakes: " << m_nFalsePV << std::endl;


}


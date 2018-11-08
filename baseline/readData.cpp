#include "readData.h"

VeloState  make_new_velostate()  {return VeloState()   ; }
XYZPoint   make_new_xyzpoint()   {return XYZPoint()    ; }
VertexBase make_new_vertexbase() {return VertexBase()  ; }
MCVertex   make_new_mcvertex()   {return MCVertex()    ; }

uint32_t get_max_tracks(){return PatPV::max_tracks;}
uint32_t get_max_seeds(){return PatPV::max_seeds;}
uint32_t get_max_vertices(){return PatPV::max_vertices;}

uint32_t readTracks(VeloState host_velo_states[], 
                    char* filename) {

  //now loop over file and convert json to arrays read by PatPV algorithm
  //keep track of current event number
  boost::property_tree::ptree tree;

  // Parse the XML into the property tree.
  boost::property_tree::read_json(filename, tree);

  auto velotrack_tree    = tree.get_child("VeloTracks"); 

  //keep track of current track number
  uint32_t current_track_number = 0;
  for(auto velotrack : velotrack_tree) {

    auto sub_tree_closest_to_beam = velotrack.second.get_child("ClosestToBeam");
    std::vector<double> track_paras;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, velotrack.second.get_child("ClosestToBeam")) {
      
      track_paras.push_back(v.second.get_value<double>()) ;
    }
     
    std::vector<double> cov_matrix;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, velotrack.second.get_child("errCTBState")) {
 
      cov_matrix.push_back(v.second.get_value<double>()) ;
    }

    VeloState velo_state;
    velo_state.x = track_paras.at(0);
    velo_state.y = track_paras.at(1);
    velo_state.z = track_paras.at(2);
    velo_state.tx = track_paras.at(3);
    velo_state.ty = track_paras.at(4);

    velo_state.c00 = cov_matrix.at(0);
    velo_state.c20 = cov_matrix.at(4);
    velo_state.c22 = cov_matrix.at(2);
    velo_state.c11 = cov_matrix.at(1);
    velo_state.c31 = cov_matrix.at(4);
    velo_state.c33 = cov_matrix.at(3);
    //fill velo state array
    host_velo_states[current_track_number] = velo_state;

    current_track_number++;
  } //loop over velo tracks

  return current_track_number;

}

uint32_t readMCVertices(MCVertex mcvertices_all[],
                        char* filename) {

  //now loop over file and convert json to arrays read by PatPV algorithm
  //keep track of current event number
  uint32_t accumulated_vertex_number = 0;

  boost::property_tree::ptree tree;

  // Parse the XML into the property tree.
  boost::property_tree::read_json(filename, tree);

  //get MC vertices
  auto mcvertices_tree    = tree.get_child("MCVertices"); 
  int counter_mcvertices = 0;
  for(auto mcvertex : mcvertices_tree) {

    std::vector<double> mcvertex_pos;
   
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, mcvertex.second.get_child("Pos")) {
      //std::cout << v.second.get_value<double>() << std::endl;
      mcvertex_pos.push_back(v.second.get_value<double>()) ;
     }  
     MCVertex mc_vertex;
     mc_vertex.x = mcvertex_pos.at(0);
     mc_vertex.y = mcvertex_pos.at(1);
     mc_vertex.z = mcvertex_pos.at(2);
     mc_vertex.numberTracks = mcvertex.second.get<int>("products");
     mcvertices_all[accumulated_vertex_number] = mc_vertex;
     ++accumulated_vertex_number;
  }

  return accumulated_vertex_number;

}

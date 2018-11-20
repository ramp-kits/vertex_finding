#ifndef RAMP_READDATA_H
#define RAMP_READDATA_H

#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <sstream>
#include <iostream>
#include <string>
#include <typeinfo>
#include "patPV_Definitions.cuh"
#include <dirent.h>
#include "json.hpp"

extern "C" {

uint32_t readTracks(VeloState VeloState[], 
                    char* filename);

uint32_t readMCVertices(MCVertex mcvertices_all[],
                        char* filename);

VeloState   make_new_velostate();
XYZPoint    make_new_xyzpoint();
VertexBase  make_new_vertexbase();
MCVertex    make_new_mcvertex();

uint32_t get_max_tracks();
uint32_t get_max_seeds();
uint32_t get_max_vertices();

}

#endif //

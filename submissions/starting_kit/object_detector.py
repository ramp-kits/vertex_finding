import numpy as np


import os,sys
import ctypes
from ctypes import *
from array import array

#does this work on different architectures?
lib = cdll.LoadLibrary('submissions/starting_kit/libRAMP.so')

class VeloState(Structure):
    _fields_=[("x",c_float),
              ("y",c_float),
              ("z",c_float),
              ("tx",c_float),
              ("ty",c_float),
              ("c00",c_float),
              ("c20",c_float),
              ("c22",c_float),
              ("c11",c_float),
              ("c31",c_float),
              ("c33",c_float),
              ("chi2",c_float),
              ("backward",c_bool)]

class XYZPoint(Structure):
    _fields_=[("x",c_float),
              ("y",c_float),
              ("z",c_float)]

class VertexBase(Structure):
    _fields_=[("x",c_double),
              ("y",c_double),
              ("z",c_double),
              ("c00",c_double),
              ("c10",c_double),
              ("c11",c_double),
              ("c20",c_double),
              ("c21",c_double),
              ("c22",c_double),
              ("chi2",c_double),
              ("ndof",c_int)]

class MCVertex(Structure):
    _fields_=[("x",c_float),
              ("y",c_float),
              ("z",c_float),
              ("numberTracks",c_int)]


# Defining pointer to the make_new_velo_state() function
make_new_velostate = lib.make_new_velostate

# Declaring the function return type - a pointer to a VeloState object - just like in the C code
make_new_velostate.restype = VeloState

# Declaring list of parameter types. In this case, the list contains zero items,
# as the function has only one parameter
make_new_velostate.argtypes = []




class ObjectDetector:
    def __init__(self):
        pass

    def fit(self, X, y):
        return self

    def predict(self, X):
        #loop over all events saved in X
        #list of of lists of all reconstructed vertices
        all_rec_vertex_list = []
        #loop over events
        for event in X:
          velo_states_in_event = []
          #loop over velo states in event
          for velo_state in event.tracks:
            #velo_state_cov = velo_state[1]
            #velo_state = velo_state[0]
            full_velo_state = make_new_velostate()
            full_velo_state.x = velo_state.x
            full_velo_state.y = velo_state.y
            full_velo_state.z = velo_state.z

            full_velo_state.tx = velo_state.tx
            full_velo_state.ty = velo_state.ty

            full_velo_state.c00 = velo_state.cov_x
            full_velo_state.c20 = velo_state.cov_xtx
            full_velo_state.c22 = velo_state.cov_tx
            full_velo_state.c11 = velo_state.cov_y
            full_velo_state.c31 = velo_state.cov_xtx
            full_velo_state.c33 = velo_state.cov_ty
            velo_states_in_event = velo_states_in_event + [full_velo_state]


                    # Get the sizes of the arrays to make
          max_tracks   = lib.get_max_tracks()
          max_seeds    = lib.get_max_seeds()
          max_vertices = lib.get_max_vertices()
          # Arrays are made available as multiplications
          # of the class by a certain number of elements
          host_velo_states = max_tracks*VeloState
          tracks2disable = max_tracks*c_bool
          seeds = max_seeds*XYZPoint
          outvtxvec = max_vertices*VertexBase
          mcvertices_all = max_vertices*MCVertex

          number_of_tracks = ctypes.c_uint(0)
          num_mc_vertices  = ctypes.c_uint(0)
          num_rec_vertices = ctypes.c_uint(0)

          # There must be a neater way to do this but I haven't found it
          # basically the above multiplication defines a construction for
          # the array, then you must instantiate it
          a  = host_velo_states()
          for i in range(0, len(velo_states_in_event)):
            a[i] = velo_states_in_event[i]
          
          

          read_mcvertices = mcvertices_all()
          read_velo_states = host_velo_states()
          read_tracks2disable = tracks2disable()
          recod_outvtxvec = outvtxvec()
          recod_seeds = seeds()


          #number_vertex = lib.reconstructMultiPVFromTracks(a, recod_outvtxvec, 
           #                                         read_tracks2disable, recod_seeds, len(velo_states_in_event))
          

          number_vertex = lib.reconstructMultiPVFromTracks(a, recod_outvtxvec, 
                                                    read_tracks2disable, recod_seeds, len(velo_states_in_event))
          
          rec_vertex_in_event = []

          for i in range (0, number_vertex):
            vertex_tuple = ( recod_outvtxvec[i].x, recod_outvtxvec[i].y, recod_outvtxvec[i].z)
          #  vertex_tuple = (1., recod_outvtxvec[i])
            rec_vertex_in_event = rec_vertex_in_event + [vertex_tuple]
          all_rec_vertex_list = all_rec_vertex_list + [rec_vertex_in_event]




        #lib.reconstructMultiPVFromTracks(read_velo_states, recod_outvtxvec, 
         #                                           read_tracks2disable, recod_seeds, number_of_tracks)
        y_pred =  all_rec_vertex_list
        # convert output into an np.array of objects
        y_pred_array = np.empty(len(y_pred), dtype=object)
        y_pred_array[:] = y_pred
        return y_pred_array

import numpy as np


import os,sys
import ctypes
from ctypes import *
from array import array

#does this work on different architectures?
import glob





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
          #random number of vertices according to Poisson distribution
          number_vertex = np.random.poisson(6)
          rec_vertex_in_event = []
          for i in range (0, number_vertex):
            x = np.random.normal(0., 0.05)
            y = np.random.normal(0., 0.05)
            z = np.random.normal(0., 50.)
            vertex_tuple = (x, y, z)
            rec_vertex_in_event = rec_vertex_in_event + [vertex_tuple]
          all_rec_vertex_list = all_rec_vertex_list + [rec_vertex_in_event]


        y_pred =  all_rec_vertex_list
        # convert output into an np.array of objects
        y_pred_array = np.empty(len(y_pred), dtype=object)
        y_pred_array[:] = y_pred
        return y_pred_array

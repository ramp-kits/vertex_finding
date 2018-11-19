from __future__ import division, print_function

import json
import os

import numpy as np
import pandas as pd

import rampwf as rw
from rampwf.prediction_types.base import BasePrediction
from rampwf.score_types import BaseScoreType

# ramp-kit for the RAPID challenge test and train data should be in data/test
# and data/train one has to implement the ObjectDetector class in
# submissions/<yoursubmission>/object_detector.py with the functions fit(X,y)
# and predict(X), where X is the input data and y the truth information the fit
# function contains the training of the model, while the predict function
# applies it


# to do: implement various scoring algorithms
# to do: decide on crossfolds, prediction

def z_dist_matched(rec_pv_z, mc_pv_z, m_distance):
    return abs(rec_pv_z - mc_pv_z) < m_distance

# define matching criteria, based on distance:


def is_matched(rec_pv_x=0., rec_pv_y=0., rec_pv_z=0., mc_pv_x=0.,
               mc_pv_y=0., mc_pv_z=0., m_distance=0.1):
    return z_dist_matched(rec_pv_z, mc_pv_z, m_distance)


# get the number of Velo tracks of the PV with key pv_key
def get_number_velo_tracks(pv_key, jdata):
    pv_key = int(pv_key)
    MCParticles = jdata['MCParticles']
    VeloTracks = jdata['VeloTracks']
    counter_velo_tracks = 0
    for mcp_key, mcp in MCParticles.items():
        is_pv = bool(mcp['ispv'])
        is_velo = bool(mcp['isvelo'])
        if pv_key == int(mcp['PV']) and is_pv and is_velo:
            counter_velo_tracks += sum([
                1 for v_track in VeloTracks.values()
                if mcp_key in v_track['MCPs']])

    return counter_velo_tracks

# class to do checking and plots


class PVChecker:
    def __init__(self):
                # print("PVChecker: checking efficiency")
                # configuration for matching
        self.m_mintracks = 2
        self.m_distance = 0.3

        self.counter_found_MC_PV = 0.
        self.counter_total_MC_PV = 0.
        self.counter_total_MC_PV_reconstructible = 0.
        self.counter_fake_PV = 0.

        self.res_x = np.empty(0)
        self.res_y = np.empty(0)
        self.res_z = np.empty(0)

    # load data
    # here we expect arr_rec_pvs to be numpy array of array[x,y,z] and
    # arr_mc_pvs to be numpy array of array[x,y,z, nTracks]

    def load_data_rec(self, arr_rec_pvs):
        self.df_rec_pvs = pd.DataFrame(arr_rec_pvs)
        self.df_rec_pvs['matched'] = 0
        # entry number of matched MC PV, -99 if not matched
        self.df_rec_pvs['matched_pv_key'] = -99
        self.df_rec_pvs.columns = ['x', 'y', 'z', 'matched', 'matched_pv_key']

    def load_data_true(self, arr_mc_pvs):
        self.df_mc_pvs = pd.DataFrame(arr_mc_pvs)
        self.df_mc_pvs.columns = ['x', 'y', 'z', 'nVeloTracks']

    def load_from_ramp(self, y_true_label_index, y_pred_label_index):

        # for MC_PVs in y_true_label_index:
        # loop over event
        for i_event in range(0, len(y_true_label_index)):
            MCPV_arr_tot = []
            RecPV_arr_tot = []
            # set-up MC PVs
            for MC_PV in y_true_label_index[i_event]:
                MCPV_arr = np.array(
                    [MC_PV.x, MC_PV.y, MC_PV.z, MC_PV.numberTracks])
                MCPV_arr_tot = MCPV_arr_tot + [MCPV_arr]
            # set-up reconstructed PVs
            for Rec_PV in y_pred_label_index[i_event]:
                RecPV_arr = np.array(Rec_PV)
                RecPV_arr_tot = RecPV_arr_tot + [RecPV_arr]
            MCPV_arr_tot = np.array(MCPV_arr_tot)
            RecPV_arr_tot = np.array(RecPV_arr_tot)
            self.load_data_true(MCPV_arr_tot)
            # catch cases where we reconstructed no PVs
            if RecPV_arr_tot.size != 0:
                self.load_data_rec(RecPV_arr_tot)
                self.check_event_df()

        # check event with previously loaded data frames

    def check_event_df(self):
        # loop over MC PVs and find rec PV with minimum z distance
        for mc_index, mc_pv in self.df_mc_pvs.iterrows():
            if mc_pv['nVeloTracks'] < self.m_mintracks:
                continue

            # loop over rec PVs
            true_z = mc_pv['z']
            min_dist = 10000.
            index_min_dist = -99
            for rec_index, rec_pv in self.df_rec_pvs.iterrows():
                rec_z = rec_pv['z']
                dist_z = abs(true_z - rec_z)
                if dist_z < min_dist:
                    min_dist = dist_z
                    index_min_dist = rec_index
            rec_z = self.df_rec_pvs['z'][index_min_dist]
            dist_z = abs(true_z - rec_z)
            # match rec and MC PVs, if the rec pv with minimum z distance
            # to MC PV fullfills matching crtierion
            crit1 = is_matched(rec_pv_z=rec_z,
                               mc_pv_z=true_z,
                               m_distance=self.m_distance)
            crit2 = self.df_rec_pvs['matched'][index_min_dist]
            if crit1 and not crit2:
                self.df_rec_pvs.loc[index_min_dist, 'matched'] = 1
                self.df_rec_pvs.loc[index_min_dist,
                                    'matched_pv_key'] = mc_index

        # test creating sub dataframes of real and fake rec pv
        df = self.df_rec_pvs[self.df_rec_pvs.matched == 1]
        self.df_fake_rec_pvs = self.df_rec_pvs[self.df_rec_pvs.matched == 0]

        df_true = pd.DataFrame(
            columns=['true_x', 'true_y', 'true_z'], dtype=float)
        for key, row in df.iterrows():
            df_true.loc[key, ['true_x', 'true_y', 'true_z']] = (
                self.df_mc_pvs.loc[
                    row['matched_pv_key'], ['x', 'y', 'z']].values
            )
        df = pd.concat([df, df_true], axis=1)
        for dim in ['x', 'y', 'z']:
            df['residual_' + dim] = df[dim] - df['true_' + dim]
        self.res_x = np.append(self.res_x, df['x'] - df['true_x'])
        self.res_y = np.append(self.res_y, df['y'] - df['true_y'])
        self.res_z = np.append(self.res_z, df['z'] - df['true_z'])
        self.df_true_rec_pvs = df

        self.counter_found_MC_PV = (
            self.counter_found_MC_PV + self.df_true_rec_pvs.index.size
        )
        self.counter_total_MC_PV = (
            self.counter_total_MC_PV + self.df_mc_pvs.index.size
        )
        self.counter_total_MC_PV_reconstructible = (
            self.counter_total_MC_PV_reconstructible + self.df_mc_pvs[
                self.df_mc_pvs.nVeloTracks >= self.m_mintracks
            ].index.size
        )
        self.counter_fake_PV = (
            self.counter_fake_PV + self.df_fake_rec_pvs.index.size
        )

    def calculate_eff(self):

        self.total_efficiency = (
            self.counter_found_MC_PV / self.counter_total_MC_PV
        )
        self.total_fake_rate = (
            self.counter_fake_PV / (
                self.counter_found_MC_PV + self.counter_fake_PV
            )
        )
        self.reconstructible_efficiency = (
            self.counter_found_MC_PV / (
                self.counter_total_MC_PV_reconstructible
            )
        )

    def calculate_resolution(self):
        self.sigma_x = self.res_x.std()
        self.sigma_y = self.res_y.std()
        self.sigma_z = self.res_z.std()

    def effective_efficiency(self):
        self.effective_eff = (
            self.reconstructible_efficiency * (
                1. - 2. * self.total_fake_rate) ** 2
        )
    # function to get determine total score

    def final_score(self):
        # critertia: efficiency, fake rate, sigma of residuals,
        # means of residuals?
        self.fin_score = self.reconstructible_efficiency * \
            (1. - 2. * self.total_fake_rate) ** 2 / \
            self.sigma_x / self.sigma_y / self.sigma_z / 100000.

        # print("the final score is", self.fin_score, "!")


class PVScore_total(BaseScoreType):
    is_lower_the_better = False
    minimum = 0.0
    maximum = 1.0

    def __init__(self, mode, name='total score', precision=2):
        self.name = name
        self.precision = precision
        self.mode = mode

    def __call__(self, y_true_label_index, y_pred_label_index):
        # we can us the python PVChecker -> need to transform data for it
        checker = PVChecker()
        checker.load_from_ramp(y_true_label_index, y_pred_label_index)
        checker.calculate_eff()
        checker.calculate_resolution()
        checker.final_score()
        checker.effective_efficiency()
        if self.mode == "total":
            return checker.fin_score
        if self.mode == "eff":
            return checker.reconstructible_efficiency
        if self.mode == "fake":
            return checker.total_fake_rate
        if self.mode == "effective_eff":
            return checker.effective_eff
        if self.mode == "resolution":
            return checker.sigma_z

    def check_y_pred_dimensions(self, y_true, y_pred):
        if len(y_true) != len(y_pred):
            msg = ('Wrong y_pred dimensions: y_pred should have {} '
                   'instances, instead it has {} instances')
            raise ValueError(msg.format(len(y_true), len(y_pred)))


class PVPredictions(BasePrediction):
    def __init__(self, y_pred=None, y_true=None, n_samples=None):
        if y_pred is not None:
            self.y_pred = y_pred
        elif y_true is not None:
            self.y_pred = y_true
        elif n_samples is not None:
            self.y_pred = np.empty(n_samples, dtype=object)
        else:
            raise ValueError(
                'Missing init argument: y_pred, y_true, or n_samples')

    def __str__(self):
        return 'y_pred = {}'.format(self.y_pred)

    @classmethod
    # combination at the moment dummy implementation
    def combine(cls, predictions_list, index_list=None):
        # if index_list is None:  # we combine the full list
        #     index_list = range(len(predictions_list))
        # y_comb_list = [predictions_list[i].y_pred for i in index_list]

        # n_preds = len(y_comb_list)
        # y_preds_combined = np.empty(n_preds, dtype=object)
        # combined_predictions = cls(y_pred=y_preds_combined)
        combined_predictions = cls(y_pred=predictions_list[0].y_pred)
        return combined_predictions

    @property
    def valid_indexes(self):
        return self.y_pred != np.empty(len(self.y_pred), dtype=np.object)
        # return True


problem_title = 'RAPID challenge'
# A type (class) which will be used to create wrapper objects for y_pred

Predictions = PVPredictions
# An object implementing the workflow
workflow = rw.workflows.ObjectDetector()

# The overlap between adjacent patches is 56 pixels
# The scoring region is chosen so that despite the overlap,
# no crater is scored twice, hence the boundaries of
# 28 = 56 / 2 and 196 = 224 - 56 / 2


score_types = [
    PVScore_total(name="total", mode="total"),
    PVScore_total(name="efficiency", mode="eff"),
    PVScore_total(name="fake rate", mode="fake"),
    PVScore_total(name="z resolution", mode="resolution"),
    PVScore_total(name="effective efficiency", mode="effective_eff")
]


def get_cv(X, y):
        # 3 quadrangles for training have not exactly the same size,
        # but for simplicity just cut in 3
        # for each fold use one quadrangle as test set,
        # the other two as training
    n_tot = len(X)
    # first entry in tuple is for training, the second for testing
    # number of tuples gives number of crossfolds
    return [
        (np.r_[0:n_tot], np.r_[0:n_tot])]


class MCVertex:
    def __init__(self, x, y, z, numberTracks):
        self.x = x
        self.y = y
        self.z = z
        self.numberTracks = numberTracks

    def __repr__(self):
            # return '{0}, {1}, {2}'.format(self.x, self.y, self.z)
        return 'MCVertex'

    def __str__(self):
        return 'MCVertex'


# helper calss to hold Velo state + covariance matrix
class VeloState:
    def __init__(self, x, y, z, tx, ty, pq,
                 cov_x, cov_y, cov_tx, cov_ty, cov_xtx):
        self.x = x
        self.y = y
        self.z = z
        self.tx = tx
        self.ty = ty
        self.pq = pq
        self.cov_x = cov_x
        self.cov_y = cov_y
        self.cov_tx = cov_tx
        self.cov_ty = cov_ty
        self.cov_xtx = cov_xtx

    def __repr__(self):
        return ('VeloState({0},{1},{2},{3},{4})'
                .format(self.x, self.y, self.z, self.tx, self.ty))

    def __str__(self):
        return ('VeloState({0},{1},{2},{3},{4})'
                .format(self.x, self.y, self.z, self.tx, self.ty))

# helper class to hold Velo hits


class VeloHit:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return 'VeloHit({0}, {1}, {2})'.format(self.x, self.y, self.z)

    def __str__(self):
        return 'VeloHit({0}, {1}, {2})'.format(self.x, self.y, self.z)


# class to hold tracks and hits of an event
class EventData:
    def __init__(self, list_tracks, list_hits):
        self.tracks = list_tracks
        self.hits = list_hits

    def __repr__(self):
        return 'EventData'

    def __str__(self):
        return 'EventData'


def _read_data(path, type):
    """
    Read and process data and labels.

    Parameters
    ----------
    path : path to directory that has 'data' subdir
    typ : {'train', 'test'}

    Returns
    -------
    X, y data
    X: np array of EventData, where EventData conists of list of VeloStates
       and list of VeloHits for a event

    Y: np array of lists of MCVertex, where a list contains all MC vertices of
       an event, and the array contains the lists of all events
    """
    # loop over all json files:

    list_y = []
    list_x = []
    # default path is .
    # have to set it for reading
    test = os.getenv('RAMP_TEST_MODE', 0)
    path = path + '/data/{0}/'.format(type)
    list_of_files = os.listdir(path)
    # for testing, run over subset of files
    if test:
        list_of_files = list_of_files[:100]
    for file in list_of_files:
        if not file.endswith('.json'):
            continue
        file_path = path + file
        jdata = json.load(open(file_path))
        MCVertices = jdata['MCVertices']

        mc_pvs = [
            MCVertex(h['Pos'][0], h['Pos'][1], h['Pos'][2],
                     get_number_velo_tracks(key, jdata))
            for key, h in MCVertices.items()
        ]
        list_y = list_y + [mc_pvs]

        VeloTracks = jdata['VeloTracks']
        VeloHits = jdata['VPClusters']
        # velo_states = [tuple(h['ClosestToBeam']) for key,h in
        # VeloTracks.items() ]
        velo_states = [
            VeloState(h['ClosestToBeam'][0],
                      h['ClosestToBeam'][1],
                      h['ClosestToBeam'][2],
                      h['ClosestToBeam'][3],
                      h['ClosestToBeam'][4],
                      h['ClosestToBeam'][5],
                      h['errCTBState'][0],
                      h['errCTBState'][1],
                      h['errCTBState'][2],
                      h['errCTBState'][3],
                      h['errCTBState'][4])
            for key, h in VeloTracks.items()]
        velo_hits = [VeloHit(h['x'], h['y'], h['z'])
                     for key, h in VeloHits.items()]
        event = EventData(velo_states, velo_hits)

        list_x = list_x + [event]

    y_array = np.empty(len(list_y), dtype=object)
    y_array[:] = list_y

    x_array = np.empty(len(list_x), dtype=object)
    x_array[:] = list_x
    x_array = np.array(x_array)

    return x_array, y_array


def get_test_data(path):
    return _read_data(path, 'test')
    # return np.array([1]),np.array([2])


def get_train_data(path):
    return _read_data(path, 'train')

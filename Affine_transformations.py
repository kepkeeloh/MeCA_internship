import numpy as np
import pandas as pd
import sys
from read_file import *
import os


def SquareToSphere(dimRectP1, dimRectP2, sulciP1, sulciP2, latP1, latP2):
    """
    As the system of coordinates is not the same in both the rectangle and the sphere for the models we are considering,
    we need to compute the position of the axes of each species (Primate1 and Primate2) on the sphere (i.e. the cortical
    surface), given the respective dimensions of the rectangles of each species and the wanted dimensions of the sphere
    (which might be different from ([0,360]*[0,180]) as we do not take into account the intervals where the poles are
    on the latitudes: with the default value, below 30 and above 150 are respectively the insular and cingular poles).
    :param latP1: [m,n] such that the m first degrees above 0 of latitudes correspond to the insular pole
    and the n first degrees from 180 downwards correspond to the cingular pole, for Primate1
    :param latP2: same thing as latP1 for Primate2
    :param dimRectP1: the dimensions of the rectangle for the Primate1 as [Longitude,latitude]
    :param dimRectP2: same thing as dimP1 for Primate2
    :param sulciP1: a tuple of two lists of floats, respectively, the coordinates of the longitudinal and the
     latitudinal sulci for Primate1 on the rectangle
    :param sulciP2: same thing as sulciP2 for Primate2
    :return: four lists of floats, respectively, the coordinates of the longitudinal and the latitudinal
    sulci for Primate1 and the coordinates of the longitudinal and the latitudinal sulci for Primate2, on the sphere.
    """

    # extract cardinals
    Nlon_P1 = len(sulciP1[0])
    Nlat_P1 = len(sulciP1[1])
    Nlon_P2 = len(sulciP2[0])
    Nlat_P2 = len(sulciP2[1])

    # extract rectangle dimensions and latitudes' boundaries
    LP1, lP1 = dimRectP1
    LP2, lP2 = dimRectP2
    insP1, cinP1 = latP1[0], 180 - latP1[1]
    insP2, cinP2 = latP2[0], 180 - latP2[1]

    # initialise lists
    newLon_P1 = np.copy(sulciP1[0])
    newLat_P1 = np.copy(sulciP1[1])
    newLon_P2 = np.copy(sulciP2[0])
    newLat_P2 = np.copy(sulciP2[1])

    for i in range(Nlon_P1):
        newLon_P1[i] *= 360 / LP1
        newLon_P1[i] += (sulciP1[0][i] < 0) * 360

    for i in range(Nlat_P1):
        if insP1 < sulciP1[1][i] < cinP1:
            newLat_P1[i] *= (cinP1 - insP1) / lP1
            newLat_P1[i] += insP1

    for i in range(Nlon_P2):
        newLon_P2[i] *= 360 / LP2
        newLon_P2[i] += (sulciP2[0][i] < 0) * 360

    for i in range(Nlat_P2):
        if insP2 < sulciP2[1][i] < cinP2:
            newLat_P2[i] *= (cinP2 - insP2) / lP2
            newLat_P2[i] += insP2

    return [newLon_P1, newLat_P1], [newLon_P2, newLat_P2]


def Affine_Transform(sulciP1, sulciP2, long_corr, lat_corr, latP1, latP2):
    """
    Given the coordinates of the sulcal lines on the sphere (i.e. cortical surface) for each species (Primate1 and
     Primate2), computes and returns the affine transformations on each interval between the corresponding axes
     for longitudinal and latitudinal sulci from Primate1 to Primate2 (i.e. meant to modify the Primate2's coordinate's
     system)
    :param latP1: [m,n] such that the m first degrees above 0 of latitudes correspond to the insular pole
    and the n first degrees from 180 downwards correspond to the cingular pole, for Primate1
    :param latP2: same thing as latP1 for Primate2
    :param sulciP1: a list of two lists of floats, respectively, the coordinates of the longitudinal and the
     latitudinal sulci for Primate1 on the sphere
    :param sulciP2: same thing as sulciP1 for Primate2
    :param long_corr: a numpy array of shape (N,2) of integers giving the indices of the N corresponding sulci on the
    longitudes, with the first index corresponding to Primate1 and the second for Primate2
    :param lat_corr: same thing as long_corr for the latitudinal sulci
    :return: a list of the longitudinal coordinates that define the intervals on which we computed the affine
    transformations, including the boundaries, a list of couples of floats (a,b) that correspond to the affine
    transformations y = ax + b on the longitudinal intervals, and respectively the same lists for the latitudes.
    """

    # extract cardinals
    Nlong = len(long_corr)
    Nlat = len(lat_corr)

    # extract latitudes' boundaries
    insP1, cinP1 = latP1[0], 180 - latP1[1]
    insP2, cinP2 = latP2[0], 180 - latP2[1]

    # initialization of the lists
    long_transform = np.zeros((Nlong + 1, 2))
    lat_transform = np.zeros((Nlat + 1, 2))

    # make the lists of the axes that have a correspondence (and therefore define the intervals)
    longP1 = np.sort(np.concatenate(([0], sulciP1[0][long_corr[:, 0]], [360])))
    latP1 = np.sort(np.concatenate(([insP1], sulciP1[1][lat_corr[:, 0]], [cinP1])))
    longP2 = np.sort(np.concatenate(([0], sulciP2[0][long_corr[:, 1]], [360])))
    latP2 = np.sort(np.concatenate(([insP2], sulciP2[1][lat_corr[:, 1]], [cinP2])))

    for i in range(Nlong + 1):
        long_transform[i][0] = (longP1[i + 1] - longP1[i]) / (longP2[i + 1] - longP2[i])
        long_transform[i][1] = longP1[i] - longP2[i] * long_transform[i][0]

    for i in range(Nlat + 1):
        lat_transform[i][0] = (latP1[i + 1] - latP1[i]) / (latP2[i + 1] - latP2[i])
        lat_transform[i][1] = latP1[i] - latP2[i] * lat_transform[i][0]

    return longP2, long_transform, latP2, lat_transform


####################################################################
#
# main function
#
# python Affine_transformations.py Primate1 Primate2 side
#
####################################################################

def main(Primate1, Primate2, side):
    dir = Primate1 + '_' + Primate2
    if not os.path.exists(dir):
        os.mkdir(dir)
        print("directory ", dir, " created ")

    print('reading models\' informations')

    modelP1F = 'model_' + Primate1 + '_' + side + '.txt'
    modelP2F = 'model_' + Primate2 + '_' + side + '.txt'
    modelP1 = read_model(modelP1F)
    modelP2 = read_model(modelP2F)
    dimRect_P1, poles_lat_P1, longID_P1, latID_P1, sulci_lon_P1, sulci_lat_P1, lon_coor_P1, lat_coor_P1 = modelP1
    dimRect_P2, poles_lat_P2, longID_P2, latID_P2, sulci_lon_P2, sulci_lat_P2, lon_coor_P2, lat_coor_P2 = modelP2

    print('reading correspondences\' table')

    name_corr = Primate1 + '_' + Primate2 + '_' + 'corr.txt'
    corrTable = read_corr(name_corr)

    print('rescaling square coordinates to sphere coordinates')

    sulciP1, sulciP2 = SquareToSphere(dimRect_P1, dimRect_P2, [lon_coor_P1, lat_coor_P1],
                                      [lon_coor_P2, lat_coor_P2], poles_lat_P1, poles_lat_P2)

    print('extracting correspondences')

    assert (len(corrTable[0]) == len(corrTable[1]) and len(corrTable[2]) == len(corrTable[3])), \
        "Error in the dimensions of the correspondences' table."

    Ncorr_lon = len(corrTable[0])
    Ncorr_lat = len(corrTable[2])
    long_corr = np.zeros((Ncorr_lon, 2)).astype('int')
    lat_corr = np.zeros((Ncorr_lat, 2)).astype('int')
    for i in range(Ncorr_lon):
        long_corr[i][0] = longID_P1[sulci_lon_P1[corrTable[0][i]]]
        long_corr[i][1] = longID_P2[sulci_lon_P2[corrTable[1][i]]]
    for i in range(Ncorr_lat):
        lat_corr[i][0] = latID_P1[sulci_lat_P1[corrTable[2][i]]]
        lat_corr[i][1] = latID_P2[sulci_lat_P2[corrTable[3][i]]]

    print('computing affine transformations from ' + Primate1 + ' to ' + Primate2)
    P1toP2 = Affine_Transform(sulciP1, sulciP2, long_corr, lat_corr, poles_lat_P1, poles_lat_P2)

    print('writing it down')
    nameP1toP2_long_int = Primate1 + '_to_' + Primate2 + '_' + side + 'lon_intervals.csv'
    nameP1toP2_long = Primate1 + '_to_' + Primate2 + '_' + side + 'lon_AffineTrans.csv'
    nameP1toP2_lat_int = Primate1 + '_to_' + Primate2 + '_' + side + 'lat_intervals.csv'
    nameP1toP2_lat = Primate1 + '_to_' + Primate2 + '_' + side + 'lat_AffineTrans.csv'
    pd.DataFrame(P1toP2[0]).to_csv(os.path.join(dir, nameP1toP2_long_int), header=False, index=False)
    pd.DataFrame(P1toP2[1]).to_csv(os.path.join(dir, nameP1toP2_long), header=False, index=False)
    pd.DataFrame(P1toP2[2]).to_csv(os.path.join(dir, nameP1toP2_lat_int), header=False, index=False)
    pd.DataFrame(P1toP2[3]).to_csv(os.path.join(dir, nameP1toP2_lat), header=False, index=False)

    print('computing affine transformations from ' + Primate2 + ' to ' + Primate1)
    # swap the indices for the correspondences' table
    long_corr_inv = np.zeros(np.shape(long_corr), dtype='int')
    long_corr_inv[:, 0] = np.copy(long_corr[:, 1])
    long_corr_inv[:, 1] = np.copy(long_corr[:, 0])
    lat_corr_inv = np.zeros(np.shape(lat_corr), dtype='int')
    lat_corr_inv[:, 0] = np.copy(lat_corr[:, 1])
    lat_corr_inv[:, 1] = np.copy(lat_corr[:, 0])

    P2toP1 = Affine_Transform(sulciP2, sulciP1, long_corr_inv, lat_corr_inv, poles_lat_P2, poles_lat_P1)

    print('writing it down')
    nameP2toP1_long_int = Primate2 + '_to_' + Primate1 + '_' + side + 'lon_intervals.csv'
    nameP2toP1_long = Primate2 + '_to_' + Primate1 + '_' + side + 'lon_AffineTrans.csv'
    nameP2toP1_lat_int = Primate2 + '_to_' + Primate1 + '_' + side + 'lat_intervals.csv'
    nameP2toP1_lat = Primate2 + '_to_' + Primate1 + '_' + side + 'lat_AffineTrans.csv'
    pd.DataFrame(P2toP1[0]).to_csv(os.path.join(dir, nameP2toP1_long_int), header=False, index=False)
    pd.DataFrame(P2toP1[1]).to_csv(os.path.join(dir, nameP2toP1_long), header=False, index=False)
    pd.DataFrame(P2toP1[2]).to_csv(os.path.join(dir, nameP2toP1_lat_int), header=False, index=False)
    pd.DataFrame(P2toP1[3]).to_csv(os.path.join(dir, nameP2toP1_lat), header=False, index=False)

    print('done')

    return None


if __name__ == '__main__':
    Primate1, Primate2, side = sys.argv[1:]
    main(Primate1, Primate2, side)

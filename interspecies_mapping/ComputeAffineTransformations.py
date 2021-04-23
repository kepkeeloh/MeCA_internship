#edited from jeanne's affine_transformation.py script

import numpy as np
import sys
from read_file import *
import os


def SquareToSphere(dimRect, sulci, lat):
    """
    As the system of coordinates is not the same in both the rectangle and the sphere for the models we are considering,
    we need to compute the position of the axes of each species on the sphere (i.e. the cortical
    surface), given the respective dimensions of the rectangle of each species and the wanted dimensions of the sphere
    (which might be different from ([0,360]*[0,180]) as we do not take into account the intervals where the poles are
    on the latitudes: with the default value, below 30 and above 150 are respectively the insular and cingular poles).
    :param lat: [m,n] such that the m first degrees above 0 of latitudes correspond to the insular pole
    and the n first degrees from 180 downwards correspond to the cingular pole, for Primate1
    :param dimRect: the dimensions of the rectangle for the Primate1 as [Longitude,latitude]
    :param sulci: a tuple of two lists of floats, respectively, the coordinates of the longitudinal and the
     latitudinal sulci for Primate1 on the rectangle
    :return: two lists of floats, respectively, the coordinates of the longitudinal and the latitudinal
    sulci for the given Primate model, on the sphere.
    """

    # extract cardinals
    Nlon = len(sulci[0])
    Nlat = len(sulci[1])

    # extract rectangle dimensions and latitudes' boundaries
    L, l = dimRect
    ins, cin = lat[0], 180 - lat[1]

    # initialise lists
    newLon = np.copy(sulci[0])
    newLat = np.copy(sulci[1])

    for i in range(Nlon):
        newLon[i] *= 360 / L
        newLon[i] += (sulci[0][i] < 0) * 360

    for i in range(Nlat):
        if ins < sulci[1][i] < cin:
            newLat[i] *= (cin - ins) / l
            newLat[i] += ins

    return newLon, newLat


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
    :return: the two lists of coordinates on the longitudes and latitudes of Primate2 and the two lists of affine
    transformations under the form [a,b] for y = ax+b for the respective longitudinal and latitudinal intervals defined
    by the coordinates given by the first two lists
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

    return longP2, latP2, long_transform, lat_transform


################################################################################
#
# main function
#
# python ComputeAffineTransformations.py P1_model P2_model P1_P2_corr
#
# output: affine_trans_P1_to_P2_side.txt, affine_trans_P2_to_P1_side.txt
# note: P1_to_P2 = transforming coordinates of P2 to P1!!
################################################################################

def main(P1_model, P2_model, P1_P2_corr):

    outputPath = os.path.split(P1_P2_corr)[0]
    side = (os.path.split(P1_model)[1].split("_")[2])
    Primate1 = (os.path.split(P1_model)[1].split("_")[1])
    Primate2 = (os.path.split(P2_model)[1].split("_")[1])

    print('reading models\' informations')

    modelP1F = P1_model
    modelP2F = P2_model
    modelP1 = read_model(modelP1F)
    modelP2 = read_model(modelP2F)
    dimRect_P1, poles_lat_P1, longID_P1, latID_P1, sulci_lon_P1, sulci_lat_P1, lon_coor_P1, lat_coor_P1 = modelP1
    dimRect_P2, poles_lat_P2, longID_P2, latID_P2, sulci_lon_P2, sulci_lat_P2, lon_coor_P2, lat_coor_P2 = modelP2

    print('reading correspondences\' table')

    corrTable = read_corr(P1_P2_corr)

    print('rescaling square coordinates to sphere coordinates')

    sulciP1 = SquareToSphere(dimRect_P1, [lon_coor_P1, lat_coor_P1], poles_lat_P1)
    sulciP2 = SquareToSphere(dimRect_P2, [lon_coor_P2, lat_coor_P2], poles_lat_P2)

    print('extracting correspondences')

    assert (len(corrTable['lon_' + Primate1]) == len(corrTable['lon_' + Primate2]) and
            len(corrTable['lat_' + Primate1]) == len(corrTable['lat_' + Primate2])), \
        "Number of corresponding sulci do not match in corr text file."

    Ncorr_lon = len(corrTable['lon_' + Primate1])
    Ncorr_lat = len(corrTable['lat_' + Primate1])
    long_corr = np.zeros((Ncorr_lon, 2)).astype('int')
    lat_corr = np.zeros((Ncorr_lat, 2)).astype('int')
    for i in range(Ncorr_lon):
        long_corr[i][0] = longID_P1[sulci_lon_P1[corrTable['lon_' + Primate1][i]]]
        long_corr[i][1] = longID_P2[sulci_lon_P2[corrTable['lon_' + Primate2][i]]]
    for i in range(Ncorr_lat):
        lat_corr[i][0] = latID_P1[sulci_lat_P1[corrTable['lat_' + Primate1][i]]]
        lat_corr[i][1] = latID_P2[sulci_lat_P2[corrTable['lat_' + Primate2][i]]]

    print('computing affine transformations from ' + Primate1 + ' to ' + Primate2)
    P1toP2 = Affine_Transform(sulciP1, sulciP2, long_corr, lat_corr, poles_lat_P1, poles_lat_P2)

    print('writing it down')

    if not os.path.exists(outputPath + Primate1 + '_to_' + Primate2):
        os.mkdir(outputPath + Primate1 + '_to_' + Primate2)

    f1 = open(os.path.join(outputPath + Primate1 + '_to_' + Primate2,
                           'affine_trans_' + Primate1 + '_to_' + Primate2 + '_' + side), 'w+')
    f1.write(Primate1 + ' to ' + Primate2 + '\n')

    f1.write('int_lon_' + Primate2 + ':')
    for i in range(len(P1toP2[0]) - 1):
        f1.write(str(P1toP2[0][i]) + ',')
    f1.write(str(P1toP2[0][-1]) + '\n')

    f1.write('int_lat_' + Primate2 + ':')
    for i in range(len(P1toP2[1]) - 1):
        f1.write(str(P1toP2[1][i]) + ',')
    f1.write(str(P1toP2[1][-1]) + '\n')

    f1.write('long_transform:')
    for i in range(len(P1toP2[2]) - 1):
        f1.write(str(P1toP2[2][i][0]) + ' ' + str(P1toP2[2][i][1]) + ',')
    f1.write(str(P1toP2[2][-1][0]) + ' ' + str(P1toP2[2][-1][1]) + '\n')

    f1.write('lat_transform:')
    for i in range(len(P1toP2[3]) - 1):
        f1.write(str(P1toP2[3][i][0]) + ' ' + str(P1toP2[3][i][1]) + ',')
    f1.write(str(P1toP2[3][-1][0]) + ' ' + str(P1toP2[3][-1][1]) + '\n')

    f1.close()

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

    if not os.path.exists(outputPath + Primate2 + '_to_' + Primate1):
        os.mkdir(outputPath + Primate2 + '_to_' + Primate1)

    f2 = open(os.path.join(outputPath + Primate2 + '_to_' + Primate1,
                           'affine_trans_' + Primate2 + '_to_' + Primate1 + '_' + side), 'w+')
    f2.write(Primate2 + ' to ' + Primate1 + '\n')

    f2.write('int_lon_' + Primate1 + ':')
    for i in range(len(P2toP1[0]) - 1):
        f2.write(str(P2toP1[0][i]) + ',')
    f2.write(str(P2toP1[0][-1]) + '\n')

    f2.write('int_lat_' + Primate1 + ':')
    for i in range(len(P2toP1[1]) - 1):
        f2.write(str(P2toP1[1][i]) + ',')
    f2.write(str(P2toP1[1][-1]) + '\n')

    f2.write('long_transform:')
    for i in range(len(P2toP1[2]) - 1):
        f2.write(str(P2toP1[2][i][0]) + ' ' + str(P2toP1[2][i][1]) + ',')
    f2.write(str(P2toP1[2][-1][0]) + ' ' + str(P2toP1[2][-1][1]) + '\n')

    f2.write('lat_transform:')
    for i in range(len(P2toP1[3]) - 1):
        f2.write(str(P2toP1[3][i][0]) + ' ' + str(P2toP1[3][i][1]) + ',')
    f2.write(str(P2toP1[3][-1][0]) + ' ' + str(P2toP1[3][-1][1]) + '\n')

    f2.close()

    print('done')

    return None


if __name__ == '__main__':
    P1_model, P2_model, P1_P2_corr = sys.argv[1:]
    main(P1_model, P2_model, P1_P2_corr)

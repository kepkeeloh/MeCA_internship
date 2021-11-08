import numpy as np
import sys
import os
from ComputeAffineTransformations import SquareToSphere
from read_file import *


def inv_affine(param, y):
    """
    Computes the inverse image of a point through a given affine transformation.
    :param param: parameters of the affine transformation [a,b] such that y = ax+b
    :return: x = (y-b)/a
    """
    a = param[0]
    b = param[1]
    return (y - b) / a


def Affine_transform(corraxesP1, corraxesP2):
    """
    Given the coordinates of the sulcal lines on the sphere (i.e. cortical surface) for each species (Primate1 and
     Primate2), computes and returns the affine transformations on each interval between the corresponding axes
     for longitudinal and latitudinal sulci from Primate1 to Primate2 (i.e. meant to modify the Primate2's coordinate's
     system)
    :param corraxesP1: a list of two lists of floats, respectively, the coordinates of the longitudinal and the
     latitudinal sulci for Primate1 on the sphere which define the intervals (i.e. have a correspondence)
    :param corraxesP2: same thing as corraxesP1 for Primate2
    (corraxesP1[i] and corraxesP2[i] must have same length for i=0,1)
    :return: the two lists of affine transformations under the form [a,b] for y = ax+b for the respective longitudinal
    and latitudinal intervals
    """

    # extract cardinals
    Nlong = len(corraxesP1[0])
    Nlat = len(corraxesP1[1])

    # initialization of the lists
    long_transform = np.zeros((Nlong - 1, 2))
    lat_transform = np.zeros((Nlat - 1, 2))

    longP1, latP1 = corraxesP1
    longP2, latP2 = corraxesP2

    for i in range(Nlong - 1):
        long_transform[i][0] = (longP1[i + 1] - longP1[i]) / (longP2[i + 1] - longP2[i])
        long_transform[i][1] = longP1[i] - longP2[i] * long_transform[i][0]

    for i in range(Nlat - 1):
        lat_transform[i][0] = (latP1[i + 1] - latP1[i]) / (latP2[i + 1] - latP2[i])
        lat_transform[i][1] = latP1[i] - latP2[i] * lat_transform[i][0]

    return long_transform, lat_transform


def Affine_composition(modelP1, modelP2, modelP3, corrTableP12, corrTableP23, corrTableP13):
    """
    Enables to reparametrize the Primate3 coordinate system to the Primate1's through a composition with a third model
    Primate2 in order to refine the reparametrization with intermediary correspondences that might exist between
    Primate1 and Primate2 or Primate2 and Primate3 but not between Primate1 and Primate3
    @:param modelP1: the hip-hop sulci model for Primate1  given by the function read_model in read_file.py
    @:param modelP2: same thing as modelP1 for Primate2
    @:param modelP3: same thing as modelP1 for Primate3
    @:param corrTableP12: dictionary of the corresponding longitudinal and latitudinal sulci between Primate1 and
    Primate2 given by the function read_corr from read_file.py
    @:param corrTableP23: same thing as corrTableP12 between Primate2 and Primate3
    @:param corrTableP13: same thing as corrTableP13 between Primate1 and Primate3
    :return: the two lists of coordinates that define the intervals on the longitudes and latitudes of Primate3 and the
    two lists of affine transformations under the form [a,b] for y = ax+b for the respective longitudinal and
    latitudinal intervals
    """

    Primate1 = (os.path.split(P1_model)[1].split("_")[1])
    Primate2 = (os.path.split(P2_model)[1].split("_")[1])
    Primate3 = (os.path.split(P3_model)[1].split("_")[1])

    dimRect_P1, poles_lat_P1, longID_P1, latID_P1, sulci_lon_P1, sulci_lat_P1, lon_coor_P1, lat_coor_P1 = modelP1
    dimRect_P2, poles_lat_P2, longID_P2, latID_P2, sulci_lon_P2, sulci_lat_P2, lon_coor_P2, lat_coor_P2 = modelP2
    dimRect_P3, poles_lat_P3, longID_P3, latID_P3, sulci_lon_P3, sulci_lat_P3, lon_coor_P3, lat_coor_P3 = modelP3
    # extract latitudes' boundaries
    insP1, cinP1 = poles_lat_P1[0], 180 - poles_lat_P1[1]
    insP2, cinP2 = poles_lat_P2[0], 180 - poles_lat_P2[1]
    insP3, cinP3 = poles_lat_P3[0], 180 - poles_lat_P3[1]

    print('rescaling square coordinates to sphere coordinates')

    sulciP1 = np.array(SquareToSphere(dimRect_P1, [lon_coor_P1, lat_coor_P1], poles_lat_P1))
    sulciP2 = np.array(SquareToSphere(dimRect_P2, [lon_coor_P2, lat_coor_P2], poles_lat_P2))
    sulciP3 = np.array(SquareToSphere(dimRect_P3, [lon_coor_P3, lat_coor_P3], poles_lat_P3))

    print('computing affine transformations from ' + Primate1 + ' to ' + Primate2)

    Ncorr_lonP12 = len(corrTableP12['lon_' + Primate1])
    Ncorr_latP12 = len(corrTableP12['lat_' + Primate1])
    long_corrP12 = np.zeros((Ncorr_lonP12, 2)).astype('int')
    lat_corrP12 = np.zeros((Ncorr_latP12, 2)).astype('int')
    for i in range(Ncorr_lonP12):
        long_corrP12[i][0] = longID_P1[sulci_lon_P1[corrTableP12['lon_' + Primate1][i]]]
        long_corrP12[i][1] = longID_P2[sulci_lon_P2[corrTableP12['lon_' + Primate2][i]]]
    for i in range(Ncorr_latP12):
        lat_corrP12[i][0] = latID_P1[sulci_lat_P1[corrTableP12['lat_' + Primate1][i]]]
        lat_corrP12[i][1] = latID_P2[sulci_lat_P2[corrTableP12['lat_' + Primate2][i]]]

    # make the lists of the axes that have a correspondence (and therefore define the intervals)
    longP12 = np.sort(np.concatenate(([0], sulciP1[0][long_corrP12[:, 0]], [360])))
    latP12 = np.sort(np.concatenate(([insP1], sulciP1[1][lat_corrP12[:, 0]], [cinP1])))
    longP21 = np.sort(np.concatenate(([0], sulciP2[0][long_corrP12[:, 1]], [360])))
    latP21 = np.sort(np.concatenate(([insP2], sulciP2[1][lat_corrP12[:, 1]], [cinP2])))

    P1toP2 = Affine_transform([longP12, latP12], [longP21, latP21])

    print('computing affine transformations from ' + Primate2 + ' to ' + Primate3)

    Ncorr_lonP23 = len(corrTableP23['lon_' + Primate2])
    Ncorr_latP23 = len(corrTableP23['lat_' + Primate2])
    long_corrP23 = np.zeros((Ncorr_lonP23, 2)).astype('int')
    lat_corrP23 = np.zeros((Ncorr_latP23, 2)).astype('int')
    for i in range(Ncorr_lonP23):
        long_corrP23[i][0] = longID_P2[sulci_lon_P2[corrTableP23['lon_' + Primate2][i]]]
        long_corrP23[i][1] = longID_P3[sulci_lon_P3[corrTableP23['lon_' + Primate3][i]]]
    for i in range(Ncorr_latP23):
        lat_corrP23[i][0] = latID_P2[sulci_lat_P2[corrTableP23['lat_' + Primate2][i]]]
        lat_corrP23[i][1] = latID_P3[sulci_lat_P3[corrTableP23['lat_' + Primate3][i]]]

    # make the lists of the axes that have a correspondence (and therefore define the intervals)
    longP23 = np.sort(np.concatenate(([0], sulciP2[0][long_corrP23[:, 0]], [360])))
    latP23 = np.sort(np.concatenate(([insP2], sulciP2[1][lat_corrP23[:, 0]], [cinP2])))
    longP32 = np.sort(np.concatenate(([0], sulciP3[0][long_corrP23[:, 1]], [360])))
    latP32 = np.sort(np.concatenate(([insP3], sulciP3[1][lat_corrP23[:, 1]], [cinP3])))

    P2toP3 = Affine_transform([longP23, latP23], [longP32, latP32])

    print('computing the composition from ' + Primate1 + ' to ' + Primate3)

    Ncorr_lonP13 = len(corrTableP13['lon_' + Primate1])
    Ncorr_latP13 = len(corrTableP13['lat_' + Primate1])
    long_corrP13 = np.zeros((Ncorr_lonP13, 2)).astype('int')
    lat_corrP13 = np.zeros((Ncorr_latP13, 2)).astype('int')
    for i in range(Ncorr_lonP13):
        long_corrP13[i][0] = longID_P1[sulci_lon_P1[corrTableP13['lon_' + Primate1][i]]]
        long_corrP13[i][1] = longID_P3[sulci_lon_P3[corrTableP13['lon_' + Primate3][i]]]
    for i in range(Ncorr_latP13):
        lat_corrP13[i][0] = latID_P1[sulci_lat_P1[corrTableP13['lat_' + Primate1][i]]]
        lat_corrP13[i][1] = latID_P3[sulci_lat_P3[corrTableP13['lat_' + Primate3][i]]]

    # make the lists of the axes that have a correspondence (and therefore define the intervals)
    longP13 = np.sort(np.concatenate(([0], sulciP1[0][long_corrP13[:, 0]], [360])))
    latP13 = np.sort(np.concatenate(([insP1], sulciP1[1][lat_corrP13[:, 0]], [cinP1])))
    longP31 = np.sort(np.concatenate(([0], sulciP3[0][long_corrP13[:, 1]], [360])))
    latP31 = np.sort(np.concatenate(([insP3], sulciP3[1][lat_corrP13[:, 1]], [cinP3])))

    for i in range(Ncorr_lonP12):
        if corrTableP12['lon_' + Primate1][i] not in corrTableP13['lon_' + Primate1]:
            j = 0
            while sulciP2[0][long_corrP12[i][1]] >= longP23[j + 1]:
                j += 1
            longP13 = np.append(longP13, sulciP1[0][long_corrP12[i][1]])
            longP31 = np.append(longP31, inv_affine(P2toP3[0][j], sulciP2[0][long_corrP12[i][1]]))
    for i in range(Ncorr_latP12):
        if corrTableP12['lat_' + Primate1][i] not in corrTableP13['lat_' + Primate1]:
            j = 0
            while sulciP2[1][lat_corrP12[i][1]] >= latP23[j + 1]:
                j += 1
            latP13 = np.append(latP13, sulciP1[1][lat_corrP12[i][1]])
            latP31 = np.append(latP31, inv_affine(P2toP3[1][j], sulciP2[1][lat_corrP12[i][1]]))

    for i in range(Ncorr_lonP23):
        if corrTableP23['lon_' + Primate3][i] not in corrTableP13['lon_' + Primate3]:
            j = 0
            while sulciP2[0][long_corrP23[i][0]] >= longP12[j + 1]:
                j += 1
            longP13 = np.append(longP13, P1toP2[0][j][0] * sulciP2[0][long_corrP23[i][0]] + P1toP2[0][j][1])
            longP31 = np.append(longP31, sulciP3[0][long_corrP23[i][1]])
    for i in range(Ncorr_latP23):
        if corrTableP23['lat_' + Primate3][i] not in corrTableP13['lat_' + Primate3]:
            j = 0
            while sulciP2[1][lat_corrP23[i][0]] >= latP12[j + 1]:
                j += 1
            latP13 = np.append(latP13, P1toP2[1][j][0] * sulciP2[1][lat_corrP23[i][0]] + P1toP2[1][j][1])
            latP31 = np.append(latP31, sulciP3[1][lat_corrP23[i][1]])

    longP13, latP13, longP31, latP31 = np.sort(longP13), np.sort(latP13), np.sort(longP31), np.sort(latP31)

    long_transform, lat_transform = Affine_transform([longP13, latP13], [longP31, latP31])

    return longP31, latP31, long_transform, lat_transform

################################################################################
#
# main function
#
# python ComputeAffineTransformations_3species.py P1_model P2_model P3_model P1_P2_corr, P2_P3_corr, P1_P3_corr
#
# output: affine_trans_P1_to_P2_side.txt, affine_trans_P2_to_P1_side.txt
# note: P1_to_P2 = transforming coordinates of P2 to P1!!
################################################################################


def main(P1_model, P2_model, P3_model, P1_P2_corr, P2_P3_corr, P1_P3_corr):
    
    outputPath = os.path.split(P1_P2_corr)[0]
    side = (os.path.split(P1_model)[1].split("_")[2])
    Primate1 = (os.path.split(P1_model)[1].split("_")[1])
    Primate2 = (os.path.split(P2_model)[1].split("_")[1])
    Primate3 = (os.path.split(P3_model)[1].split("_")[1])
    
    print('reading models\' informations')

    modelP1F = P1_model
    modelP2F = P2_model
    modelP3F = P3_model
    modelP1 = read_model(modelP1F)
    modelP2 = read_model(modelP2F)
    modelP3 = read_model(modelP3F)

    print('reading correspondences\' table')

    corrTableP12 = read_corr(P1_P2_corr)
    corrTableP23 = read_corr(P2_P3_corr)
    corrTableP13 = read_corr(P1_P3_corr)

    P1toP3 = Affine_composition(modelP1, modelP2, modelP3, corrTableP12, corrTableP23, corrTableP13)

    if not os.path.exists(outputPath + Primate1 + '_to_' + Primate3 + '_via' + Primate2):
        os.mkdir(outputPath + Primate1 + '_to_' + Primate3 + '_via' + Primate2)

    dir = outputPath + Primate1 + '_to_' + Primate3 + '_via' + Primate2

    file = open(os.path.join(dir, 'affine_trans_' + Primate1 + '_to_' + Primate3 + '_via' + Primate2 + '_' + side + '.txt'), 'w+')
    file.write(Primate1 + ' to ' + Primate3 + ' via ' + Primate2 + '\n')

    file.write('int_lon_' + Primate3 + ':')
    for i in range(len(P1toP3[0]) - 1):
        file.write(str(P1toP3[0][i]) + ',')
    file.write(str(P1toP3[0][-1]) + '\n')

    file.write('int_lat_' + Primate3 + ':')
    for i in range(len(P1toP3[1]) - 1):
        file.write(str(P1toP3[1][i]) + ',')
    file.write(str(P1toP3[1][-1]) + '\n')

    file.write('long_transform:')
    for i in range(len(P1toP3[2]) - 1):
        file.write(str(P1toP3[2][i][0]) + ' ' + str(P1toP3[2][i][1]) + ',')
    file.write(str(P1toP3[2][-1][0]) + ' ' + str(P1toP3[2][-1][1]) + '\n')

    file.write('lat_transform:')
    for i in range(len(P1toP3[3]) - 1):
        file.write(str(P1toP3[3][i][0]) + ' ' + str(P1toP3[3][i][1]) + ',')
    file.write(str(P1toP3[3][-1][0]) + ' ' + str(P1toP3[3][-1][1]) + '\n')

    file.close()

    print('done')

    return None


if __name__ == '__main__':
    P1_model, P2_model, P3_model, P1_P2_corr, P2_P3_corr, P1_P3_corr = sys.argv[1:]
    main(P1_model, P2_model, P3_model, P1_P2_corr, P2_P3_corr, P1_P3_corr)

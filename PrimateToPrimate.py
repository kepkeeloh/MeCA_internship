import numpy as np
from soma import aims
import sys
from read_file import *

"""
This code aims to rescale the axes' coordinates of Primate2 in order to be able to map the Primate1's texture onto 
the Primate2 surface.
"""


def SquareToSphere(dimRectP1, dimRectP2, sulciP1, sulciP2):
    """
    As the system of coordinates is not the same in both the rectangle and the sphere for the models we are considering,
    we need to compute the position of the axes of each species (Primate1 and Primate2) on the sphere (i.e. the cortical
    surface), given the respective dimensions of the rectangles of each species and the wanted dimensions of the sphere
    (which might be different from ([0,360]*[0,180]) as we do not take into account the intervals where the poles are
    on the latitudes: with the default value, below 30 and above 150 are respectively the insular and cingular poles).
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

    # extract rectangle dimensions
    LP1, lP1 = dimRectP1
    LP2, lP2 = dimRectP2

    # initialise lists
    newLon_P1 = np.copy(sulciP1[0])
    newLat_P1 = np.copy(sulciP1[1])
    newLon_P2 = np.copy(sulciP2[0])
    newLat_P2 = np.copy(sulciP2[1])

    for i in range(Nlon_P1):
        newLon_P1[i] *= 360 / LP1
        newLon_P1[i] += (sulciP1[0][i] < 0) * 360

    for i in range(Nlat_P1):
        if 30 < sulciP1[1][i] < 150:
            newLat_P1[i] *= 120 / lP1
            newLat_P1[i] += 30

    for i in range(Nlon_P2):
        newLon_P2[i] *= 360 / LP2
        newLon_P2[i] += (sulciP2[0][i] < 0) * 360

    for i in range(Nlat_P2):
        if 30 < sulciP1[1][i] < 150:
            newLat_P2[i] *= 120 / lP2
            newLat_P2[i] += 30

    return [newLon_P1, newLat_P1], [newLon_P2, newLat_P2]


def Affine_Transform(sulciP1, sulciP2, long_corr, lat_corr):
    """
    Given the coordinates of the sulcal lines on the sphere (i.e. cortical surface) for each species (Primate1 and
     Primate2), computes and returns the affine transformations on each interval between the corresponding axes
     for longitudinal and latitudinal sulci from Primate1 to Primate2 (i.e. meant to modify the Primate2's coordinate's
     system)
    :param sulciP1: a list of two lists of floats, respectively, the coordinates of the longitudinal and the
     latitudinal sulci for Primate1 on the sphere
    :param sulciP2: same thing as sulciP1 for Primate2
    :param long_corr: a list of couples of integers giving the indices of the corresponding axes on the longitudes,
    with the first index corresponding to Primate1 and the second for Primate2
    :param lat_corr: same thing as long_corr for the latitudinal sulci
    :return: two lists of couples of floats corresponding to the parameters (a,b) of the affine transformations
    y = ax + b on each interval, respectively for the longitudinal and latitudinal sulci.
    """

    # extract cardinals
    Nlong = len(long_corr)
    Nlat = len(lat_corr)

    # initialization of the lists
    long_transform = np.zeros((Nlong + 1, 2))
    lat_transform = np.zeros((Nlat + 1, 2))

    # make the lists of the axes that have a correspondence (and therefore define the intervals)
    longP1 = np.sort(np.concatenate(([0], sulciP1[0][long_corr[:, 0]], [360])))
    latP1 = np.sort(np.concatenate(([30], sulciP1[1][lat_corr[:, 0]], [150])))
    longP2 = np.sort(np.concatenate(([0], sulciP2[0][long_corr[:, 1]], [360])))
    latP2 = np.sort(np.concatenate(([30], sulciP2[1][lat_corr[:, 1]], [150])))

    for i in range(Nlong + 1):
        long_transform[i][0] = (longP1[i + 1] - longP1[i]) / (longP2[i + 1] - longP2[i])
        long_transform[i][1] = longP1[i] - longP2[i] * long_transform[i][0]

    for i in range(Nlat + 1):
        lat_transform[i][0] = (latP1[i + 1] - latP1[i]) / (latP2[i + 1] - latP2[i])
        lat_transform[i][1] = latP1[i] - latP2[i] * lat_transform[i][0]

    return long_transform, lat_transform


def rescale(sulci, affine, intervals):
    """
    Updates the sulci coordinates on every interval between the corresponding axes thanks to the given affine
    transformations. It can be used for either longitudinal or latitudinal sulcal lines and is meant to be used on
    Primate2 for a mapping from Primate1 to Primate2.
    :param sulci: the list of coordinates of the sulcal lines
    :param affine: the list of affine transformations under the form (a,b) for y = ax + b
    :param intervals: the list of coordinates that define the intervals of transformation (i.e. the list of coordinates
    which have a correspondence with the other primate species we are comparing it to)
    :return: the updated coordinates of the sulcal lines
    """

    N = len(sulci)  # there are N sulci, N+1 intervals, hence N+1 affine transformations
    rescaled = np.copy(sulci)  # initialize the list

    for j in range (N):  # first N intervals
        for i in range (len(intervals)-1):
            if intervals[i] < sulci[j] < intervals[i+1]:
                rescaled[j] *= affine[i][0]
                rescaled[j] += affine[i][1]

    return rescaled


####################################################################
#
# main function
#
# python PrimateToPrimate.py Primate1 Primate2 side
#
####################################################################

def main(Primate1, Primate2, side):
    nameLon = Primate2 + '_' + side + 'white_lon.gii'
    nameLat = Primate2 + '_' + side + 'white_lat.gii'

    print('reading models\' informations')

    modelP1F = 'model_' + Primate1 + '_' + side + '.txt'
    modelP2F = 'model_' + Primate2 + '_' + side + '.txt'
    modelP1 = read_model(modelP1F)
    modelP2 = read_model(modelP2F)
    dimRect_P1, longID_P1, latID_P1, sulci_lon_P1, sulci_lat_P1, lon_coor_P1, lat_coor_P1 = modelP1
    dimRect_P2, longID_P2, latID_P2, sulci_lon_P2, sulci_lat_P2, lon_coor_P2, lat_coor_P2 = modelP2

    print('reading correspondences\' table')

    name_corr = Primate1 + '_' + Primate2 + '_' + 'corr.txt'
    corrTable = read_corr(name_corr)

    print('reading coordinates')

    r = aims.Reader()
    texLatF = r.read(nameLat)
    texLonF = r.read(nameLon)
    texLat = np.array(texLatF[0])
    texLon = np.array(texLonF[0])

    print('rescaling square coordinates to sphere coordinates')

    sulciP1, sulciP2 = SquareToSphere(dimRect_P1, dimRect_P2, [lon_coor_P1, lat_coor_P1],
                                      [lon_coor_P2, lat_coor_P2])

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

    print('computing affine transformations')

    long_transform, lat_transform = Affine_Transform(sulciP1, sulciP2, long_corr, lat_corr)

    print('processing longitude')

    intervals_lon = np.concatenate(([0], sulciP2[0][long_corr[:, 1]], [360]))
    intervals_lon = np.sort(intervals_lon)
    newLon = rescale(texLon, long_transform, intervals_lon)

    print('processing latitude')

    intervals_lat = np.concatenate(([30], sulciP2[1][lat_corr[:, 1]], [150]))
    intervals_lat = np.sort(intervals_lat)
    newLat = rescale(texLat, lat_transform, intervals_lat)

    print('writing textures')

    nv = texLat.size
    newLatT = aims.TimeTexture_FLOAT(1, nv)
    newLonT = aims.TimeTexture_FLOAT(1, nv)

    for i in range(nv):
        newLatT[0][i] = newLat[i]
        newLonT[0][i] = newLon[i]

    outLat = Primate1 + '_' + side + 'white_lat_to' + Primate2 + '.gii'
    outLon = Primate1 + '_' + side + 'white_lon_to' + Primate2 + '.gii'

    r = aims.Writer()
    r.write(newLatT, outLat)
    r.write(newLonT, outLon)
    print('done')

    return None


if __name__ == '__main__':
    Primate1, Primate2, side = sys.argv[1:]
    main(Primate1, Primate2, side)

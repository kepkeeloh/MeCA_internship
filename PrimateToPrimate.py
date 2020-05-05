import numpy as np
from soma import aims
import sys


def SquareToSphere(dimP1, dimP2, sulciP1, sulciP2, lat_sphere=[30., 150.]):
    """
    As the system of coordinates is not the same in both the rectangle and the sphere for the models we are considering,
    we need to compute the position of the axes of each species (Primate1 and Primate2) on the sphere (i.e. the cortical
    surface), given the respective dimensions of the rectangles of each species and the wanted dimensions of the sphere
    (which might be different from ([0,360]*[0,180]) as we do not take into account the intervalles where the poles are
    on the latitudes: with the default value, below 30° and above 150° are respectively the insular and cingular poles).
    :param dimP1: the dimensions of the rectangle for the Primate1 ('[Longitude,latitude]')
    :param dimP2: same thing as dimP1 for Primate2
    :param sulciP1: a tuple of two lists of floats, respectively, the coordinates of the longitudinal and the
     latitudinal sulci for Primate1 on the rectangle
    :param sulciP2: same thing as sulciP2 for Primate2
    :param dim_sphere: the dimensions of the sphere without the poles' intervalles.
    :return: a tuple of four lists of floats, respectively, the coordinates of the longitudinal and the latitudinal
    sulci for Primate1 and the coordinates of the longitudinal and the latitudinal sulci for Primate2, on the sphere.
    """

    # extract cardinals
    Nlong_P1 = len(sulciP1[0])
    Nlat_P1 = len(sulciP1[1])
    Nlong_P2 = len(sulciP2[0])
    Nlat_P2 = len(sulciP2[0])

    # extract dimensions
    LP1, lP1 = dimP1
    LP2, lP2 = dimP2
    lat_min, lat_max = lat_sphere
    range_lat = lat_max - lat_min

    # initialization of the lists
    longS_P1 = np.copy(sulciP1[0])
    latS_P1 = np.copy(sulciP1[1])
    longS_P2 = np.copy(sulciP2[0])
    latS_P2 = np.copy(sulciP2[1])

    for i in range(Nlong_P1):
        longS_P1[i] *= 360 / LP1
        longS_P1[i] += (sulciP1[0][i] < 0) * 360

    for i in range(Nlat_P1):
        if (sulciP1[1][i] > lat_min and sulciP1[1][i] < lat_max):
            latS_P1[i] *= range_lat / lP1
            latS_P1[i] += 30

    for i in range(Nlong_P2):
        longS_P2[i] *= 360 / LP2
        longS_P2[i] += (sulciP2[0][i] < 0) * 360

    for i in range(Nlat_P2):
        if (sulciP2[1][i] > lat_min and sulciP2[1][i] < lat_max):
            latS_P2[i] *= range_lat / lP2
            latS_P2[i] += 30

    return (longS_P1, latS_P1, longS_P2, latS_P2)


def Affine_Transform(sulciP1, sulciP2, long_corr, lat_corr):
    """
    Given the coordinates of the sulcal lines on the sphere (i.e. cortical surface) for each species (Primate1 and
     Primate2), computes and returns the affine transformations on each interval between the corresponding axes
     for longitudinal and latitudinal sulci from Primate1 to Primate2.
    :param sulciP1: a tuple of two lists of floats, respectively, the coordinates of the longitudinal and the
     latitudinal sulci for Primate1 on the sphere
    :param sulciP2: same thing as sulciP1 for Primate2
    :param long_corr: list of couples of integers giving the indices of the corresponding axes on the longitudes,
    with the first index corresponding to Primate1 and the second for Primate2
    :param lat_corr: same thing as long_corr for the latitudinal sulci
    :return: two lists of couples of floats corresponding to the parameters (a,b) of the affine transformations
    y = ax + b on each intervalle, respectively for the longitudinal and latitudinal sulci.
    """

    # extract cardinals
    Nlong = len(long_corr)
    Nlat = len(lat_corr)

    # initialization of the lists
    long_transform = np.zeros((Nlong, 2))
    lat_transform = np.zeros((Nlat, 2))

    # make the lists of the axes that have a correspondence (and therefore define the intervals)
    longP1 = sulciP1[0][long_corr[:, 0]]
    latP1 = sulciP1[1][lat_corr[:, 0]]
    longP2 = sulciP1[0][long_corr[:, 1]]
    latP2 = sulciP2[1][lat_corr[:, 1]]

    # the transformation on the first interval is linear
    long_transform[0][0] = sulciP2[0][0] / sulciP1[0][0]
    lat_transform[0][0] = sulciP2[1][0] / sulciP1[1][0]

    for i in range(1, Nlong):
        long_transform[i][0] = (longP2[i + 1] - longP2[i]) / (longP1[i + 1] - longP1[i])
        long_transform[i][1] = longP2[i] - longP1[i] * long_transform[i][0]

    for i in range(1, Nlat):
        lat_transform[i][0] = (latP2[i + 1] - latP2[i]) / (latP1[i + 1] - latP1[1][i])
        lat_transform[i][1] = latP2[i] - latP1[i] * lat_transform[i][0]

    return long_transform, lat_transform


def rescale(sulci, affine, intervals):
    """
    Updates the sulci coordinates on every interval between the corresponding axes
    thanks to the given affine transformations.
    Remark: can be used for either longitudinal or latitudinal sulcal lines.
    :param lat: Boolean indicating whether we are looking at the latitudes (TRUE) or the longitudes (FALSE)
    :param sulci: the list of coordinates of the sulcal lines
    :param affine: the list of affine transformations under the form (a,b) for y = ax + b
    :param intervals: the list of coordinates that define the intervals of transformation (i.e. the list of coordinates
    which have a correspondence with the other primate species we are comparing it to)
    :return: the updated coordinates of the sulci
    """

    assert (len(affine) == len(intervals) + 1), "The lengths of the affine transformations and the intervals " \
                                                "do not match."

    new_intervals = np.concatenate([0], intervals, [360]) # we add the boundaries to the list of coordinates
    N = len(sulci)
    rescaled = np.zeros(N)
    i = 0
    j = 0
    while i < N:
        while (sulci[j] >= new_intervals[i]) and (sulci[j] < new_intervals[i+1]):
            rescaled[j] *= affine[i][0]
            rescaled[j] += affine[i][1]
            j += 1
        i += 1
    return rescaled


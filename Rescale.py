import numpy as np
import sys
from soma import aims
from read_file import *


def rescale(texture, affine, intervals):
    """
    Updates the axes' coordinates on every interval between the corresponding axes thanks to the given affine
    transformations. It can be used for either longitudinal or latitudinal sulcal lines and is meant to be used on
    Primate2 for a mapping from Primate1 to Primate2.
    :param texture: numpy array of coordinates of the (longitudinal or latitudinal) axes
    :param affine: the list of affine transformations under the form (a,b) for y = ax + b
    :param intervals: the list of coordinates that define the intervals of transformation (i.e. the list of coordinates
    which have a correspondence with the other primate species we are comparing it to)
    :return: the updated coordinates of the axes
    """

    N = len(texture)
    M = len(intervals) - 1
    rescaled = np.copy(texture)  # initialize the list

    for j in range(N):
        for i in range(M):
            if intervals[i] < texture[j] < intervals[i + 1]:
                rescaled[j] *= affine[i][0]
                rescaled[j] += affine[i][1]

    return rescaled

####################################################################
#
# main function
#
# python Rescale.py Primate1 Primate2 side
#
####################################################################

def main(Primate1, Primate2, side):
    """
    Main function that uses the text file returned by the Affine_transformations.py code in order to rescale the
    Primate2's texture so we can map Primate1's textures onto Primate1's brain surface.
    :param Primate1: string - name of Primate1
    :param Primate2: string - name of Primate2
    :param side: string - hemisphere 'L' or 'R'
    :return: the rescaled Primate2's longitudinal and latitudinal textures
    """

    nameLon = Primate2 + '_' + side + 'white_lon.gii'
    nameLat = Primate2 + '_' + side + 'white_lat.gii'

    print('reading coordinates')

    r = aims.Reader()
    texLatF = r.read(nameLat)
    texLonF = r.read(nameLon)
    texLat = np.array(texLatF[0])
    texLon = np.array(texLonF[0])

    print('reading affine transformations')

    affine_model = 'affine_trans_' + Primate1 + '_to_' + Primate2 + '_' + side + '.txt'
    int_lon, int_lat, lon_transform, lat_transform = read_affine(affine_model)

    print('processing longitude')

    newLon = rescale(texLon, lon_transform, int_lon)

    print('processing latitude')

    newLat = rescale(texLat, lat_transform, int_lat)

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

if __name__ == '__main__':
    Primate1, Primate2, side = sys.argv[1:]
    main(Primate1, Primate2, side)
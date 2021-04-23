#Edited from jeanne's rescale.py

import os
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
# python RescaleLongLatTextures.py affine_transformation_P1_to_P2 LonTex_P2 LatTex_P2
#
####################################################################

def main(affine_transformation, LonTex, LatTex):
    """
    Main function that uses the affine_transformation text file (P1_to_P2) to rescale the
    P2's texture so we can map P1's textures onto P2's brain surface.
    :param affine_transformation: transformation from P1 to P2
    :param longTex: P2 long texture that you want to rescale to P1
    :param latTex: P2 lat texture that you want to rescale to P1
    :return: the rescaled P2's longitudinal and latitudinal textures
    """
    
    outputPath = os.path.split(os.path.abspath(affine_transformation))[0]
    side = os.path.split(os.path.abspath(affine_transformation))[1].split("_")[5].split(".")[0]
    P1 = (os.path.split(affine_transformation)[1].split("_")[2])
    P2 = (os.path.split(affine_transformation)[1].split("_")[4])

    print(outputPath)
    print(side)
    print(P2)

    print('reading coordinates')

    r = aims.Reader()
    texLatF = r.read(LatTex)
    texLonF = r.read(LonTex)
    texLat = np.array(texLatF[0])
    texLon = np.array(texLonF[0])

    print('reading affine transformations')

    affine_model = affine_transformation
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

    outLat = os.path.join(outputPath, P1 + '_to_' + P2 + '_lat_' + side + '.gii')
    outLon = os.path.join(outputPath, P1 + '_to_' + P2 + '_lon_' + side + '.gii')

    r = aims.Writer()
    r.write(newLatT, outLat)
    r.write(newLonT, outLon)
    print('done')

if __name__ == '__main__':
    affine_transformation, LonTex, LatTex = sys.argv[1:]
    main(affine_transformation, LonTex, LatTex)

import numpy as np
import os
import sys
from read_file import *
from soma import aims
from Rescale import rescale

def iterate (Primate1, Primate2, side):
    """
    Meant to be used on a database to iterate the reparametrization on a list of individuals of one species to the
    coordinate system of Primate1.
    :param Primate1: string - name of Primate1
    :param Primate2: string - name of Primate2
    :param side: string - hemisphere 'L' or 'R'
    :return: the rescaled longitudinal and latitudinal textures for each individual of Primate2
    """

    individuals = [i for i in os.listdir(Primate2) if '.' not in i]

    print('reading affine transformations')

    affine_model = os.path.join(Primate1 + '_to_' + Primate2, 'affine_trans_' + Primate1 + '_to_' + Primate2 + '_' + side + '.txt')
    int_lon, int_lat, lon_transform, lat_transform = read_affine(affine_model)

    for ind in individuals:

        nameLon = os.path.join(Primate2, ind, ind + '_' + side + 'white_lon.gii')
        nameLat = os.path.join(Primate2, ind, ind + '_' + side + 'white_lat.gii')

        # print('reading coordinates')

        r = aims.Reader()
        texLatF = r.read(nameLat)
        texLonF = r.read(nameLon)
        texLat = np.array(texLatF[0])
        texLon = np.array(texLonF[0])

        # print('processing longitude')

        newLon = rescale(texLon, lon_transform, int_lon)

        # print('processing latitude')

        newLat = rescale(texLat, lat_transform, int_lat)

        # print('writing textures')

        nv = texLat.size
        newLatT = aims.TimeTexture_FLOAT(1, nv)
        newLonT = aims.TimeTexture_FLOAT(1, nv)

        for i in range(nv):
            newLatT[0][i] = newLat[i]
            newLonT[0][i] = newLon[i]

        outLat = os.path.join(Primate1 + '_to_' + Primate2, Primate1 + '_' + side + 'white_lat_to' + ind + '.gii')
        outLon = os.path.join(Primate1 + '_to_' + Primate2, Primate1 + '_' + side + 'white_lon_to' + ind + '.gii')

        r = aims.Writer()
        r.write(newLatT, outLat)
        r.write(newLonT, outLon)

    print('done')

if __name__ == '__main__':
    Primate1, Primate2, side = sys.argv[1:]
    iterate(Primate1, Primate2, side)

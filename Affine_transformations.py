import numpy as np
from soma import aims
import sys
from read_file import *
from PrimateToPrimate import *

####################################################################
#
# main function
#
# python PrimateToPrimate.py Primate1 Primate2 side
#
####################################################################

def main(Primate1, Primate2, side):

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

    print('writing affine transformations in text file')

    nameP1toP2 = Primate1 + 'To' + Primate2 + '_AffineTrans.txt'
    nameP2toP1 = Primate2 + 'To' + Primate1 + '_AffineTrans.txt'
    file = open(nameP1toP2)
    

    print('done')


    return None

if __name__ == '__main__':
    Primate1, Primate2, side = sys.argv[1:]
    main(Primate1, Primate2, side)
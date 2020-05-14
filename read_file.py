import numpy as np


def read_model(file):
    data = open(file, 'r').read()
    data = data.split('\n')
    for i in range(len(data)):
        data[i] = data[i].split(' ')

    dimRect = [round(float(data[2][1]) - float(data[1][1]), ndigits=1),
               round(float(data[3][1]) - float(data[4][1]), ndigits=1)]

    poles_lat = [float(data[5][1]), float(data[6][1])]

    for i in range(7, 13):
        data[i][1] = data[i][1].split(',')

    # Dictionaries for the axes' ID
    longID = {}
    latID = {}
    for i in range(len(data[7][1])):
        longID[data[7][1][i]] = i
    for i in range(len(data[10][1])):
        latID[data[10][1][i]] = i

    # Dictionaries for the sulci
    sulci_lon = {}
    sulci_lat = {}
    for i in range(len(data[9][1])):
        data[9][1][i] = data[9][1][i].split(':')
        data[9][1][i][1] = data[9][1][i][1][1:-1].split(';')
        sulci_lon[data[9][1][i][0]] = data[9][1][i][1][0]
    for i in range(len(data[12][1])):
        data[12][1][i] = data[12][1][i].split(':')
        data[12][1][i][1] = data[12][1][i][1][1:-1].split(';')
        sulci_lat[data[12][1][i][0]] = data[12][1][i][1][0]

    # Lists of the coordinates
    lon_coor = np.zeros(len(data[8][1]))
    lat_coor = np.zeros(len(data[11][1]))
    for i in range(len(data[8][1])):
        if data[8][1][i] != 'None':
            lon_coor[i] = float(data[8][1][i])
    for i in range(len(data[11][1])):
        if data[11][1][i] != 'None':
            lat_coor[i] = float(data[11][1][i])

    return dimRect, poles_lat, longID, latID, sulci_lon, sulci_lat, lon_coor, lat_coor


def read_corr(file):
    corrTable = open(file, 'r').read()
    corrTable = corrTable.split('\n')
    for i in range(4):
        corrTable[i] = corrTable[i].split(',')

    return corrTable


def read_affine(file):
    data = open(file, 'r').read()
    data = data.split('\n')
    for i in range(1, 5):
        data[i] = data[i].split(':')
    int_lon = np.fromstring(data[1][1], dtype='float', sep=',')
    int_lat = np.fromstring(data[2][1], dtype='float', sep=',')
    Nlon = len(int_lon) - 1
    Nlat = len(int_lat) - 1
    data[3][1] = data[3][1].split(',')
    data[4][1] = data[4][1].split(',')
    lon_transform = np.zeros((Nlon, 2))
    lat_transform = np.zeros((Nlat, 2))
    for i in range(Nlon):
        lon_transform[i] = [data[3][1][i][0], data[3][1][i][1]]
    for i in range(Nlat):
        lat_transform[i] = [data[4][1][i][0], data[4][1][i][1]]
    lon_transform = np.array(lon_transform).astype('float')
    lat_transform = np.array(lat_transform).astype('float')
    return int_lon, int_lat, lon_transform, lat_transform

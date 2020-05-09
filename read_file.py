import numpy as np

def read_model(file):

    data = open(file, 'r').read()
    data = data.split('\n')
    for i in range(len(data)):
        data[i] = data[i].split(' ')

    dimRect = [round(float(data[2][1]) - float(data[1][1]), ndigits=1),
               round(float(data[3][1]) - float(data[4][1]), ndigits=1)]

    for i in range (7,13):
        data[i][1] = data[i][1].split(',')

    # Dictionaries for the axes' ID
    longID = {}
    latID = {}
    for i in range (len(data[7][1])):
        longID[data[7][1][i]] = i
    for i in range (len(data[10][1])):
        latID[data[10][1][i]] = i

    # Dictionaries for the sulci
    sulci_lon = {}
    sulci_lat = {}
    for i in range (len(data[9][1])):
        data[9][1][i] = data[9][1][i].split(':')
        data[9][1][i][1] = data[9][1][i][1][1:-1].split(';')
        sulci_lon[data[9][1][i][0]] = data[9][1][i][1][0]
    for i in range (len(data[12][1])):
        data[12][1][i] = data[12][1][i].split(':')
        data[12][1][i][1] = data[12][1][i][1][1:-1].split(';')
        sulci_lat[data[12][1][i][0]] = data[12][1][i][1][0]

    # Lists of the coordinates
    lon_coor = np.zeros(len(data[8][1]))
    lat_coor = np.zeros(len(data[11][1]))
    for i in range (len(data[8][1])):
        lon_coor[i] = float(data[8][1][i])
    for i in range (len(data[11][1])):
        lat_coor[i] = float(data[11][1][i])

    return dimRect, longID, latID, sulci_lon, sulci_lat, lon_coor, lat_coor

def read_corr(file):

    corrTable = open(file, 'r').read()
    corrTable = corrTable.split('\n')
    for i in range (4):
        corrTable[i] = corrTable[i].split(',')

    return corrTable
def read_model(file):
    """
    Extracts information from a text file that must must have the exact same pattern (and more particularly
    for the separators and order of information) as this macaque model example:

    "modelVersion 3
    left -241.470347818
    right 228.529651456
    top 100.0
    bottom 0.0
    insularPoleBoundaryCoord 30.0
    cingularPoleBoundaryCoord 30.0
    longitudeAxisID 1,3,5,6,8,11,12,16
    longitudeAxisCoord 9.04765792564e-18,19.6207569386,-17.3644152163,-60.0518962536,-43.1257004743,-54.884437685,-51.4820049146,44.8763290846
    longitudeAxisSulci ROS.:(16;10.0),IOS.:(6;10.0),LuS.:(12;10.0),CaS.:(11;10.0),PoCS.:(5;10.0),PrOS.:(8;10.0),S.C.:(1;10.0),IPS-v.:(5;10.0),IArS.:(3;10.0),PrCS.:(3;10.0)
    latitudeAxisID 2,6,10,12,13,14
    latitudeAxisCoord 73.0011706043,54.6248021011,77.4143512621,90.0938482673,83.2011641003,67.9623418772
    latitudeAxisSulci CgS.:(12;10.0),OTS-a.:(10;10.0),PS.:(2;50.0),PmTS.:(14;10.0),PSPD.:(10;10.0),ITS.:(2;10.0),MOS-m.:(13;10.0),OTS-p.:(10;10.0),STS.:(6;10.0),IPD.:(6;10.0)
    "
    :param file: a text file as explained above
    :return: the dimensions of the rectangle, the sphere's latitude range without the poles, the dictionary of ID of the
    axes, the list of coordinates of the longitudinal and latitudinal sulci, and two dictionaries for the longitudes and
    latitudes whose keys are the sulci's names and values are lists of their ID
    and weight.
    """
    model = open(file, 'r')
    data = model.read()
    data = data.split('\n')
    if data[-1] == '':
        data = data[:-1]
    for i in range(len(data)):
        data[i] = data[i].split()

    dimRect = [round(float(data[2][1]) - float(data[1][1]), ndigits=1),
               round(float(data[3][1]) - float(data[4][1]), ndigits=1)]
    lat_sphere = [float(data[5][1]), float(data[6][1])]

    for j in range(7, len(data)):
        data[j] = [data[j][0], data[j][1].split(',')]

    longitudeAxisID = {}
    latitudeAxisID = {}
    for i in range(len(data[7][1])):
        longitudeAxisID[data[7][1][i]] = i
    for i in range(len(data[10][1])):
        latitudeAxisID[data[10][1][i]] = i
    # print(longitudeAxisID)
    # print(latitudeAxisID)

    longitudeAxisCoord = data[8][1]
    latitudeAxisCoord = data[11][1]
    for i in range(len(data[8][1])):
        if data[8][1][i] != 'None':
            longitudeAxisCoord[i] = float(data[8][1][i])
    for i in range(len(data[11][1])):
        if data[11][1][i] != 'None':
            latitudeAxisCoord[i] = float(data[11][1][i])
    # print('longitudeAxisCoord', longitudeAxisCoord)
    # print('latitudeAxisCoord', latitudeAxisCoord)

    longitudeAxisSulci = {}
    latitudeAxisSulci = {}

    for i in range(len(data[9][1])):
        data[9][1][i] = data[9][1][i].split(':')
        longitudeAxisSulci[data[9][1][i][0]] = data[9][1][i][1]
    for i in range(len(data[12][1])):
        data[12][1][i] = data[12][1][i].split(':')
        latitudeAxisSulci[data[12][1][i][0]] = data[12][1][i][1]

    for key in longitudeAxisSulci.keys():
        longitudeAxisSulci[key] = [float(i) for i in longitudeAxisSulci[key][1:-1].split(';')]
    for key in latitudeAxisSulci.keys():
        latitudeAxisSulci[key] = [float(i) for i in latitudeAxisSulci[key][1:-1].split(';')]

    # print('longitudeAxisSulci', longitudeAxisSulci)
    # print('latitudeAxisSulci', latitudeAxisSulci)

    return dimRect, lat_sphere, longitudeAxisID, latitudeAxisID, longitudeAxisCoord, latitudeAxisCoord, \
           longitudeAxisSulci, latitudeAxisSulci


def read_corr(file):
    """
    monkey long, human long, monkey lat, human lat
    :param file: the corresponding sulci between the two primate species as a .txt file that contains 4 lines for
    respectively the Primate1 longitudinal sulci, Primate2 longitudinal sulci, and the last two lines similarly for
    Primate2. (separator '\n' between the lines and ',' with no space between the names of the sulci)
    :return: two lists for the longitudinal and latitudinal correspondences, each having as many rows as correspondences
    and two columns for respectively Primate1, Primate2.
    """

    corr = open(file, 'r').read()
    corr = corr.split('\n')
    if corr[-1] == '':
        corr = corr[:-1]
    for i in range(len(corr)):
        corr[i] = corr[i].split(',')

    return corr


read_model('sulcal_coordinates_model_left.txt')

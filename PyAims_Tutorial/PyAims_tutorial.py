from __future__ import print_function
import os
from soma import aims
import numpy
import copy

src = '/Users/jeanneabitbol/Documents/MeCA_internship/demo_data/'


obj = aims.read(os.path.join(src,'data_for_anatomist/subject01/subject01.nii'))
print("obj size: ",obj.getSize())
obj2 = aims.read(os.path.join(src,'data_for_anatomist/subject01/Audio-Video_T_map.nii'))
print("obj2 size: ",obj2.getSize())
obj3 = aims.read(os.path.join(src,'data_for_anatomist/subject01/subject01_Lhemi.mesh'))
print('obj3 length(vertex(0)): ',len(obj3.vertex(0)))
assert(obj.getSize() == [256, 256, 124, 1])
assert(obj2.getSize() == [53, 63, 46, 1])
assert(obj3.size() == 1 and len(obj3.vertex(0)) == 33837)

# aims.write(obj2, 'Audio-Video_T_map.ima')
# aims.write(obj3, 'subject01_Lhemi.gii')

# create a 3D volume of signed 16-bit ints, of size 192x256x128
vol = aims.Volume(192, 256, 128, dtype='int16')
# fill it with zeros
vol.fill(0)
# set value 12 at voxel (100, 100, 60)
vol.setValue(12, 100, 100, 60)
# get value at the same position
x = vol.value(100, 100, 60)
# print(x)
assert(x == 12)

# set the voxels size
vol.header()['voxel_size'] = [0.9, 0.9, 1.2, 1.]
print('vol.header: ',vol.header())
assert(vol.header() == {'sizeX': 192, 'sizeY': 256, 'sizeZ': 128, 'sizeT': 1, 'voxel_size': [0.9, 0.9, 1.2, 1]})


# multiplication, addition etc
vol *= 2
vol2 = vol * 3 + 12
# print(vol2.value(100, 100, 60))
vol /= 2
vol3 = vol2 - vol - 12
# print(vol3.value(100, 100, 60))
vol4 = vol2 * vol / 6
# print(vol4.value(100, 100, 60))
assert(vol2.value(100, 100, 60) == 84)
assert(vol3.value(100, 100, 60) == 60)
assert(vol4.value(100, 100, 60) == 168)

# fill the volume with the distance to voxel (100, 100, 60)
vs = vol.header()['voxel_size']
pos0 = (100 * vs[0], 100 * vs[1], 60 * vs[2]) # in millimeters
for z in range(vol.getSizeZ()):
    for y in range(vol.getSizeY()):
        for x in range(vol.getSizeX()):
            # get current position in an aims.Point3df structure, in mm
            p = aims.Point3df(x * vs[0], y * vs[1], z * vs[2])
            # get relative position to pos0, in voxels
            p -= pos0
            # distance: norm of vector p
            dist = p.norm()
            # set it into the volume
            vol.setValue(dist, x, y, z)
# print(vol.value(100, 100, 60))
# save the volume
aims.write(vol, 'distance.nii')
assert(vol.value(100, 100, 60) == 0)

vol = aims.read(os.path.join(src,'data_for_anatomist/subject01/Audio-Video_T_map.nii'))
# print(vol.value(20, 20, 20) < 3. and vol.value(20, 20, 20) != 0.)
assert(vol.value(20, 20, 20) < 3. and vol.value(20, 20, 20) != 0.)
for z in range(vol.getSizeZ()):
    for y in range(vol.getSizeY()):
        for x in range(vol.getSizeX()):
            if vol.value(x, y, z) < 3.:
                vol.setValue(0, x, y, z)
# print(vol.value(20, 20, 20))
aims.write(vol, 'Audio-Video_T_thresholded.nii')
assert(vol.value(20, 20, 20) == 0.)

vol = aims.read(os.path.join(src,'data_for_anatomist/subject01/subject01.nii'))
# allocate a new volume with half dimensions
vol2 = aims.Volume(vol.getSizeX() / 2, vol.getSizeY() / 2, vol.getSizeZ() / 2, dtype='DOUBLE')
# print(vol2.getSizeX())
assert(vol2.getSizeX() == 128)
# set the voxel size to twice it was in vol
vs = vol.header()['voxel_size']
vs2 = [x * 2 for x in vs]
vol2.header()['voxel_size'] = vs2
for z in range(vol2.getSizeZ()):
    for y in range(vol2.getSizeY()):
        for x in range(vol2.getSizeX()):
            vol2.setValue(vol.value(x*2, y*2, z*2), x, y, z)
# print(vol.value(100, 100, 40))
# print(vol2.value(50, 50, 20))
aims.write(vol2, 'resampled.nii')
assert(vol.value(100, 100, 40) == 775)
assert(vol2.value(50, 50, 20) == 775.)

vol.fill(0)
arr = numpy.asarray(vol)
# set value 100 in a whole sub-volume
arr[60:120, 60:120, 40:80] = 100
# note that arr is a shared view to the volume contents,
# modifications will also affect the volume
# print(vol.value(65, 65, 42))
# print(vol.value(65, 65, 30))
aims.write(vol, "cube.nii")
assert(vol.value(65, 65, 42) == 100)
assert(vol.value(65, 65, 30) == 0)

vol.fill(0)
# vol[60:120, 60:120, 40:80] = 100 ERROR 'Volume_S16' object does not support item assignment
# print(vol.value(65, 65, 42))
# print(vol.value(65, 65, 30))
# assert(numpy.all(vol[65, 65, 42, 0] == 100)) ERROR  'Volume_S16' object is unsubscriptable
# assert(numpy.all(vol[65, 65, 30, 0] == 0))

## Let's rewrite the tresholding example using numpy

vol = aims.read(os.path.join(src,'data_for_anatomist/subject01/Audio-Video_T_map.nii'))
arr = numpy.asarray(vol)
arr[numpy.where(arr < 3.)] = 0.
# print(vol.value(20, 20, 20))
aims.write(vol, 'Audio-Video_T_thresholded2.nii')
assert(vol.value(20, 20, 20) == 0)

vol = aims.Volume(192, 256, 128, 'S16')
vol.header()['voxel_size'] = [0.9, 0.9, 1.2, 1.]
vs = vol.header()['voxel_size']
pos0 = (100 * vs[0], 100 * vs[1], 60 * vs[2]) # in millimeters
arr = numpy.asarray(vol)
# build arrays of coordinates for x, y, z
x, y, z = numpy.ogrid[0.:vol.getSizeX(), 0.:vol.getSizeY(), 0.:vol.getSizeZ()] # Is it equivalent to x = numpy.arange(0:vol.getSizeX()), and similarly for y and z ?
# get coords in millimeters
x *= vs[0]
y *= vs[1]
z *= vs[2]
# relative to pos0
x -= pos0[0]
y -= pos0[1]
z -= pos0[2]
# get norm, using numpy arrays broadcasting
arr[:, :, :, 0] = numpy.sqrt(x**2 + y**2 + z**2)

# print(vol.value(100, 100, 60))
assert(vol.value(100, 100, 60) == 0)

# and save result
aims.write(vol, 'distance2.nii')

## Deep-copy of a volume

vol2 = aims.Volume(vol)
vol2.setValue(12, 100, 100, 60)
# now vol and vol2 have different values
# print('vol.value(100, 100, 60):', vol.value(100, 100, 60))
assert(vol.value(100, 100, 60) == 0)
# print('vol2.value(100, 100, 60):', vol2.value(100, 100, 60))
assert(vol2.value(100, 100, 60) == 12)

vol2 = aims.Volume(vol.getSizeX(), vol.getSizeY(), vol.getSizeZ(), vol.getSizeT(), 'FLOAT')
vol2.header().update(vol.header())
# print(vol2.header())
assert(vol2.header() == {'sizeX': 192, 'sizeY': 256, 'sizeZ': 128, 'sizeT': 1, 'voxel_size': [0.9, 0.9, 1.2, 1]})

arr = numpy.array(numpy.diag(range(40)), dtype=numpy.float32).reshape(40, 40, 1) \
    + numpy.array(range(20), dtype=numpy.float32).reshape(1, 1, 20)
# WARNING: the array must be in Fortran ordering for AIMS, at least at the moment
# whereas the numpy addition always returns a C-ordered array
arr = numpy.array(arr, order='F')
arr[10, 12, 3] = 25
vol = aims.Volume(arr)
print('vol.value(10, 12, 3):', vol.value(10, 12, 3))
assert(vol.value(10, 12, 3) == 25.)

# data are shared with arr
vol.setValue(35, 10, 15, 2)
print('arr[10, 15, 2]:', arr[10, 15, 2])
assert(arr[10, 15, 2] == 35.0)
arr[12, 15, 1] = 44
print('vol.value(12, 15, 1):', vol.value(12, 15, 1))
assert(vol.value(12, 15, 1) == 44.0)

# create a 4D volume of signed 16-bit ints, of size 30x30x30x4
vol = aims.Volume(30, 30, 30, 4, 'S16')
# fill it with zeros
vol.fill(0)
# set value 12 at voxel (10, 10, 20, 2)
vol.setValue(12, 10, 10, 20, 2)
# get value at the same position
x = vol.value(10, 10, 20, 2)
print(x)
assert(x == 12)
# set the voxels size
vol.header()['voxel_size'] = [0.9, 0.9, 1.2, 1.]
print(vol.header())
assert(vol.header() == {'sizeX': 30, 'sizeY': 30, 'sizeZ': 30, 'sizeT': 4, 'voxel_size': [0.9, 0.9, 1.2, 1]})

mesh = aims.read(os.path.join(src,'data_for_anatomist/subject01/subject01_Lhemi.mesh'))
print(mesh.header().items())
vert = mesh.vertex()
print('vertices:', len(vert))
assert(len(vert) == 33837)
poly = mesh.polygon()
print('polygons:', len(poly))
assert(len(poly) == 67678)
norm = mesh.normal()
print('normals:', len(norm))
assert(len(norm) == 33837)

mesh = aims.AimsTimeSurface(3)
print(mesh.header().items())
# a mesh has a header
mesh.header()['toto'] = 'a message in the header'
mesh.header()['header item']='another message in the header'
print(mesh.header().items())

vert = mesh.vertex()
poly = mesh.polygon()
x = numpy.cos(numpy.ogrid[0.: 20] * numpy.pi / 10.) * 100
y = numpy.sin(numpy.ogrid[0.: 20] * numpy.pi / 10.) * 100
z = numpy.zeros(20)
c = numpy.vstack((x, y, z)).transpose()
vert.assign(numpy.vstack((numpy.array([(0., 0., -40.), (0., 0., 40.)]), c))) # ERROR \
# vector_POINT3DF.assign(): argument 1 has unexpected type 'numpy.ndarray'
pol = numpy.vstack((numpy.zeros(20, dtype=numpy.int32), numpy.ogrid[3: 23], numpy.ogrid[2: 22])).transpose()
pol[19, 1] = 2
pol2 = numpy.vstack((numpy.ogrid[2: 22], numpy.ogrid[3: 23], numpy.ones(20, dtype=numpy.int32))).transpose()
pol2[19, 1] = 2
poly.assign(numpy.vstack((pol, pol2)))
# write result
aims.write(mesh, 'saucer.gii')
# automatically calculate normals
mesh.updateNormals()

## Modifying a mesh
#slightly inflate a mesh

mesh = aims.read(os.path.join(src,'data_for_anatomist/subject01/subject01_Lwhite.mesh'))
vert = mesh.vertex() # ERROR vector_POINT3DF.assign(): argument 1 has unexpected type 'numpy.ndarray' (same as before)
varr = numpy.array(vert)
norm = numpy.array(mesh.normal())
varr += norm * 2 # push vertices 2mm away along normal
vert.assign(varr)
mesh.updateNormals()
aims.write(mesh, 'subject01_Lwhite_semiinflated.mesh')

#alternative without numpy

mesh = aims.read(os.path.join(src,'data_for_anatomist/subject01/subject01_Lwhite.mesh'))
vert = mesh.vertex()
norm = mesh.normal()
for v, n in zip(vert, norm):
    v += n * 2
mesh.updateNormals()
aims.write(mesh, 'subject01_Lwhite_semiinflated.mesh')

## In AIMS, meshes are actually time-indexed dictionaries of meshes.
## This way a deforming mesh can be stored in the same object. To copy a timestep to another, use the following:

mesh = aims.read(os.path.join(src,'data_for_anatomist/subject01/subject01_Lwhite.mesh'))
# mesh.vertex() is equivalent to mesh.vertex(0)
mesh.vertex(1).assign(mesh.vertex(0))
# same for normals and polygons
mesh.normal(1).assign(mesh.normal(0))
mesh.polygon(1).assign(mesh.polygon(0))
print('number of time steps:', mesh.size())
assert(mesh.size() == 2)

## Exercise: Make a deforming mesh that goes from the original mesh to 5mm away, by steps of 0.5 mm

mesh = aims.read(os.path.join(src,'data_for_anatomist/subject01/subject01_Lwhite.mesh'))
norm = mesh.normal()
for i in range (1,10):
    mesh.normal(i).assign(mesh.normal())
    mesh.polygon(i).assign(mesh.polygon())
    mesh.vertex(i).assign(mesh.vertex())
    for nv,n in zip(mesh.vertex(i),norm):
        nv += n*i*0.5
print('number of time steps:', mesh.size())
assert(mesh.size() == 10)
mesh.updateNormals()
aims.write(mesh,'subject01_Lwhite_semiinflated_time.mesh')

# Correction (with numpy)

mesh = aims.read(os.path.join(src,'data_for_anatomist/subject01/subject01_Lwhite.mesh'))
vert = mesh.vertex()
varr = numpy.array(vert)
norm = numpy.array(mesh.normal())
for i in range(1, 10):
    mesh.normal(i).assign(mesh.normal())
    mesh.polygon(i).assign(mesh.polygon())
    varr += norm * 0.5
    mesh.vertex(i).assign(varr) # TypeError: vector_POINT3DF.assign(): argument 1 has unexpected type 'numpy.ndarray'
print('number of time steps:', mesh.size())
assert(mesh.size() == 10)
mesh.updateNormals()
aims.write(mesh, 'subject01_Lwhite_semiinflated_time.mesh')
"""

## TEXTURE

tex = aims.TimeTexture('FLOAT')
t = tex[0] # time index, inserts on-the-fly
t.reserve(10) # pre-allocates memory
for i in range(10):
    t.append(i / 10.)
# print(tex.size())
assert(len(tex) == 1)
# print(tex[0].size())
assert(len(tex[0]) == 10)
# print(tex[0][5])
assert(tex[0][5] == 0.5)

# Make a time-texture, with at each time/vertex of the previous mesh,
# sets the value of the underlying volume data_for_anatomist/subject01/subject01.nii

mesh = aims.read('subject01_Lwhite_semiinflated_time.mesh')
vol = aims.read(os.path.join(src,'data_for_anatomist/subject01/subject01.nii'))
tex = aims.TimeTexture('FLOAT')
vs = vol.header()['voxel_size']
for i in range(mesh.size()):
    t = tex[i]
    vert = mesh.vertex(i)
    t.reserve(len(vert))
    for p in vert:
        t.append(vol.value(*[int(round(x / y)) for x, y in zip(p, vs)]))
aims.write(tex, 'subject01_Lwhite_semiinflated_texture.tex')


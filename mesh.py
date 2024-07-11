import argparse
import jsonref as json
import warnings 

import gmsh
import numpy as np

from auxiliaryFunctions import *

parser = argparse.ArgumentParser(
          prog='(py)gmsh based tool for the simpliest turbomachinery mesh',
          description='For more information see README.md.')

parser.add_argument('path', metavar='path', type=str,
                    help='path to json with description of the case')
args = parser.parse_args()

featuresFile = open(args.path)
f = json.load(featuresFile)

coor = dict()

geometry = f['geometry']
profile = geometry['profile']
mesh = f['mesh']
domain = f['domain']

coor_suction  = np.genfromtxt( f['working directory']+profile['name of suction' ], delimiter=profile['delimiter'] )
coor_pressure = np.genfromtxt( f['working directory']+profile['name of pressure'], delimiter=profile['delimiter'] )
if profile['scaled by']:
  coor_suction = coor_suction*profile['scaled by']
  coor_pressure = coor_pressure*profile['scaled by']

if np.size(coor_suction, axis=0) != np.size(coor_pressure, axis=0):
  warnings.warn( 'Numbers of points on pressure and suction side are not equal. It is implemened but the implementation is not tested. Might not work well.' )
  resol = np.max([np.size(coor_suction, axis=0), np.size(coor_pressure, axis=0)])
  x_interpol_0 = np.max( [coor_suction[:,0].min() , coor_pressure[:,0].min()] )
  x_interpol_l = np.min( [coor_suction[:,0].max() , coor_pressure[:,0].max()] )
  x_interpol = np.linspace(x_interpol_0, x_interpol_l, resol)

  for coor in [coor_suction, coor_pressure]:
    y_interpol = np.interp(x_interpol, coor[:,0], coor[:,1])
    coor = np.hstack( (np.array([x_interpol]).T, np.array([y_interpol]).T)  )

gmsh.initialize()
gmsh.model.add( f['name'] )
gmsh.model.mesh.setOrder(1)

numProfilePoints = np.size(coor_suction, axis=0)

points_pressure = np.zeros(numProfilePoints  , dtype=int)
points_suction  = np.zeros(numProfilePoints  , dtype=int)
points_midline  = np.zeros(numProfilePoints+2, dtype=int)

for i in range(numProfilePoints):

  points_pressure[i] = gmsh.model.occ.addPoint(coor_pressure[ i, 0 ], coor_pressure[ i, 1 ], 0)
  points_suction[i]  = gmsh.model.occ.addPoint( coor_suction[ i, 0 ],  coor_suction[ i, 1 ], 0)
  points_midline[i+1]  = gmsh.model.occ.addPoint(
    (coor_pressure[ i, 0 ]+coor_suction[ i, 0 ])/2, 
    (coor_pressure[ i, 1 ]+coor_suction[ i, 1 ])/2,
    0
    )

xLen_inlet, yLen_inlet = lengths(geometry['inlet flow angle'],domain['length of inlet'])

points_midline[0] =   gmsh.model.occ.addPoint(
    (coor_pressure[ 0, 0 ]+coor_suction[ 0, 0 ])/2 - xLen_inlet, 
    (coor_pressure[ 0, 1 ]+coor_suction[ 0, 1 ])/2 - yLen_inlet,
    0
    )

xLen_outlet, yLen_outlet = lengths(-geometry['outlet flow angle'],domain['length of outlet'])
points_midline[-1] =   gmsh.model.occ.addPoint(
    (coor_pressure[ -1, 0 ]+coor_suction[ -1, 0 ])/2 + xLen_outlet, 
    (coor_pressure[ -1, 1 ]+coor_suction[ -1, 1 ])/2 + yLen_outlet,
    0
    )

line_midline = gmsh.model.occ.add_bspline(points_midline)
line_suction = gmsh.model.occ.add_bspline(points_suction)
line_pressure = gmsh.model.occ.add_bspline(points_pressure)

yLen_pitch, xLen_pitch = lengths(geometry['stagger angle'], geometry['pitch']/2)

line_upperPeriodicity = gmsh.model.occ.copy([(1, line_midline)])
points_upperPeriodicity = gmsh.model.occ.copy([(0, points_midline[0]), (0, points_midline[-1])])
gmsh.model.occ.translate(line_upperPeriodicity, xLen_pitch, yLen_pitch, 0)
gmsh.model.occ.translate(points_upperPeriodicity, xLen_pitch, yLen_pitch, 0)

line_lowerPeriodicity = gmsh.model.occ.copy([(1, line_midline)])
points_lowerPeriodicity = gmsh.model.occ.copy([(0, points_midline[0]), (0, points_midline[-1])])
gmsh.model.occ.translate(line_lowerPeriodicity, -xLen_pitch, -yLen_pitch, 0)
gmsh.model.occ.translate(points_lowerPeriodicity, -xLen_pitch, -yLen_pitch, 0)

line_inlet = gmsh.model.occ.add_line(points_lowerPeriodicity[0][1], points_upperPeriodicity[0][1])
line_outlet = gmsh.model.occ.add_line(points_lowerPeriodicity[-1][1], points_upperPeriodicity[-1][1])

loop_outer = gmsh.model.occ.add_curve_loop([line_inlet, line_upperPeriodicity[0][1], -line_outlet, -line_lowerPeriodicity[0][1]])
loop_blade = gmsh.model.occ.add_curve_loop([line_pressure, -line_suction])

if mesh['boundary layer'] != 'transfinite':

  surface_whole = gmsh.model.occ.addPlaneSurface([loop_outer, loop_blade])

  volume = gmsh.model.occ.extrude([(2,surface_whole)], 0, 0, domain['thickness'], [mesh['n_layers_z']], recombine=True)

  gmsh.model.occ.synchronize()

  gmsh.model.addPhysicalGroup(2, [surface_whole], name = "back" )
  gmsh.model.addPhysicalGroup(2, [volume[0][1]], name = "front" )
  gmsh.model.addPhysicalGroup(2, [volume[2][1]], name = "inlet" )
  gmsh.model.addPhysicalGroup(2, [volume[4][1]], name = "outlet" )
  gmsh.model.addPhysicalGroup(2, [volume[3][1]], name = "periodicity_suction" )
  gmsh.model.addPhysicalGroup(2, [volume[5][1]], name = "periodicity_pressure" )
  gmsh.model.addPhysicalGroup(2, [volume[6][1]], name = "blade_pressure" )
  gmsh.model.addPhysicalGroup(2, [volume[7][1]], name = "blade_suction" )

  gmsh.model.addPhysicalGroup(3, [volume[1][1]], name = "fluid" )
  
  surface_periodicity_suction  = volume[3][1]
  surface_periodicity_pressure = volume[5][1]

else:
  blProperties = mesh['BL properties']
  coor_pressureInflated = inflateSide(coor_pressure, coor_suction, blProperties['thickness'], a =  1)
  coor_suctionInflated  = inflateSide(coor_suction, coor_pressure, blProperties['thickness'], a = -1)

  points_pressureInflated = np.zeros(numProfilePoints  , dtype=int)
  points_suctionInflated  = np.zeros(numProfilePoints  , dtype=int)
  
  for i in range(numProfilePoints):
    points_pressureInflated[i] = gmsh.model.occ.addPoint(coor_pressureInflated[ i, 0 ], coor_pressureInflated[ i, 1 ], 0)
    points_suctionInflated[i]  = gmsh.model.occ.addPoint( coor_suctionInflated[ i, 0 ],  coor_suctionInflated[ i, 1 ], 0)

  line_suctionInflated = gmsh.model.occ.add_bspline(  points_suctionInflated)
  line_pressureInflated = gmsh.model.occ.add_bspline(points_pressureInflated)
  line_midlineInflated_i = gmsh.model.occ.add_line(points_pressure[ 0], points_pressureInflated[ 0])
  line_midlineInflated_o = gmsh.model.occ.add_line(points_pressure[-1], points_pressureInflated[-1])

'''

if mesh['periodicities internal match']:
  translation = [1, 0, 0, 2*xLen_pitch, 
                 0, 1, 0, 2*yLen_pitch, 
                 0, 0, 1, 0, 
                 0, 0, 0, 1]
  gmsh.model.mesh.setPeriodic(2, [surface_periodicity_suction], [surface_periodicity_pressure], translation)

refinementFields  = []

if mesh['refine wake']:
  wake = mesh['wake']
  xLen_wake, yLen_wake = lengths(-geometry['outlet flow angle'],wake['length'])
  point_wake = gmsh.model.occ.addPoint(
    (coor_pressure[ -1, 0 ]+coor_suction[ -1, 0 ])/2 + xLen_wake, 
    (coor_pressure[ -1, 1 ]+coor_suction[ -1, 1 ])/2 + yLen_wake,
    0
    )

  line_wake = gmsh.model.occ.add_line(points_midline[-2], point_wake)
  gmsh.model.geo.synchronize()
  gmsh.model.occ.synchronize()

  gmsh.model.mesh.field.add("Distance", 1)
  gmsh.model.mesh.field.setNumbers(1, "CurvesList", [line_wake])
  gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
  
  gmsh.model.mesh.field.add("Threshold", 2)
  gmsh.model.mesh.field.setNumber(2, "InField", 1)
  gmsh.model.mesh.field.setNumber(2, "SizeMin", wake['size'])
  gmsh.model.mesh.field.setNumber(2, "SizeMax", mesh['max size'])
  gmsh.model.mesh.field.setNumber(2, "DistMin", wake['thickness'])
  gmsh.model.mesh.field.setNumber(2, "DistMax", wake['diffuse']*wake['thickness'])

  refinementFields = [2]

gmsh.model.mesh.field.add("Distance", 3)
gmsh.model.mesh.field.setNumbers(3, "CurvesList", [line_suction, line_pressure])
gmsh.model.mesh.field.setNumber(3, "Sampling", 100)
  
gmsh.model.mesh.field.add("Threshold", 4)
gmsh.model.mesh.field.setNumber(4, "InField", 3)
gmsh.model.mesh.field.setNumber(4, "SizeMin", mesh['baseline size'])
gmsh.model.mesh.field.setNumber(4, "SizeMax", mesh['max size'])
gmsh.model.mesh.field.setNumber(4, "DistMin", geometry['pitch']/4)
gmsh.model.mesh.field.setNumber(4, "DistMax", geometry['pitch']/2)

refinementFields = np.append(refinementFields, 4)

gmsh.model.mesh.field.add("Min", 7)
gmsh.model.mesh.field.setNumbers(7, "FieldsList", refinementFields)
gmsh.model.mesh.field.setAsBackgroundMesh(7)

if mesh['boundary layer'] == 'extruded':

  extrudedBL = gmsh.model.mesh.field.add('BoundaryLayer')
  gmsh.model.mesh.field.setNumbers(extrudedBL, 'CurvesList', [line_suction, line_pressure])
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Size', mesh['BL properties']['size']) 
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Ratio', mesh['BL properties']['ratio']) 
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Quads', 1)
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Thickness', mesh['BL properties']['thickness'])
  gmsh.option.setNumber('Mesh.BoundaryLayerFanElements', 0)
  gmsh.model.mesh.field.setAsBoundaryLayer(extrudedBL)

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", mesh['mesh size from curvature'])


gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh["max size"] )
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh["min size"]  )
gmsh.option.setNumber("Mesh.RecombineAll", 1)

if f['create mesh']: gmsh.model.mesh.generate(3)
'''
gmsh.model.geo.synchronize()
gmsh.model.occ.synchronize()
if f['save']:    gmsh.write(f['working directory']+f['name']+f['format'])
if f['run GUI']: gmsh.fltk.run()

import argparse
import jsonref as json
import warnings 

import gmsh
import numpy as np

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

points_pressure = np.zeros(numProfilePoints)
points_suction  = np.zeros(numProfilePoints)
points_midline  = np.zeros(numProfilePoints+2)

for i in range(numProfilePoints):

  points_pressure[i] = gmsh.model.occ.addPoint(coor_pressure[ i, 0 ], coor_pressure[ i, 1 ], 0)
  points_suction[i]  = gmsh.model.occ.addPoint( coor_suction[ i, 0 ],  coor_suction[ i, 1 ], 0)
  points_midline[i+1]  = gmsh.model.occ.addPoint(
    (coor_pressure[ i, 0 ]+coor_suction[ i, 0 ])/2, 
    (coor_pressure[ i, 1 ]+coor_suction[ i, 1 ])/2,
    0
    )

xLen_inlet = np.cos(np.deg2rad(geometry['inlet flow angle']))*domain['length of inlet']
yLen_inlet = np.sin(np.deg2rad(geometry['inlet flow angle']))*domain['length of inlet']
points_midline[0] =   gmsh.model.occ.addPoint(
    (coor_pressure[ 0, 0 ]+coor_suction[ 0, 0 ])/2 - xLen_inlet, 
    (coor_pressure[ 0, 1 ]+coor_suction[ 0, 1 ])/2 - yLen_inlet,
    0
    )

xLen_outlet = np.cos(np.deg2rad(-geometry['outlet flow angle']))*domain['length of outlet']
yLen_outlet = np.sin(np.deg2rad(-geometry['outlet flow angle']))*domain['length of outlet']
points_midline[-1] =   gmsh.model.occ.addPoint(
    (coor_pressure[ -1, 0 ]+coor_suction[ -1, 0 ])/2 + xLen_outlet, 
    (coor_pressure[ -1, 1 ]+coor_suction[ -1, 1 ])/2 + yLen_outlet,
    0
    )

line_midline = gmsh.model.occ.add_bspline(points_midline)
line_suction = gmsh.model.occ.add_bspline(points_suction)
line_pressure = gmsh.model.occ.add_bspline(points_pressure)

xLen_pitch = np.sin(np.deg2rad(geometry['stagger angle']))*geometry['pitch']/2
yLen_pitch = np.cos(np.deg2rad(geometry['stagger angle']))*geometry['pitch']/2

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

gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh["max size"] )
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh["min size"]  )
gmsh.option.setNumber("Mesh.RecombineAll", 1)

gmsh.model.geo.synchronize()
gmsh.model.occ.synchronize()
if f['save']:    gmsh.write(f['working directory']+f['name']+f['format'])
if f['run GUI']: gmsh.fltk.run()


import gmsh
import numpy as np

import auxiliaryFunctions as aux

f = aux.initialize()

geometry = f['geometry']
profile  = geometry['profile']
mesh     = f['mesh']
domain   = f['domain']

coor_suction, coor_pressure = aux.loadProfile(f)

gmsh.initialize()
gmsh.model.add( f['name'] )
gmsh.model.mesh.setOrder(1)

numProfilePoints = np.size(coor_suction, axis=0)

points_pressure = np.zeros(numProfilePoints  , dtype=int)
points_suction  = np.zeros(numProfilePoints  , dtype=int)

xLen_inlet, yLen_inlet = aux.lengths(
  geometry['inlet flow angle'],
  domain['length of inlet']
  )
xLen_outlet, yLen_outlet = aux.lengths(
  -geometry['outlet flow angle'],
  domain['length of outlet']
  )

midline_x = (coor_pressure[:, 0] + coor_suction[:, 0])/2
midline_y = (coor_pressure[:, 1] + coor_suction[:, 1])/2

midline_x = np.append(
  (coor_pressure[ 0, 0 ]+coor_suction[ 0, 0 ])/2 - xLen_inlet, 
  midline_x
)
midline_x = np.append(
  midline_x,
  (coor_pressure[ -1, 0 ]+coor_suction[ -1, 0 ])/2 + xLen_outlet
)

midline_y = np.append(
  (coor_pressure[ 0, 1 ]+coor_suction[ 0, 1 ])/2 - yLen_inlet, 
  midline_y
)
midline_y = np.append(
  midline_y,
  (coor_pressure[ -1, 1 ]+coor_suction[ -1, 1 ])/2 + yLen_outlet
)

if 'smooth midline' in mesh.keys():
  midline_xSmooth = np.linspace(midline_x[1], midline_x[-2], mesh['level of sharpness'])
  midline_xSmooth = np.append(
    (coor_pressure[ 0, 0 ]+coor_suction[ 0, 0 ])/2 - xLen_inlet, 
    midline_xSmooth
  )
  midline_xSmooth = np.append(
    midline_xSmooth,
    (coor_pressure[ -1, 0 ]+coor_suction[ -1, 0 ])/2 + xLen_outlet
  )

  midline_ySmooth = np.interp(
    midline_xSmooth.astype(float), 
    midline_x.astype(float), 
    midline_y.astype(float))

  midline_x = midline_xSmooth
  midline_y = midline_ySmooth

points_midline  = np.zeros(len(midline_x), dtype=int)

for i in range(numProfilePoints):

  points_pressure[i] = gmsh.model.occ.addPoint(
    coor_pressure[i, 0], 
    coor_pressure[i, 1],
    -domain['thickness']/2)
  points_suction[i]  = gmsh.model.occ.addPoint(
    coor_suction[i, 0], 
    coor_suction[i, 1],
    -domain['thickness']/2)
  
for i in range(len(midline_x)):
  points_midline[i]  = gmsh.model.occ.addPoint(
    midline_x[i],
    midline_y[i],
    -domain['thickness']/2
    )

for i in [0, -1]:
  points_suction[i] = points_pressure[i]

points_midline[1] = points_pressure[0]
points_midline[-2] = points_pressure[-1]

line_midline  = gmsh.model.occ.add_bspline(points_midline)
line_suction  = gmsh.model.occ.add_spline(points_suction)
line_pressure = gmsh.model.occ.add_spline(points_pressure)

yLen_pitch, xLen_pitch = aux.lengths(
  geometry['stagger angle']*(not 'stagger included in definition' in profile.keys()), 
  geometry['pitch']/2
  )

line_upperPeriodicity   = gmsh.model.occ.copy([
  (1, line_midline)
  ])
points_upperPeriodicity = gmsh.model.occ.copy([
  (0, points_midline[0]), 
  (0, points_midline[-1])
  ])

line_lowerPeriodicity = gmsh.model.occ.copy([
  (1, line_midline)
  ])
points_lowerPeriodicity = gmsh.model.occ.copy([
  (0, points_midline[0]), 
  (0, points_midline[-1])
  ])

centerOfProfile = [(coor_pressure[int(numProfilePoints/2), 0]+coor_suction[int(numProfilePoints/2), 0])/2,
                   (coor_pressure[int(numProfilePoints/2), 1]+coor_suction[int(numProfilePoints/2), 1])/2]

if 'dilate' in geometry.keys():
  gmsh.model.occ.dilate(  line_upperPeriodicity ,centerOfProfile[0], centerOfProfile[1], 0, 1,1-geometry['dilate'],1)
  gmsh.model.occ.dilate(points_upperPeriodicity ,centerOfProfile[0], centerOfProfile[1], 0, 1,1-geometry['dilate'],1)
  gmsh.model.occ.dilate(  line_lowerPeriodicity ,centerOfProfile[0], centerOfProfile[1], 0, 1,1-geometry['dilate'],1)
  gmsh.model.occ.dilate(points_lowerPeriodicity ,centerOfProfile[0], centerOfProfile[1], 0, 1,1-geometry['dilate'],1)

gmsh.model.occ.translate(  line_upperPeriodicity, xLen_pitch, yLen_pitch, 0)
gmsh.model.occ.translate(points_upperPeriodicity, xLen_pitch, yLen_pitch, 0)

gmsh.model.occ.translate(  line_lowerPeriodicity, -xLen_pitch, -yLen_pitch, 0)
gmsh.model.occ.translate(points_lowerPeriodicity, -xLen_pitch, -yLen_pitch, 0)

line_inlet  = gmsh.model.occ.add_line(
  points_lowerPeriodicity[0][1], 
  points_upperPeriodicity[0][1]
  )
line_outlet = gmsh.model.occ.add_line(
  points_lowerPeriodicity[-1][1],
  points_upperPeriodicity[-1][1]
  )

loop_outer = gmsh.model.occ.add_curve_loop([
  line_inlet, 
  line_upperPeriodicity[0][1], 
  -line_outlet, 
  -line_lowerPeriodicity[0][1]
  ])

loop_blade = gmsh.model.occ.add_curve_loop([line_pressure, -line_suction])

surface_whole = gmsh.model.occ.addPlaneSurface([loop_outer, loop_blade])

if '2D' not in f.keys():
  extrude_z, extrude_num, extrude_heights = aux.getExtrusionParameters(f)
  volume = gmsh.model.occ.extrude([(2,surface_whole)], 
                                    0, 0, extrude_z, 
                                    extrude_num, 
                                    extrude_heights,
                                    recombine=True)

  gmsh.model.occ.synchronize(), gmsh.model.geo.synchronize()

  gmsh.model.addPhysicalGroup(2, [surface_whole], name = "back" )
  gmsh.model.addPhysicalGroup(2, [volume[0][1]], name = "front" )
  gmsh.model.addPhysicalGroup(2, [volume[2][1]], name = "inlet" )
  gmsh.model.addPhysicalGroup(2, [volume[4][1]], name = "outlet" )
  gmsh.model.addPhysicalGroup(2, [volume[3][1]], name = "periodicity_suction" )
  gmsh.model.addPhysicalGroup(2, [volume[5][1]], name = "periodicity_pressure" )
  gmsh.model.addPhysicalGroup(2, [volume[6][1]], name = "blade_pressure" )
  gmsh.model.addPhysicalGroup(2, [volume[7][1]], name = "blade_suction" )

  gmsh.model.addPhysicalGroup(3, [volume[1][1]], name = "fluid" )

  if 'periodicities internal match' in mesh.keys():

    periodicities = [3,5]
    if 'switch periodicities' in mesh.keys():
      periodicities = [5,3]  
      
    surface_periodicity_suction  = volume[periodicities[0]][1]
    surface_periodicity_pressure = volume[periodicities[1]][1]

    translation = [1, 0, 0, 2*xLen_pitch*(-1)**('switch periodicities' in mesh.keys()), 
                  0, 1, 0, 2*yLen_pitch*(-1)**('switch periodicities' in mesh.keys()), 
                  0, 0, 1, 0, 
                  0, 0, 0, 1]
    
    gmsh.model.mesh.setPeriodic(2, [surface_periodicity_suction], 
                                  [surface_periodicity_pressure], translation)

else:

  gmsh.model.occ.synchronize(), gmsh.model.geo.synchronize()

  gmsh.model.addPhysicalGroup(1, [line_inlet] , name = "inlet" )
  gmsh.model.addPhysicalGroup(1, [line_outlet], name = "outlet" )
  gmsh.model.addPhysicalGroup(1, [line_upperPeriodicity[0][1]], name = "periodicity_suction" )
  gmsh.model.addPhysicalGroup(1, [line_lowerPeriodicity[0][1]], name = "periodicity_pressure" )
  gmsh.model.addPhysicalGroup(1, [line_pressure], name = "blade_pressure" )
  gmsh.model.addPhysicalGroup(1, [line_suction], name = "blade_suction" )

  gmsh.model.addPhysicalGroup(2, [surface_whole], name = "fluid" )

  if 'periodicities internal match' in mesh.keys():

    periodicities = [line_upperPeriodicity, line_lowerPeriodicity]
    if 'switch periodicities' in mesh.keys():
      periodicities = [line_lowerPeriodicity, line_upperPeriodicity] 
      
    surface_periodicity_suction  = periodicities[0]
    surface_periodicity_pressure = periodicities[1]

    translation = [1, 0, 0, 2*xLen_pitch*(-1)**('switch periodicities' in mesh.keys()), 
                  0, 1, 0, 2*yLen_pitch*(-1)**('switch periodicities' in mesh.keys()), 
                  0, 0, 1, 0, 
                  0, 0, 0, 1]
    
    gmsh.model.mesh.setPeriodic(1, [surface_periodicity_suction [0][1]], 
                                   [surface_periodicity_pressure[0][1]], translation)


refinementFields  = []

if 'refine wake' in mesh.keys():
  wake = mesh['wake']
  xLen_wake, yLen_wake = aux.lengths(
    -geometry['outlet flow angle'],
    wake['length']
    )
  point_wake = gmsh.model.occ.addPoint(
    (coor_pressure[ -1, 0 ]+coor_suction[ -1, 0 ])/2 + xLen_wake, 
    (coor_pressure[ -1, 1 ]+coor_suction[ -1, 1 ])/2 + yLen_wake,
    -domain['thickness']/2
    )

  line_wake = gmsh.model.occ.add_line(points_midline[-2], point_wake)

  gmsh.model.occ.synchronize(), gmsh.model.geo.synchronize()

  gmsh.model.mesh.field.add("Distance", 1)
  gmsh.model.mesh.field.setNumbers(1, "CurvesList", [line_wake])
  gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
  
  gmsh.model.mesh.field.add("Threshold", 2)
  gmsh.model.mesh.field.setNumber(2, "InField", 1)
  gmsh.model.mesh.field.setNumber(2, "SizeMin", wake['size'])
  gmsh.model.mesh.field.setNumber(2, "SizeMax", mesh['max size'])
  gmsh.model.mesh.field.setNumber(2, "DistMin", wake['thickness'])
  gmsh.model.mesh.field.setNumber(2, "DistMax", wake['diffuse']*wake['thickness'])

  refinementFields.append(2)

gmsh.model.mesh.field.add("Distance", 3)
gmsh.model.mesh.field.setNumbers(3, "CurvesList", [line_suction, line_pressure])
gmsh.model.mesh.field.setNumber(3, "Sampling", 100)
  
gmsh.model.mesh.field.add("Threshold", 4)
gmsh.model.mesh.field.setNumber(4, "InField", 3)
gmsh.model.mesh.field.setNumber(4, "SizeMin", mesh['baseline size'])
gmsh.model.mesh.field.setNumber(4, "SizeMax", mesh['max size'])
gmsh.model.mesh.field.setNumber(4, "DistMin", geometry['pitch']/4)
gmsh.model.mesh.field.setNumber(4, "DistMax", geometry['pitch'])

refinementFields.append(4)

if 'refine LE' in mesh.keys():
  gmsh.model.mesh.field.add("Distance", 5)
  gmsh.model.mesh.field.setNumbers(5, "PointsList", [points_pressure[0]])
    
  gmsh.model.mesh.field.add("Threshold", 6)
  gmsh.model.mesh.field.setNumber(6, "InField", 5)
  gmsh.model.mesh.field.setNumber(6, "SizeMin", mesh['LE']['size'])
  gmsh.model.mesh.field.setNumber(6, "SizeMax", mesh['max size'])
  gmsh.model.mesh.field.setNumber(6, "DistMin", mesh['LE']['radius'])
  gmsh.model.mesh.field.setNumber(6, "DistMax", mesh['LE']['radius']*mesh['LE']['diffuse'])

  refinementFields.append(6)

if 'refine TE' in mesh.keys():
  gmsh.model.mesh.field.add("Distance", 7)
  gmsh.model.mesh.field.setNumbers(7, "PointsList", [points_pressure[-1]])
    
  gmsh.model.mesh.field.add("Threshold", 8)
  gmsh.model.mesh.field.setNumber(8, "InField", 7)
  gmsh.model.mesh.field.setNumber(8, "SizeMin", mesh['TE']['size'])
  gmsh.model.mesh.field.setNumber(8, "SizeMax", mesh['max size'])
  gmsh.model.mesh.field.setNumber(8, "DistMin", mesh['TE']['radius'])
  gmsh.model.mesh.field.setNumber(8, "DistMax", mesh['TE']['radius']*mesh['TE']['diffuse'])

  refinementFields.append(8)

gmsh.model.mesh.field.add("Min", 100)
gmsh.model.mesh.field.setNumbers(100, "FieldsList", refinementFields)
gmsh.model.mesh.field.setAsBackgroundMesh(100)

if 'boundary layer' in mesh.keys():

  extrudedBL = gmsh.model.mesh.field.add('BoundaryLayer')
  gmsh.model.mesh.field.setNumbers(extrudedBL, 'CurvesList', [line_suction, line_pressure])
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Size', mesh['BL properties']['size']) 
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Ratio', mesh['BL properties']['ratio']) 
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Quads', 1)
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Thickness', mesh['BL properties']['thickness'])
  gmsh.model.mesh.field.setAsBoundaryLayer(extrudedBL)

if 'extend from boundary' in mesh.keys():
  gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
else: gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)

gmsh.option.setNumber("General.NumThreads", 4)

gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh["max size"] )
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh["min size"] )

gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.RecombineMinimumQuality",   0.3)
gmsh.option.setNumber("Mesh.RecombineOptimizeTopology",25)

gmsh.model.mesh.setCompound(1, [line_pressure, line_suction])

if 'recombination algorithm' in mesh.keys():
  gmsh.option.setNumber("Mesh.RecombinationAlgorithm", mesh["recombination algorithm"] )

if 'mesh algorithm' in mesh.keys():
  gmsh.option.setNumber("Mesh.Algorithm", mesh["mesh algorithm"] )

if 'mesh size from curvature' in mesh.keys():
  gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", mesh['mesh size from curvature'])

if 'create mesh' in f.keys(): gmsh.model.mesh.generate(2 if '2D' in f.keys() else 3)
if 'version'     in f.keys(): gmsh.option.setNumber("Mesh.MshFileVersion",f['version'])

gmsh.model.occ.synchronize(), gmsh.model.geo.synchronize()

if 'save' in f.keys(): gmsh.write(f['working directory']+f['name']+f['format'])
if 'run GUI' in f.keys(): gmsh.fltk.run()

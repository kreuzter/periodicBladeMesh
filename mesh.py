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
points_midline  = np.zeros(numProfilePoints+2, dtype=int)

for i in range(numProfilePoints):

  points_pressure[i] = gmsh.model.occ.addPoint(coor_pressure[i, 0], coor_pressure[i, 1], 0)
  points_suction[i]  = gmsh.model.occ.addPoint( coor_suction[i, 0],  coor_suction[i, 1], 0)
  points_midline[i+1]  = gmsh.model.occ.addPoint(
    (coor_pressure[ i, 0 ]+coor_suction[ i, 0 ])/2, 
    (coor_pressure[ i, 1 ]+coor_suction[ i, 1 ])/2,
    0
    )

for i in [0, -1]:
  points_suction[i] = points_pressure[i]

points_midline[1] = points_pressure[0]
points_midline[-2] = points_pressure[-1]
  
xLen_inlet, yLen_inlet = aux.lengths(
  geometry['inlet flow angle'],
  domain['length of inlet']
  )

points_midline[0] =   gmsh.model.occ.addPoint(
    (coor_pressure[ 0, 0 ]+coor_suction[ 0, 0 ])/2 - xLen_inlet, 
    (coor_pressure[ 0, 1 ]+coor_suction[ 0, 1 ])/2 - yLen_inlet,
    0
    )

xLen_outlet, yLen_outlet = aux.lengths(
  -geometry['outlet flow angle'],
  domain['length of outlet']
  )
points_midline[-1] =   gmsh.model.occ.addPoint(
    (coor_pressure[ -1, 0 ]+coor_suction[ -1, 0 ])/2 + xLen_outlet, 
    (coor_pressure[ -1, 1 ]+coor_suction[ -1, 1 ])/2 + yLen_outlet,
    0
    )

line_midline = gmsh.model.occ.add_bspline(points_midline)
line_suction = gmsh.model.occ.add_bspline(points_suction)
line_pressure = gmsh.model.occ.add_bspline(points_pressure)

yLen_pitch, xLen_pitch = aux.lengths(
  geometry['stagger angle']*(not profile['stagger included in definition']), 
  geometry['pitch']/2
  )

line_upperPeriodicity   = gmsh.model.occ.copy([
  (1, line_midline)
  ])
points_upperPeriodicity = gmsh.model.occ.copy([
  (0, points_midline[0]), 
  (0, points_midline[-1])
  ])

gmsh.model.occ.translate(line_upperPeriodicity, xLen_pitch, yLen_pitch, 0)
gmsh.model.occ.translate(points_upperPeriodicity, xLen_pitch, yLen_pitch, 0)

line_lowerPeriodicity = gmsh.model.occ.copy([
  (1, line_midline)
  ])
points_lowerPeriodicity = gmsh.model.occ.copy([
  (0, points_midline[0]), 
  (0, points_midline[-1])
  ])

gmsh.model.occ.translate(line_lowerPeriodicity, -xLen_pitch, -yLen_pitch, 0)
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

#gmsh.option.setNumber("Mesh.Algorithm", 8) 

if mesh['boundary layer'] != 'transfinite':

  loop_blade = gmsh.model.occ.add_curve_loop([line_pressure, -line_suction])

  surface_whole = gmsh.model.occ.addPlaneSurface([loop_outer, loop_blade])

  volume = gmsh.model.occ.extrude([(2,surface_whole)], 
                                  0, 0, domain['thickness'], 
                                  [mesh['n_layers_z']], 
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
  
  surface_periodicity_suction  = volume[3][1]
  surface_periodicity_pressure = volume[5][1]

else:
  blProperties = mesh['BL properties']
  coor_pressureInflated = aux.inflateSide(coor_pressure, coor_suction, 
                                          blProperties['thickness'], a =  1)
  coor_suctionInflated  = aux.inflateSide(coor_suction, coor_pressure, 
                                          blProperties['thickness'], a = -1)

  points_pressureInflated = np.zeros(numProfilePoints  , dtype=int)
  points_suctionInflated  = np.zeros(numProfilePoints  , dtype=int)
  
  for i in range(numProfilePoints):
    points_pressureInflated[i] = gmsh.model.occ.addPoint(
      coor_pressureInflated[ i, 0 ], 
      coor_pressureInflated[ i, 1 ], 
      0
      )
    points_suctionInflated[i]  = gmsh.model.occ.addPoint( 
      coor_suctionInflated[ i, 0 ],  
      coor_suctionInflated[ i, 1 ], 
      0
      )

  for i in [0,-1]:
    points_suctionInflated[i] = points_pressureInflated[i]

  line_suctionInflated  = gmsh.model.occ.add_bspline(  points_suctionInflated)
  line_pressureInflated = gmsh.model.occ.add_bspline(points_pressureInflated)
  line_midlineInflated_i = gmsh.model.occ.add_line(
    points_pressure[ 0], 
    points_pressureInflated[ 0]
    )
  line_midlineInflated_o = gmsh.model.occ.add_line(
    points_pressure[-1], 
    points_pressureInflated[-1]
    )

  loop_inflated = gmsh.model.occ.add_curve_loop([
    line_pressureInflated, 
    -line_suctionInflated
    ])
  loop_pressureBlade = gmsh.model.occ.add_curve_loop([
    line_pressure, 
    line_midlineInflated_o, 
    -line_pressureInflated, 
    -line_midlineInflated_i
    ])
  loop_suctionBlade = gmsh.model.occ.add_curve_loop([
    line_suction, 
    line_midlineInflated_o, 
    -line_suctionInflated, 
    -line_midlineInflated_i
    ])
  
  surface_pressure = gmsh.model.occ.addPlaneSurface([loop_pressureBlade])
  surface_suction  = gmsh.model.occ.addPlaneSurface([loop_suctionBlade])
  surface_outer    = gmsh.model.occ.addPlaneSurface([loop_outer, loop_inflated])

  gmsh.model.occ.synchronize(), gmsh.model.geo.synchronize()

  for curve in [line_pressure, line_suction]:
    gmsh.model.mesh.setTransfiniteCurve(curve, blProperties['cells on blade'], 
                                        'Bump', 
                                        coef=1/blProperties['refinement of edges'])

  for curve in [line_pressureInflated, line_suctionInflated]:
    gmsh.model.mesh.setTransfiniteCurve(curve, blProperties['cells on blade'], 
                                        'Bump', 
                                        coef=1/blProperties['refinement of inflated edges'])

  for curve in [line_midlineInflated_i, line_midlineInflated_o]:
    gmsh.model.mesh.setTransfiniteCurve(curve, blProperties['num points'], 
                                        coef=blProperties['ratio'])

  gmsh.model.occ.synchronize(), gmsh.model.geo.synchronize()

  for surface in [surface_pressure, surface_suction]:
    gmsh.model.mesh.setTransfiniteSurface(surface, cornerTags=[points_pressure[0], 
                                                               points_pressure[-1],
                                                               points_pressureInflated[-1], 
                                                               points_pressureInflated[0]])
    gmsh.model.mesh.setRecombine(2, surface)
    gmsh.model.mesh.setSmoothing(2, surface, 5)
    
  gmsh.model.occ.synchronize(), gmsh.model.geo.synchronize()
  volume = gmsh.model.occ.extrude([(2,surf) for surf in [surface_outer, 
                                                         surface_pressure, 
                                                         surface_suction] ], 
                                  0, 0, domain['thickness'], 
                                  [mesh['n_layers_z']], recombine=True)

  gmsh.model.occ.synchronize(), gmsh.model.geo.synchronize()

  gmsh.model.addPhysicalGroup(2, [surface_outer, 
                                  surface_pressure, 
                                  surface_suction], name = "back" )
  gmsh.model.addPhysicalGroup(2, [volume[8][1], 
                                  volume[14][1], 
                                  volume[0][1]], name = "front" )
  gmsh.model.addPhysicalGroup(2, [volume[2][1]], name = "inlet" )
  gmsh.model.addPhysicalGroup(2, [volume[4][1]], name = "outlet" )
  gmsh.model.addPhysicalGroup(2, [volume[3][1]], name = "periodicity_suction" )
  gmsh.model.addPhysicalGroup(2, [volume[5][1]], name = "periodicity_pressure" )
  gmsh.model.addPhysicalGroup(2, [volume[6][1]], name = "blade_pressure" )
  gmsh.model.addPhysicalGroup(2, [volume[16][1]], name = "blade_suction" )

  gmsh.model.addPhysicalGroup(3, [ent[1] for ent in  gmsh.model.getEntities(3)], name = "fluid" )

  surface_periodicity_suction  = volume[3][1]
  surface_periodicity_pressure = volume[5][1]

if mesh['periodicities internal match']:
  translation = [1, 0, 0, 2*xLen_pitch, 
                 0, 1, 0, 2*yLen_pitch, 
                 0, 0, 1, 0, 
                 0, 0, 0, 1]
  gmsh.model.mesh.setPeriodic(2, [surface_periodicity_suction], 
                                 [surface_periodicity_pressure], translation)

refinementFields  = []

if mesh['refine wake']:
  wake = mesh['wake']
  xLen_wake, yLen_wake = aux.lengths(
    -geometry['outlet flow angle'],
    wake['length']
    )
  point_wake = gmsh.model.occ.addPoint(
    (coor_pressure[ -1, 0 ]+coor_suction[ -1, 0 ])/2 + xLen_wake, 
    (coor_pressure[ -1, 1 ]+coor_suction[ -1, 1 ])/2 + yLen_wake,
    0
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

  refinementFields = [2]

gmsh.model.mesh.field.add("Distance", 3)
gmsh.model.mesh.field.setNumbers(3, "CurvesList", [line_suction, line_pressure])
gmsh.model.mesh.field.setNumber(3, "Sampling", 100)
  
gmsh.model.mesh.field.add("Threshold", 4)
gmsh.model.mesh.field.setNumber(4, "InField", 3)
gmsh.model.mesh.field.setNumber(4, "SizeMin", mesh['baseline size'])
gmsh.model.mesh.field.setNumber(4, "SizeMax", mesh['max size'])
gmsh.model.mesh.field.setNumber(4, "DistMin", geometry['pitch']/4)
gmsh.model.mesh.field.setNumber(4, "DistMax", geometry['pitch'])

refinementFields = np.append(refinementFields, 4)

if mesh['refine LE']:
  gmsh.model.mesh.field.add("Distance", 5)
  gmsh.model.mesh.field.setNumbers(5, "PointsList", [points_pressure[0]])
    
  gmsh.model.mesh.field.add("Threshold", 6)
  gmsh.model.mesh.field.setNumber(6, "InField", 5)
  gmsh.model.mesh.field.setNumber(6, "SizeMin", mesh['LE']['size'])
  gmsh.model.mesh.field.setNumber(6, "SizeMax", mesh['max size'])
  gmsh.model.mesh.field.setNumber(6, "DistMin", mesh['LE']['radius'])
  gmsh.model.mesh.field.setNumber(6, "DistMax", mesh['LE']['radius']*mesh['LE']['diffuse'])

  refinementFields = np.append(refinementFields, 6)

if mesh['refine TE']:
  gmsh.model.mesh.field.add("Distance", 7)
  gmsh.model.mesh.field.setNumbers(7, "PointsList", [points_pressure[-1]])
    
  gmsh.model.mesh.field.add("Threshold", 8)
  gmsh.model.mesh.field.setNumber(8, "InField", 7)
  gmsh.model.mesh.field.setNumber(8, "SizeMin", mesh['TE']['size'])
  gmsh.model.mesh.field.setNumber(8, "SizeMax", mesh['max size'])
  gmsh.model.mesh.field.setNumber(8, "DistMin", mesh['TE']['radius'])
  gmsh.model.mesh.field.setNumber(8, "DistMax", mesh['TE']['radius']*mesh['TE']['diffuse'])

  refinementFields = np.append(refinementFields, 8)

gmsh.model.mesh.field.add("Min", 100)
gmsh.model.mesh.field.setNumbers(100, "FieldsList", refinementFields)
gmsh.model.mesh.field.setAsBackgroundMesh(100)

if mesh['boundary layer'] == 'extruded':

  extrudedBL = gmsh.model.mesh.field.add('BoundaryLayer')
  gmsh.model.mesh.field.setNumbers(extrudedBL, 'CurvesList', [line_suction, line_pressure])
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Size', mesh['BL properties']['size']) 
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Ratio', mesh['BL properties']['ratio']) 
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Quads', 1)
  gmsh.model.mesh.field.setNumber(extrudedBL, 'Thickness', mesh['BL properties']['thickness'])
  gmsh.model.mesh.field.setAsBoundaryLayer(extrudedBL)

gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", mesh['mesh size from curvature'])

gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh["max size"] )
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh["min size"]  )
gmsh.option.setNumber("Mesh.RecombineAll", 1)

if f['create mesh']: gmsh.model.mesh.generate(3)
if f['version']:     gmsh.option.setNumber("Mesh.MshFileVersion",f['version'])   
gmsh.model.occ.synchronize(), gmsh.model.geo.synchronize()

if f['save']:    gmsh.write(f['working directory']+f['name']+f['format'])
if f['run GUI']: gmsh.fltk.run()

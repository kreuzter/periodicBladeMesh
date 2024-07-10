import argparse
import json
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
profile = f['profile']

coor_suction  = np.genfromtxt( f['working directory']+profile['name of suction' ], delimiter=profile['delimiter'] )
coor_pressure = np.genfromtxt( f['working directory']+profile['name of pressure'], delimiter=profile['delimiter'] )
if profile['scaled']:
  coor_suction = coor_suction*geometry[profile['scaled by']]
  coor_pressure = coor_pressure*geometry[profile['scaled by']]

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

if f['save']:
  gmsh.write(f['working directory']+f['name']+f['format'])

if f['run GUI']:
  gmsh.fltk.run()


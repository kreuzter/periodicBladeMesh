import numpy as np

import warnings 
import argparse
import jsonref as json

def initialize():
  parser = argparse.ArgumentParser(
          prog='(py)gmsh based tool for the simpliest turbomachinery mesh',
          description='For more information see README.md.')

  parser.add_argument('path', metavar='path', type=str,
                      help='path to json with description of the case')
  args = parser.parse_args()

  featuresFile = open(args.path)
  f = json.load(featuresFile)
  
  return f

def loadProfile(features):
  profile = features['geometry']['profile']
  coor_suction  = np.genfromtxt(features['working directory']+profile['name of suction' ], 
                                delimiter=profile['delimiter'])
  coor_pressure = np.genfromtxt(features['working directory']+profile['name of pressure'], 
                                delimiter=profile['delimiter'])
  
  if'reverse suction' in profile.keys(): 
    coor_suction = np.flip(coor_suction, axis=0)
  if'reverse pressure' in profile.keys(): 
    coor_pressure = np.flip(coor_pressure, axis=0)
  if 'scaled by' in profile:
    coor_suction = coor_suction*profile['scaled by']
    coor_pressure = coor_pressure*profile['scaled by']

  if np.size(coor_suction, axis=0) != np.size(coor_pressure, axis=0):
    coor_suction, coor_pressure = interpolate(coor_suction, coor_pressure)
  
  return coor_suction, coor_pressure

def lengths(angle, length): 
  return np.cos(np.deg2rad(angle))*length, np.sin(np.deg2rad(angle))*length

def interpolate(coor_suction, coor_pressure):
  warnings.warn("""Numbers of points on pressure and suction side are not equal. 
        It is implemened but the implementation is not tested. Might not work well.""")
  
  resol = np.max([np.size(coor_suction, axis=0), np.size(coor_pressure, axis=0)])
  x_interpol_0 = np.max([coor_suction[:,0].min(), coor_pressure[:,0].min()])
  x_interpol_l = np.min([coor_suction[:,0].max(), coor_pressure[:,0].max()])
  x_interpol = np.linspace(x_interpol_0, x_interpol_l, resol)

  y_interpol = np.interp(x_interpol, coor_suction[:,0], coor_suction[:,1])
  coor_suction = np.hstack((np.array([x_interpol]).T, np.array([y_interpol]).T))

  y_interpol = np.interp(x_interpol, coor_pressure[:,0], coor_pressure[:,1])
  coor_pressure = np.hstack((np.array([x_interpol]).T, np.array([y_interpol]).T))

  return coor_suction, coor_pressure

def elementSize(f):
  domain = f['domain']
  mesh = f['mesh']
  iterator = 0
  arr = [mesh["SW BL properties"]["size first"]]

  while np.sum(arr) <= domain['thickness']/2:
    arr = np.append(
      arr, 
      np.min([
        mesh["SW BL properties"]["size first"]*mesh["SW BL properties"]['ratio']**iterator, 
        mesh["SW BL properties"]["size last"]
        ])
      )
    iterator +=1
  
  return arr

def getExtrusionParameters(f):
  domain = f['domain']
  mesh = f['mesh']

  extrude_z = domain['thickness']

  if mesh['n layers in z'] == 1:
    extrude_num = [1]
    extrude_heights = []
  
  else:
    if 'side walls boundary layer' not in mesh.keys():
      extrude_num = [mesh['n layers in z']]
      extrude_heights = []
  
    else:
      extrude_size = elementSize(f)
      extrude_num = np.ones(len(extrude_size)*2)
      print(f'Number of cells in z is {len(extrude_size)*2}.')
      extrude_size = np.append(extrude_size, np.flip(extrude_size))
      extrude_heights = np.ones(len(extrude_size))
      extrude_heights[0] = mesh["SW BL properties"]["size first"]
      for i in range(1, len(extrude_heights)):
        extrude_heights[i] = extrude_heights[i-1]+extrude_size[i]
      extrude_heights = extrude_heights/extrude_heights[-1]

  return extrude_z, extrude_num, extrude_heights

if __name__ == "__main__":
  print('I do nothing, I am just a storage of functions.')
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
  if profile['reverse suction']: coor_suction = np.flip(coor_suction, axis=0)
  if profile['reverse pressure']: coor_pressure = np.flip(coor_pressure, axis=0)
  if profile['scaled by']:
    coor_suction = coor_suction*profile['scaled by']
    coor_pressure = coor_pressure*profile['scaled by']

  if np.size(coor_suction, axis=0) != np.size(coor_pressure, axis=0):
    coor_suction, coor_pressure = interpolate(coor_suction, coor_pressure)
  
  return coor_suction, coor_pressure

def lengths(angle, length): 
  return np.cos(np.deg2rad(angle))*length, np.sin(np.deg2rad(angle))*length

def inflateSide(sideToInflate, theOtherSide, thickness, a):
  numPoints = np.size(sideToInflate, axis=0)

  assert (numPoints == np.size(theOtherSide, axis=0))
  for idx in [0,-1]: assert (sideToInflate[idx,:] == theOtherSide[idx,:]).all()

  longerSuction = np.vstack((theOtherSide[2,:], sideToInflate))
  longerSuction = np.vstack((longerSuction, theOtherSide[-2,:]))

  normal = np.zeros(2)
  sideToInflateInflated  = np.zeros((numPoints, 2))

  for i in range(numPoints):
    tangent = longerSuction[i+2, :] - longerSuction[i, :]
    tangent = tangent/np.linalg.norm(tangent)
    normal = np.array([tangent[1], -tangent[0]])

    sideToInflateInflated[i, :] = sideToInflate[i, :]+a*normal*thickness   

  return sideToInflateInflated

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

def getExtrusionParameters(f):
  domain = f['domain']
  mesh = f['mesh']

  extrude_z = domain['thickness']

  if mesh['n layers in z'] == 1:
    extrude_num = [1]
    extrude_heights = []
  
  else:
    if not mesh['side walls boundary layer']:
      extrude_num = [mesh['n layers in z']]
      extrude_heights = []

    else:
      extrude_num = np.ones(mesh['n layers in z'])
  
      extrude_field = [
        mesh["SW BL properties"]["size first"], 
        mesh["SW BL properties"]["size first"]*(1+mesh["SW BL properties"]["ratio"])
        ]
  
      while np.abs(extrude_field[-1] - extrude_field[-2]) < mesh["SW BL properties"]["size last"]:
        extrude_field = np.append(
          extrude_field, 
          extrude_field[-1]+np.abs(extrude_field[-1]-extrude_field[-2])*mesh["SW BL properties"]["ratio"]
          )
      if len(extrude_field) > mesh['n layers in z']/2:
        extrude_field = np.delete(
          extrude_field, 
          range(int(len(extrude_field) - mesh['n layers in z']/2), len(extrude_field))
          )
      else:
        while len(extrude_field) < mesh['n layers in z']/2:
          extrude_field = np.append(
            extrude_field,
            extrude_field[-1] + mesh["SW BL properties"]["size last"])
          
      extrude_field2 = np.zeros(len(extrude_field))
      extrude_field2[0] = extrude_field[-1] + np.abs(extrude_field[-1] - extrude_field[-2] )
      for i in range(1,len(extrude_field)-1):
        extrude_field2[i] = extrude_field2[i-1] + np.abs(extrude_field[-(i+1)] - extrude_field[-(i+2)])
      extrude_field2[-1] = extrude_field2[-2] + mesh["SW BL properties"]["size first"]
      extrude_field3 = np.append(extrude_field, extrude_field2)
      
      print(extrude_field3[-1])
      extrude_heights = [extrude_field3[i]/extrude_field3[-1] for i in range(len(extrude_field3))]
  return extrude_z, extrude_num, extrude_heights

if __name__ == "__main__":
  print('I do nothing, I am just a storage of functions.')
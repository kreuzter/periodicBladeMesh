import numpy as np

def lengths(angle, length): return np.cos(np.deg2rad(angle))*length, np.sin(np.deg2rad(angle))*length

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

if __name__ == "__main__":
  print('I do nothing, I am just a storage of functions.')
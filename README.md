# periodicBladeMesh
(py)gmsh based tool for the simpliest turbomachinery mesh

All parameters of the geometry and mesh can be defined in json file, theoretically, there is no need to touch the python code. Practically, the code is definitely missing something, but it might be a good starting point.

Internal references can be used in json files. E.g.:
```
  "scaled by"  : {"$ref": "#/geometry/true chord"},
```

## Predefined refinements 

- `baseline size` is applied on the surface of the blade
- if allowed, leading and trailing edges are (independently) refined to `(L/T)E : size` with a given `(L/T)E : radius`
- if allowed, wake region - rectangle aligned with outlet flow with given `wake : length` and `wake : thickness` - is refined to `wake : size`

## Boundary layer meshing

Two setups of boundary layer are available:
  - extruded BL (preferred, gmsh internal tool)
  - transfinite mesh around blade (useless without manual tuning)

Inputs for **extruded** boundary layer mesh are:
- size of first layer, 
- growth ratio and 
- total thickness.

For the **transfinite** one:
- number of cells, 
- growth ratio and 
- total thickness

To translate between number of cells and size of first layer [online calculators](https://caefn.com/calculator/boundary-layer-mesh) can be used.

## Test cases:

### [SPLEEN](https://doi.org/10.5281/zenodo.7264761) turbine blade

run with 
```
  python3 mesh.py data/spleen.json
```

# periodicBladeMesh
(py)gmsh based tool for the simpliest turbomachinery mesh

Internal references can be used in json files. E.g.:
```
  "scaled by"  : {"$ref": "#/geometry/true chord"},
```

## Boundary layer meshing

Two setups of boundary layer available:
  - extruded BL (gmsh internal tool)
  - transfinite mesh around blade

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

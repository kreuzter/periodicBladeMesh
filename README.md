# periodicBladeMesh
(py)gmsh based tool for the simpliest turbomachinery mesh

Internal references can be used in json files. E.g.:
```
  "scaled by"  : {"$ref": "#/geometry/true chord"},
```

Two setups of boundary layer available:
  - extruded BL (gmsh internal tool)
  - transfinite mesh around blade (theoretically more stable and versatile)

## Test cases:

### [SPLEEN](https://doi.org/10.5281/zenodo.7264761) turbine blade

- run with 
```
  python3 mesh.py data/spleen.json
```
- LE and TE definition included in the suction and pressure side
- transfinite mesh around blade
- periodicities matching in geometry and mesh (the solver does not need to perform any kind of interpolation)

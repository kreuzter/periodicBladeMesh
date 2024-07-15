# periodicBladeMesh
(py)gmsh based tool for the simpliest turbomachinery mesh

All parameters of the geometry and mesh can be defined in json file, theoretically, there
is no need to touch the python code. Practically, the code is definitely missing something, 
but it might be a good starting point.

Internal references can be used in json files. E.g.:
```
 "scaled by"  : {"$ref": "#/geometry/true chord"},
```

json files cannot handle mathematics. Following will **not** work:

<strike>

```
 "length of inlet"  : 3*{"$ref": "#/geometry/true chord"}, 
```

</strike>

## Refinements 

- `baseline size` is applied on the **surface of the blade**
- `mesh size from curvature` is quite a tricky parameter. It could take care of the edges 
and there would be no need for tuning of parameters of their refinement.
But this parameter also refines the periodicities, which can lead to lower quality of mesh.
Use it carefully.
- if more refinements are affecting a region (e.g. on trailing edge refinement of itself, 
blade surface and wake take place), minimum required cell size is taken into account

### Predefined refinements

- if allowed, **leading and trailing edges** are (independently) refined to `(L/T)E : size` 
with a given `(L/T)E : radius`
- if allowed, **wake region** - rectangle aligned with outlet flow with given 
`wake : length` and `wake : thickness` - is refined to `wake : size`
- all the predefined refinements have a `diffuse` parameter, it sets a distance in which 
the mesh is not affected by given refinement

Example of refinement setup in a json file:

```
"refine TE": true,   | region around trailing edge will be refined 
"TE" : {
  "size"   : 1e-4,   | cells in the closest region will have size 1e-4 m
  "radius" : 1e-3,   | the closest region is a circle with radius 1e-3 m
  "diffuse": 2       | the transition between given size and rest of the 
                     | mesh will happen in a circle with radius 
                     | (1+diffuse)*radius = 3e-3 m
},
```

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

To translate between number of cells and size of first layer 
[online calculators](https://caefn.com/calculator/boundary-layer-mesh) can be used.

## Test cases:

### [SPLEEN](https://doi.org/10.5281/zenodo.7264761) turbine blade

run with 
```
  python3 mesh.py data/spleen.json
```

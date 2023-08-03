kanon
======

Program to compute chirality indices and assess molecular symmetry 

## Installing

### Building from source

To build kanon from source code you need to have:
- a Fortran compiler
- build systems: make or fpm
- [fpm](https://github.com/fortran-lang/fpm) version 0.2.0 or newer

#### Building with make

```
cp app/kanon.f90 src/
cp Makefile src/
cd src/
make clean
make
make clean
```

#### Building with fpm

```
fpm build
```

#### How to use kanon?

kanon uses a XYZFile as input

```
natom
comment
label1 xcoor1 ycoord1 zcoord1
label2 xcoor2 ycoord2 zcoord2
...
```

Examples
```
3
Water H2O.xyz
O     0.000000     0.000000     0.000000
H     0.000000     0.000000     0.947000
H     0.895670     0.000000    -0.316663

```

```
kanon --pattern H2O.xyz

---------------------------------------------------------------------------------
Cartesian coordinates input file: H2O.xyz
Number atoms:     3
 O    0.00000   0.00000   0.00000
 H    0.00000   0.00000   0.94700
 H    0.89567   0.00000  -0.31666

-----------------------------------------------------------------------------
Diameter of molecule:    1.54889 ang.
--------------------------------------------------------------------------------
Center of masses
   0.05011   0.00000   0.03527
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
Principal moments of inertia
     Ia           Ib           Ic
   0.53688      1.20907      1.74595
--------------------------------------------------------------------------------
Inertia axis
   0.00000   0.81430  -0.58045
   1.00000   0.00000   0.00000
   0.00000   0.58045   0.81430
--------------------------------------------------------------------------------
Symmetry elements. Tolerance:    0.05000
--------------------------------------------------------------------------------

```
```
kanon --single <H2O>.xyz

-----------------------------------------------------------------------------
Diameter of molecule:    1.54889 ang.
-----------------------------------------------------------------------------
Center of masses:    0.05011   0.00000   0.03527
------------------------------------------------------------------------------
Moments of inertia:       1.74595      1.20907      0.53688
Eigenvectors
     0.000     0.814    -0.580
     1.000     0.000     0.000
     0.000     0.580     0.814
---------------------------------------------------------------------------------
Symmetrical matrix in Diamond`s method to optimal rotation: B
---------------------------------------------------------------------------------
     0.000     0.000     0.000     0.000
     0.000     3.353    -0.000    -0.000
     0.000    -0.000     2.399    -0.001
     0.000    -0.000    -0.001     0.954
---------------------------------------------------------------------------------
 Eigenvalues   -   Eigenvectors
---------------------------------------------------------------------------------
     0.954     0.000     0.000     0.000     1.000
     0.000     1.000     0.000     0.000     0.000
     3.353     0.000     1.000    -0.000    -0.000
     2.399     0.000     0.000     1.000    -0.000
--------------------------------------------------------------------------------
Check files: Pattern.xyz / Rotate_image.xyz
---------------------------------------------------------------------------------
Root-mean-square deviation
RMSD:   0.00000000
--------------------------------------------------------------------------------
Hausdorff-derived chirality Index
CHI:    0.00000000
--------------------------------------------------------------------------------
```

Comparing an object with any image

```
kanon --pattern <object>.xyz --image <image>.xyz

```

Delaunay triangulation of surfaces

```
kanon --pattern <file>.xyz --delaunay <output_grid>.dat

```

```
output_grid.dat

GRID       1286
REF  6 2 1
            6         2         1
            6         3         2
            4        70         3
            4        69        70
           70        69        74
           71        70        74
...

```


 

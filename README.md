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
make
make clean
```

#### Building with fpm

```
fpm build
```

# Instructions
FFN MHD solvers

1. Where the files go
Inside your foamForNuclear clone, alongside the compressibleInterFoam module:
```
foamForNuclear/applications/modules/openFoamImported/
├── compressibleInterFoam/
│   ├── compressibleInterFoam.C
│   ├── compressibleInterFoam.H
│   └── ...
└── mhdFoam/                      ← new
    ├── mhdFoam.C
    ├── mhdFoam.H
    ├── UEqn.H
    ├── pEqn.H
    ├── BEqn.H
    ├── pBEqn.H
    └── Make/
        ├── files
        └── options
```
2. Environment :

```
source /path/to/OpenFOAM-vXXXX/etc/bashrc
# and, if foamForNuclear has its own environment script:
source /path/to/foamForNuclear/etc/bashrc 
```
3. Compiling :

From inside the module directory:
```
cd /path/to/foamForNuclear/applications/modules/openFoamImported/mhdFoam
wclean      # optional, clears any stale object files
wmake libso
```

`wmake libso` (rather than plain wmake) because the target is a shared library (LIB = ...), not an executable.
If there's a top-level Allwmake script in foamForNuclear (very likely, since GeN-Foam-style repos almost always have one), it's usually safer to use that instead, since it builds things in the right dependency order:
cd /path/to/foamForNuclear
./Allwmake

Sanity check : The file should exist `ls -la $FOAM_USER_LIBBIN/libMhdSolver.so`
Also, check if solver compiling was correct : 
```
nm -D $FOAM_USER_LIBBIN/libMhdSolver.so | grep -i mhd
```

4. Update GeN-Foam's Make/options to link the new lib
Add the include path back (this one is needed now, since mhdFoam is its own separate header tree no longer reachable via -lOFSolvers's shared include dir):
```
-I./../../modules/openFoamImported/mhdFoam/lnInclude
```

And add the library to `EXE_LIBS`, alongside the others:
```
-L$(FOAM_USER_LIBBIN) -lOFSolvers \
-L$(FOAM_USER_LIBBIN) -lMhdSolver \
```

Then rebuild GeN-Foam itself:
```
cd .../applications/solvers/GeN-Foam
wmake
```




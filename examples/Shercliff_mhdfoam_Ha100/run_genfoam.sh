#!/bin/bash
# Source the OpenFOAM profile

ulimit -s unlimited

cores=32 
solver=GeN-Foam #epotFoam_isoTh  # cEpotmhdEpotMultiRegionFoam
source ~/.openfoam_v2512_profile
source ~/.openfoam_v2512_ffn_profile

#blockMesh   2>&1 | tee log.blockMesh_fluid
blockMesh -region fluid  2>&1 | tee log.blockMesh


if [[ ! $cores -eq 1 ]]
then
   decomposePar -region fluid  2>&1 | tee log.decomposePar
   mpirun -np $cores  $solver  -parallel 2>&1 | tee log.genFoam
#   mpirun -np $cores $solver -parallel  2>&1 | tee log.$solver
   reconstructPar -region fluid 2>&1 | tee log.reconstructPar
else
        $solver 2>&1 | tee log.$solver
fi


# EDIPIC-2D
# this is an example how to install EDIPIC-2D and petsc on perseus

To compile the code in your directory: 

module load openmpi intel
rm -rf *.mod
rm -rf *.o
gcc -c WELL19937a_new.c

mpifort -c -C pic2d_Modules.f90 pic2d_Snapshots.f90 pic2d_Prepare_FFTX_SYSY.f90 \
pic2d_RandomNumberInterface.f90 pic2d_TimeDependences.f90 pic2d_VelocityDistributions.f90 \
pic2d_Setup.f90 pic2d_WallPotentials.f90 pic2d_ParticleExchange.f90 pic2d_BlockSetProc.f90 \
pic2d_IonMoments.f90 pic2d_Checkpoints.f90 pic2d_IonWallCollisions.f90 pic2d_LoadBalancing.f90 \
pic2d_Diagnostics.f90 pic2d_MainProgram.f90 pic2d_ElectricFieldCalc_FFT_X.f90 \
pic2d_Materials.f90 pic2d_ElectronDynamics.f90 pic2d_ElectronExchange.f90 \
pic2d_ElectronMoments.f90 pic2d_ElectronWallCollisions.f90 pic2d_enCollisionsGeneralProc.f90 \
pic2d_enCollisionsSpecificProc.f90 pic2d_ExternalFields.f90 pic2d_HTSetup.f90 pic2d_IonDynamics.f90

##
if you have petsc preinstalled and available as a module (it did not work on perseus),
load the module and finish compilation/linking with the commands below:

module load petsc
mpifort -c pic2d_PETSc_Solver.F90 pic2d_ElectricFieldCalc_PETSc.F90 pic2d_CurProblemValues.f90
mpifort -o ./runpic2dpetsc.out *.o -lpetsc

##
if petsc has to be installed, please use the instructions at https://www.mcs.anl.gov/petsc/documentation/installation.html
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc
module load openmpi intel
./configure --download-hypre
make PETSC_DIR=/scratch/gpfs/your_directory/petsc PETSC_ARCH=arch-linux-c-debug all
make PETSC_DIR=/scratch/gpfs/your_directory/petsc PETSC_ARCH=arch-linux-c-debug check

why hypre is needed?

and provide petsc paths to finish compilation/linking with the commands below:

export PETSC_DIR=/scratch/gpfs/your_directory/petsc
export PETSC_ARCH=arch-linux-c-debug
mpifort -c pic2d_PETSc_Solver.F90 pic2d_ElectricFieldCalc_PETSc.F90 pic2d_CurProblemValues.f90 -C -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include/petsc/finclude
mpifort -o ./runpic2dpetsc.out -C *.o -lpetsc -L${PETSC_DIR}/${PETSC_ARCH}/lib

##

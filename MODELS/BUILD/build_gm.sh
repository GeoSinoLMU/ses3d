ftn -c ../SOURCE/modules_gm.f90
ftn -c ../SOURCE/generate_models.f90
ftn -c ../MODELS_1D/ak135.f90
ftn -c ../MODELS_1D/eumod_bg.f90
ftn -c ../MODELS_1D/homogeneous.f90
ftn -c ../MODELS_1D/homogeneous_plus_Q.f90
ftn -c ../MODELS_1D/prem_iso.f90
ftn -o ../MAIN/gm modules_gm.o generate_models.o ak135.o eumod_bg.o homogeneous.o homogeneous_plus_Q.o prem_iso.o


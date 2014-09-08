ftn -c ../../SOURCE/ses3d_modules.f90
ftn -c ../SOURCE/add_perturbation.f90
ftn -o ../MAIN/add_perturbation ses3d_modules.o add_perturbation.o


# tp-mpi-idris
Ma solution aux TPs MPI de l'IDRIS ( http://www.idris.fr/formations/mpi.html ), avec quelques trucs en bonus (ex: parallélisation hybride MPI + OpenMP du solveur Poisson).

Le code est rédigé en Fortran 95, et a été testé avec l'implémentation MPI mpich2 version 3.0.4 et le compilateur GFortran version 4.8.4.

Il devrait être compatible avec tout compilateur Fortran 95 en mode de respect des standards, et toute implémentation MPI 3.0+.

Pour le support OpenMP (optionnel), l'implémentation doit être compatible avec OpenMP 3.1+.

COMPILER        = 
CFLAGS   	= -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -O3 -march=native -mtune=native -std=c99
LDFLAGS         =


#include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
include ${SLEPC_DIR}/lib/slepc/conf/slepc_common

main: main.o
	-${CLINKER} -o main main.o ${PETSC_SNES_LIB}
	${RM} main.o

main_stab: main_stab.o
	-${CLINKER} -o main_stab main_stab.o ${SLEPC_EPS_LIB} ${PETSC_SNES_LIB}
	${RM} main_stab.o

main_plot: main_plot.c plot/plot_solution.c plot/plot_volume.c
	mpicc $(CFLAGS) -o main_plot main_plot.c -lmgl-mpi -lmgl -lm

#include <stdio.h>
#include <math.h>
#include "util.h"

int file_exist(char fname[], int rank)
{
    int fexist,ierror;
    FILE *fptr;
    if (rank == 0)
    {
        fexist = ( ((fptr=fopen(fname,"rb")) != NULL) ? 1:0 );
        if (fexist == 1) { fclose(fptr); }
    }
    ierror=MPI_Bcast(&fexist,1,MPI_INT,0,MPI_COMM_WORLD);assert(ierror==MPI_SUCCESS);

    return fexist;
}




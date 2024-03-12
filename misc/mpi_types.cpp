#include "mpi_types.h"


#include <stddef.h>

MPI_Datatype mpi_body_type;
MPI_Datatype mpi_cell_type;

void create_mpi_body(){
    int blocklengths[5] = {1, DIM_SIZE, DIM_SIZE, 1, 1};
    MPI_Datatype types[5] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,  MPI_DOUBLE};
    MPI_Aint offsets[5];

    offsets[0] = offsetof(Body, idx);
    offsets[1] = offsetof(Body, pos);
    offsets[2] = offsetof(Body, vel);
    offsets[3] = offsetof(Body, m);
    offsets[4] = offsetof(Body, work);

    MPI_Type_create_struct(5, blocklengths, offsets, types, &mpi_body_type);
    MPI_Type_commit(&mpi_body_type);
}

void create_mpi_cell(){
    int blocklengths[5] = {DIM_SIZE, DIM_SIZE, 1, DIM_SIZE, 1};
    MPI_Datatype types[5] = {MPI_DOUBLE, MPI_DOUBLE,  MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Aint offsets[5];

    offsets[0] = offsetof(MPICell, min_bounds);
    offsets[1] = offsetof(MPICell, max_bounds);
    offsets[2] = offsetof(MPICell, m);
    offsets[3] = offsetof(MPICell, rm); // centre of mass
    offsets[4] = offsetof(MPICell, parent_idx);

    MPI_Type_create_struct(5, blocklengths, offsets, types, &mpi_cell_type);
    MPI_Type_commit(&mpi_cell_type);
}


void init_mpi_types(){
    create_mpi_body();
    create_mpi_cell();
}

void free_mpi_types(){
    MPI_Type_free(&mpi_body_type);
    MPI_Type_free(&mpi_cell_type);
}

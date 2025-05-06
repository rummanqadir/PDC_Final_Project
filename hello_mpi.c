#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int namelen;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Get_processor_name(hostname, &namelen);

    printf("Hello from process %d of %d on %s\n", rank, size, hostname);

    MPI_Finalize();
    return 0;
}
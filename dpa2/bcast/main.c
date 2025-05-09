#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void print_binary_exp(int n) {
    int i;
    for (i = sizeof(n) * 8 - 1; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
    }
    putchar('\n');
}

void print_checksum(int checksum, int rank) {
    const int shift = sizeof(checksum) * 8 - 3;

    checksum = (checksum << shift);
    printf("Checksum of rank %d is: ", rank);
    for (int i = 2; i >= 0; --i) {
        printf("%d", (checksum >> (shift + i) & 1));
    }
    putchar('\n');
}

int main() {
    const int n = 10000000;

    int rank, size;
    int m = 0, total_hits = 0;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int size_of_randvec = 114514;
    int *randvec = (int *)malloc(size_of_randvec * sizeof(int));

    if (rank == 0) {
        int checksum = 0;
        srand(424212);
        for (int i = 0; i < size_of_randvec; ++i) {
            randvec[i] = rand() % 2;
            checksum += randvec[i];
        }

        print_checksum(checksum, rank);

        for (int i = 1; i < size; i++) {
            MPI_Send(randvec, size_of_randvec, MPI_INT, i, 0, MPI_COMM_WORLD);
            printf("Vector sent to rank %d\n", i);
        }
    }

    else /*rank!=0*/ {
        int checksum = 0;
        MPI_Recv(randvec, size_of_randvec, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Vector received by rank %d\n", rank);
        for (int i = 0; i < size_of_randvec; ++i) {
            checksum += randvec[i];
        }
        print_checksum(checksum, rank);
    }

    MPI_Finalize();
}
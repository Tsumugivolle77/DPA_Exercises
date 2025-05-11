#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// Erroneous code

#define PRINT_RECEIVER_END

void print_binary_exp(int n) {
    int i;
    for (i = sizeof(n) * 8 - 1; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
    }
    putchar('\n');
}

void print_checksum(int checksum, int rank) {
#ifndef PRINT_RECEIVER_END
    const int shift = sizeof(checksum) * 8 - 3;

    checksum = (checksum << shift);
    printf("Checksum of rank %d is: ", rank);
    for (int i = 2; i >= 0; --i) {
        printf("%d", (checksum >> (shift + i) & 1));
    }
    putchar('\n');
#endif
}

int count_trailing_zeroes(int i) {
    int count = 0;
    while ((i & 1) == 0 && i != 0) {
        count++;
        i >>= 1;
    }
    return count;
}

void send_msg(int *randvec, int size_of_randvec, int rank, int tree_size) {
    long long trailing_zeroes = count_trailing_zeroes(rank);
    printf("Rank %d has %lld trailing zeroes\n", rank, trailing_zeroes);
    long long k = (tree_size < trailing_zeroes ? tree_size : trailing_zeroes) - 1;

    for (int i = k; i >= 0; --i) {
        printf("Vector sent to rank %d\n", i);
        MPI_Send(randvec, size_of_randvec, MPI_INT, (int)(i + pow(2, k)), 0, MPI_COMM_WORLD);
    }
}

int main() {
    const int n = 10000000;

    int rank, size;
    int m = 0, total_hits = 0;

    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int size_of_randvec = 10;
    int *randvec = (int *)malloc(size_of_randvec * sizeof(int));
    long long tree_size = ceil(log(size));
    printf("Tree size: %lld\n", tree_size);

    if (rank == 0) {
#ifndef PRINT_RECEIVER_END
        int checksum = 0;
        srand(42);
        for (int i = 0; i < size_of_randvec; ++i) {
            randvec[i] = rand() % 2;
            checksum += randvec[i];
        }

        print_checksum(checksum, rank);
#endif
        send_msg(randvec, size_of_randvec, rank, tree_size);
    }

    else /*rank!=0*/ {
        MPI_Recv(randvec, size_of_randvec, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Vector received by rank %d\n", rank);
#ifndef PRINT_RECEIVER_END
        int checksum = 0;
        for (int i = 0; i < size_of_randvec; ++i) {
            checksum += randvec[i];
        }
        print_checksum(checksum, rank);
#endif
        send_msg(randvec, size_of_randvec, rank, tree_size);
    }

    MPI_Finalize();
}
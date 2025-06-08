#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "mmio.h"

// After the build is done, run the program with:
// <executable-path>/dpa4-nr2 <matrix-market-file-path> > output.txt
// s.t. you can read the output clearly.

// define the csr-format graph struct to speed up calculation,
// as a sparse graph was given.
// See: https://netlib.org/linalg/html_templates/node91.html
typedef struct {
    int n;          // #vertices
    int m;          // #edges
    int *degree;
    int *row_ptr;
    int *col_ind;
} csr_graph;

// factory function of csr graph
csr_graph* create_graph(int n, int nz, int *I, int *J) {
    csr_graph *g = (csr_graph*)malloc(sizeof(csr_graph));
    g->n = n;
    g->m = nz;
    g->degree = (int*)calloc(n, sizeof(int));
    
    // store the size of the adjacency lists
    int *sizes = (int*)calloc(n, sizeof(int));
    
    // compute the degree of each vertex
    for (int i = 0; i < nz; i++) {
        int u = I[i], v = J[i];
        // edge like (u, u) will not be counted
        if (u != v) {
            sizes[u]++;
            sizes[v]++;
        }
    }
    
    g->row_ptr = (int*)malloc((n+1) * sizeof(int));
    g->row_ptr[0] = 0;

    // **prefix sum** of sizes[]. can be made parallel if necessary.
    for (int i = 0; i < n; i++) {
        g->row_ptr[i+1] = g->row_ptr[i] + sizes[i];
    }
    
    g->col_ind = (int*)malloc(g->row_ptr[n] * sizeof(int));
    
    int *counters = (int*)calloc(n, sizeof(int));

    // traverse the edges (I[i], J[i])
    for (int i = 0; i < nz; i++) {
        int u = I[i], v = J[i];
        if (u != v) {
            int pos_u = g->row_ptr[u] + counters[u]++;
            int pos_v = g->row_ptr[v] + counters[v]++;
            g->col_ind[pos_u] = v;
            g->col_ind[pos_v] = u;
        }
    }
    
    // prepare the degrees' information for Luby's algorithm (random priority)
    memcpy(g->degree, sizes, n * sizeof(int));
    
    free(sizes);
    free(counters);
    return g;
}

// destructor of csr graph
void free_graph(csr_graph *g) {
    if (g) {
        free(g->degree);
        free(g->row_ptr);
        free(g->col_ind);
        free(g);
    }
}

// By modifying the priority of the nodes, we can implement different algorithms for MIS.
// This condition `active[neighbor] && (priority[neighbor] > priority[i] || (priority[neighbor] == priority[i] && neighbor > i))`
// ensures that the node with the highest priority is selected, and in case of a tie (same degree or random priority), the node with the larger ID is preferred.

// Luby's algorithm for MIS
int* luby_mis_random_priority(csr_graph *g) {
    const int n = g->n;
    int *active = (int*)malloc(n * sizeof(int)); // the active nodes (nodes that are not yet in MIS)
    int *inMIS = (int*)calloc(n, sizeof(int));
    double *priority = (double*)malloc(n * sizeof(double));
    
    // init
    #pragma omp parallel
    {
        unsigned int seed = time(NULL) ^ omp_get_thread_num();
        #pragma omp for
        for (int i = 0; i < n; i++) {
            active[i] = 1;
            priority[i] = (double)rand_r(&seed) / RAND_MAX; // assign random priority
        }
    }

    int n_active = n;
    int iteration = 0;
    
    while (n_active > 0) {
        // calculate the candidate set for MIS
        int *is_candidate = (int*)calloc(n, sizeof(int));
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < n; i++) {
                if (!active[i]) continue;
                
                int is_max = 1;
                int start = g->row_ptr[i];
                int end = g->row_ptr[i+1];
                for (int j = start; j < end; j++) {
                    int neighbor = g->col_ind[j];
                    if (active[neighbor] && (priority[neighbor] > priority[i] || (priority[neighbor] == priority[i] && neighbor > i))) {
                        is_max = 0;
                        break;
                    }
                }
                // nodes with max priority are joined to the candidate set
                if (is_max) is_candidate[i] = 1;
            }
        }

        // join the candidate nodes to the MIS
        // It will not lead to conflict, as the adjacent nodes
        // will not be in the candidate set at the same time.
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            if (is_candidate[i]) {
                inMIS[i] = 1;
                active[i] = 0;
                // remove neighborhood
                int start = g->row_ptr[i];
                int end = g->row_ptr[i+1];
                for (int j = start; j < end; j++) {
                    int neighbor = g->col_ind[j];
                    active[neighbor] = 0;
                }
            }
        }

        n_active = 0;
        #pragma omp parallel for reduction(+:n_active)
        for (int i = 0; i < n; i++) {
            n_active += active[i];
        }

        free(is_candidate);
        iteration++;
        printf("Iteration %d: %d active nodes remaining\n", iteration, n_active);
    }

    free(active);
    free(priority);
    return inMIS;
}

// Bigger ID Bigger Priority for MIS
int* bibp_mis(csr_graph *g) {
    const int n = g->n;
    int *active = (int*)malloc(n * sizeof(int)); // the active nodes (nodes that are not yet in MIS)
    int *inMIS = (int*)calloc(n, sizeof(int));
    double *priority = (double*)malloc(n * sizeof(double));
    
    // init
    #pragma omp parallel
    {
        unsigned int seed = time(NULL) ^ omp_get_thread_num();
        #pragma omp for
        for (int i = 0; i < n; i++) {
            active[i] = 1;
            priority[i] = i; // assign priority as ID
        }
    }

    int n_active = n;
    int iteration = 0;
    
    while (n_active > 0) {
        // calculate the candidate set for MIS
        int *is_candidate = (int*)calloc(n, sizeof(int));
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < n; i++) {
                if (!active[i]) continue;
                
                int is_max = 1;
                int start = g->row_ptr[i];
                int end = g->row_ptr[i+1];
                for (int j = start; j < end; j++) {
                    int neighbor = g->col_ind[j];
                    if (active[neighbor] && priority[neighbor] > priority[i]) {
                        is_max = 0;
                        break;
                    }
                }
                // nodes with max priority are joined to the candidate set
                if (is_max) is_candidate[i] = 1;
            }
        }

        // join the candidate nodes to the MIS
        // It will not lead to conflict, as the adjacent nodes
        // will not be in the candidate set at the same time.
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            if (is_candidate[i]) {
                inMIS[i] = 1;
                active[i] = 0;
                // remove neighborhood
                int start = g->row_ptr[i];
                int end = g->row_ptr[i+1];
                for (int j = start; j < end; j++) {
                    int neighbor = g->col_ind[j];
                    active[neighbor] = 0;
                }
            }
        }

        n_active = 0;
        #pragma omp parallel for reduction(+:n_active)
        for (int i = 0; i < n; i++) {
            n_active += active[i];
        }

        free(is_candidate);
        iteration++;
        printf("Iteration %d: %d active nodes remaining\n", iteration, n_active);
    }

    free(active);
    free(priority);
    return inMIS;
}

// Bigger Degree Bigger Priority for MIS
int* bdbp_mis(csr_graph *g) {
    const int n = g->n;
    int *active = (int*)malloc(n * sizeof(int)); // the active nodes (nodes that are not yet in MIS)
    int *inMIS = (int*)calloc(n, sizeof(int));
    double *priority = (double*)malloc(n * sizeof(double));
    
    // init
    #pragma omp parallel
    {
        unsigned int seed = time(NULL) ^ omp_get_thread_num();
        #pragma omp for
        for (int i = 0; i < n; i++) {
            active[i] = 1;
            priority[i] = g->degree[i]; // assign priority as degree
        }
    }

    int n_active = n;
    int iteration = 0;
    
    while (n_active > 0) {
        // calculate the candidate set for MIS
        int *is_candidate = (int*)calloc(n, sizeof(int));
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < n; i++) {
                if (!active[i]) continue;
                
                int is_max = 1;
                int start = g->row_ptr[i];
                int end = g->row_ptr[i+1];
                for (int j = start; j < end; j++) {
                    int neighbor = g->col_ind[j];
                    if (active[neighbor] && (priority[neighbor] > priority[i] || (priority[neighbor] == priority[i] && neighbor > i))) {
                        is_max = 0;
                        break;
                    }
                }
                // nodes with max priority are joined to the candidate set
                if (is_max) is_candidate[i] = 1;
            }
        }

        // join the candidate nodes to the MIS
        // It will not lead to conflict, as the adjacent nodes
        // will not be in the candidate set at the same time.
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            if (is_candidate[i]) {
                inMIS[i] = 1;
                active[i] = 0;
                // remove neighborhood
                int start = g->row_ptr[i];
                int end = g->row_ptr[i+1];
                for (int j = start; j < end; j++) {
                    int neighbor = g->col_ind[j];
                    active[neighbor] = 0;
                }
            }
        }

        n_active = 0;
        #pragma omp parallel for reduction(+:n_active)
        for (int i = 0; i < n; i++) {
            n_active += active[i];
        }

        free(is_candidate);
        iteration++;
        printf("Iteration %d: %d active nodes remaining\n", iteration, n_active);
    }

    free(active);
    free(priority);
    return inMIS;
}

int main(int argc, char *argv[]) {
    // Read the matrix market file. Copied from the example code
    if (argc < 2) {
        fprintf(stderr, "Usage: %s [matrix-market-filename]\n", argv[0]);
        exit(1);
    }

    char *filename = argv[1];
    
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("Error opening file");
        exit(1);
    }
    
    MM_typecode matcode;
    if (mm_read_banner(f, &matcode) != 0) {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }
    
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
        printf("Complex matrices are not supported.\n");
        exit(1);
    }
    
    int M, N, nz;
    if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0) {
        printf("Error reading matrix size.\n");
        exit(1);
    }
    
    int *I = (int *)malloc(nz * sizeof(int));
    int *J = (int *)malloc(nz * sizeof(int));
    
    for (int i = 0; i < nz; i++) {
        int row, col;
        double val; // value was not even given in the matrix market file.
        if (mm_is_pattern(matcode)) {
            if (fscanf(f, "%d %d\n", &row, &col) != 2) {
                printf("Error reading data at line %d\n", i+1);
                exit(1);
            }
        } else {
            if (fscanf(f, "%d %d %lf\n", &row, &col, &val) != 3) {
                printf("Error reading data at line %d\n", i+1);
                exit(1);
            }
        }
        // 1-based to 0-based (indices start from 1/0)
        I[i] = row - 1;
        J[i] = col - 1;
    }
    fclose(f);
    
    csr_graph *g = create_graph(M, nz, I, J);
    printf("csr graph created with %d nodes and %d edges\n", g->n, g->m);
    
    // compute mip and measure the time elapsed
    double start_time = omp_get_wtime();
    int *mis = luby_mis_random_priority(g);
    // int *mis = bdbp_mis(g);
    // int *mis = bibp_mis(g);
    double end_time = omp_get_wtime();
    
    int mis_size = 0;
    for (int i = 0; i < g->n; i++) {
        if (mis[i]) mis_size++;
    }
    
    printf("MIS size: %d (%.2f%% of nodes)\n", mis_size, 100.0 * mis_size / g->n);
    printf("Computation time: %.4f seconds\n", end_time - start_time);
    
    printf("Verifying MIS...\n");
    int valid = 1;
    #pragma omp parallel for
    for (int i = 0; i < g->n; i++) {
        if (mis[i]) {
            int start = g->row_ptr[i];
            int end = g->row_ptr[i+1];
            for (int j = start; j < end; j++) {
                int neighbor = g->col_ind[j];
                if (mis[neighbor]) {
                    #pragma omp critical
                    {
                        printf("Conflict between %d and %d\n", i, neighbor);
                        valid = 0;
                    }
                }
            }
        }
    }
    
    if (valid) {
        printf("MIS is valid: no adjacent nodes in the set\n");
    } else {
        printf("MIS is NOT valid!\n");
    }
    
    free(I);
    free(J);
    free(mis);
    free_graph(g);
    
    return 0;
}
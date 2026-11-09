#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Matrixs.h"
#include "Richardson.h"

#ifndef PROJECT_ROOT_DIR
#define PROJECT_ROOT_DIR "."
#endif

#define MATRIX_SIZE 10
#define MAX_CYCLES 100000
#define FIXED_EPS 1e-10
#define DECREASE_FACTOR 10.0
#define RANDOM_SEED 1

static Matrix* createZeroVector(int n)
{
    Matrix* v = newMatrix(n, 1);
    if (v == NULL) return NULL;
    for (int i = 0; i < n; ++i) v->d[i] = 0.0;
    return v;
}

static int buildSystem(const Matrix* Q, const Matrix* D, const Matrix* x_star,
                       Matrix** A_out, Matrix** b_out)
{
    Matrix* Qt = transposeMatrix((Matrix*)Q);
    Matrix* QD = multiplyMatrix((Matrix*)Q, (Matrix*)D);
    Matrix* A = (QD != NULL && Qt != NULL) ? multiplyMatrix(QD, Qt) : NULL;
    Matrix* b = (A != NULL) ? multiplyMatrix(A, (Matrix*)x_star) : NULL;
    freeMatrix(Qt);
    freeMatrix(QD);
    if (A == NULL || b == NULL) {
        freeMatrix(A);
        freeMatrix(b);
        return 0;
    }
    *A_out = A;
    *b_out = b;
    return 1;
}

static int openCsv(const char* relativePath, char* path, size_t pathSize, FILE** output)
{
    snprintf(path, pathSize, "%s/%s", PROJECT_ROOT_DIR, relativePath);
    // printf("%s\n", path);
    *output = fopen(path, "w");
    if (*output == NULL) {
        perror(relativePath);
        return 0;
    }
    return 1;
}

static int runOne(const Matrix* A, const Matrix* b, const Matrix* x0,
                  const Matrix* x_star, const Matrix* D, double eps,
                  const int order[3], const char* experiment,
                  RichardsonResult* result, double* det_out,
                  double* delta_out, double* Delta_out, double alpha_out[3], double* M3_out)
{
    double delta, Delta;
    if (!getSpectralBoundsD(D, &delta, &Delta)) return 0;
    double alphas[3];
    RichardsonParams(delta, Delta, alphas);
    double M3 = RichardsonM3(delta, Delta);
    if (!isfinite(M3) || M3 >= 1.0) return 0;

    *result = RichardsonMethod(A, b, x0, x_star, alphas, order, M3, eps, MAX_CYCLES);
    if (result->x == NULL) return 0;

    *det_out = calculateDetD(D);
    *delta_out = delta;
    *Delta_out = Delta;
    *M3_out = M3;
    for (int i = 0; i < 3; ++i) alpha_out[i] = alphas[i];
    (void)experiment;
    return 1;
}

static void writeCommonHeader(FILE* f)
{
    fprintf(f, "n,seed,max_cycles,eps,det_D,delta,Delta,alpha0,alpha1,alpha2,M3,order,cycles,iterations,stopping_condition_met,actual_accuracy_met,error_estimate,actual_error,residual_norm,first_cycle_error_step1,first_cycle_error_step2,first_cycle_error_step3\n");
}

static void writeCommonRow(FILE* f, double eps, double det, double delta, double Delta,
                           const double alpha[3], double M3, const int order[3],
                           const RichardsonResult* r)
{
    int actual_accuracy_met = (isfinite(r->actual_error) && r->actual_error <= eps);
    fprintf(f, "%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,\"%d-%d-%d\",%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            MATRIX_SIZE, RANDOM_SEED, MAX_CYCLES, eps, det, delta, Delta,
            alpha[0], alpha[1], alpha[2], M3, order[0], order[1], order[2],
            r->cycles, r->iterations, r->converged, actual_accuracy_met,
            r->error_estimate, r->actual_error, r->residual_norm,
            r->first_cycle_error_step1, r->first_cycle_error_step2, r->first_cycle_error_step3);
}

static int experimentEps(const Matrix* Q, const Matrix* D, const Matrix* x_star, const Matrix* x0)
{
    FILE* f;
    char path[1024];
    if (!openCsv("data/experiment_eps.csv", path, sizeof(path), &f)) return 0;
    writeCommonHeader(f);
    const int order[3] = {0, 1, 2};
    for (int exponent = 1; exponent <= 16; ++exponent) {
        double eps = pow(10.0, -exponent);
        Matrix *A = NULL, *b = NULL;
        if (!buildSystem(Q, D, x_star, &A, &b)) { fclose(f); return 0; }
        RichardsonResult r;
        double det, delta, Delta, alpha[3], M3;
        int ok = runOne(A, b, x0, x_star, D, eps, order, "eps", &r, &det, &delta, &Delta, alpha, &M3);
        if (ok) {
            writeCommonRow(f, eps, det, delta, Delta, alpha, M3, order, &r);
            freeMatrix(r.x);
        } else {
            fprintf(stderr, "Experiment eps failed for eps=%.3e\n", eps);
        }
        freeMatrix(A); freeMatrix(b);
    }
    fclose(f);
    return 1;
}

static int experimentDet(const Matrix* Q, const Matrix* D_initial, const Matrix* x_star, const Matrix* x0)
{
    FILE* f;
    char path[1024];
    if (!openCsv("data/experiment_det.csv", path, sizeof(path), &f)) return 0;
    writeCommonHeader(f);
    const int order[3] = {0, 1, 2};
    Matrix* D = copyMatrix(D_initial);
    if (D == NULL) { fclose(f); return 0; }

    for (int exponent = 0; exponent <= 16; ++exponent) {
        Matrix *A = NULL, *b = NULL;
        if (!buildSystem(Q, D, x_star, &A, &b)) {
            freeMatrix(D); fclose(f); return 0;
        }
        RichardsonResult r;
        double det, delta, Delta, alpha[3], M3;
        int ok = runOne(A, b, x0, x_star, D, FIXED_EPS, order, "det", &r,
                        &det, &delta, &Delta, alpha, &M3);
        if (ok) {
            writeCommonRow(f, FIXED_EPS, det, delta, Delta, alpha, M3, order, &r);
            freeMatrix(r.x);
        } else {
            fprintf(stderr, "Experiment det failed at target exponent %d\n", exponent);
        }
        freeMatrix(A); freeMatrix(b);

        if (exponent < 16) {
            int changed = -1;
            Matrix* nextD = decreaseDElement(D, DECREASE_FACTOR, &changed);
            if (nextD == NULL) { freeMatrix(D); fclose(f); return 0; }
            freeMatrix(D);
            D = nextD;
        }
    }
    freeMatrix(D);
    fclose(f);
    return 1;
}

static int experimentOrder(const Matrix* Q, const Matrix* D, const Matrix* x_star, const Matrix* x0)
{
    FILE* f;
    char path[1024];
    if (!openCsv("data/experiment_order.csv", path, sizeof(path), &f)) return 0;
    writeCommonHeader(f);
    const int permutations[6][3] = {
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2},
        {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
    };
    Matrix *A = NULL, *b = NULL;
    if (!buildSystem(Q, D, x_star, &A, &b)) { fclose(f); return 0; }

    for (int i = 0; i < 6; ++i) {
        RichardsonResult r;
        double det, delta, Delta, alpha[3], M3;
        int ok = runOne(A, b, x0, x_star, D, FIXED_EPS, permutations[i], "order", &r,
                        &det, &delta, &Delta, alpha, &M3);
        if (ok) {
            writeCommonRow(f, FIXED_EPS, det, delta, Delta, alpha, M3, permutations[i], &r);
            freeMatrix(r.x);
        } else {
            fprintf(stderr, "Experiment order failed for permutation %d\n", i);
        }
    }
    freeMatrix(A); freeMatrix(b);
    fclose(f);
    return 1;
}

int main(void)
{
    srand(RANDOM_SEED); /* Воспроизводимость всех случайных исходных данных */
    Matrix* Q = createQMatrix(MATRIX_SIZE);
    Matrix* D = builtDMatrix(MATRIX_SIZE, 1.0);
    Matrix* x_star = createRandomMatrix(MATRIX_SIZE, 1, -5.0, 5.0);
    // Matrix* D1 = builtDMatrix(MATRIX_SIZE, 1e-3);
    Matrix* x0 = createZeroVector(MATRIX_SIZE);

    if (Q == NULL || D == NULL || x_star == NULL || x0 == NULL) {
        fprintf(stderr, "Failed to create initial data.\n");
        freeMatrix(Q); freeMatrix(D); freeMatrix(x_star); freeMatrix(x0);
        return EXIT_FAILURE;
    }

    printf("Initial det(D) = %.17g\n", calculateDetD(D));
    // printf("Initial det(D1) = %.17g\n", calculateDetD(D1));
    printf("Running experiment 1: epsilon...\n");
    int ok1 = experimentEps(Q, D, x_star, x0);
    printf("Running experiment 2: determinant...\n");
    int ok2 = experimentDet(Q, D, x_star, x0);
    printf("Running experiment 3: parameter order...\n");
    int ok3 = experimentOrder(Q, D, x_star, x0);

    freeMatrix(Q); freeMatrix(D); freeMatrix(x_star); freeMatrix(x0);

    if (!ok1 || !ok2 || !ok3) {
        fprintf(stderr, "One or more experiments could not be completed.\n");
        return EXIT_FAILURE;
    }

    printf("Experiments completed. CSV files: experiment_eps.csv, experiment_det.csv, experiment_order.csv\n");
    return EXIT_SUCCESS;
}

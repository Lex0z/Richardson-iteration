//
// Created by Алексей on 11.10.2026.
//

#include "Matrixs.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DEBUG 0

Matrix* copyMatrix(const Matrix* source) {
    if (source == NULL || source->n <= 0 || source->m <= 0 || source->d == NULL) {
        return NULL;
    }
    Matrix* copy = newMatrix(source->n, source->m);
    if (copy == NULL || copy->d == NULL) {
        if (copy != NULL) freeMatrix(copy);
        return NULL;
    }
    for (int i = 0; i < source->n * source->m; ++i) {
        copy->d[i] = source->d[i];
    }
    return copy;
}

Matrix* newMatrix(int n, int m) {
    if (n <= 0 || m <= 0) return NULL;
    Matrix *matrix = (Matrix*)malloc(sizeof(Matrix));
    if (matrix == NULL) return NULL;
    matrix->n = n;
    matrix->m = m;
    matrix->d = (double *)malloc(sizeof(double) * (size_t)m * (size_t)n);
    if (matrix->d == NULL) {
        free(matrix);
        return NULL;
    }
    return matrix;
}

Matrix* createMatrix(int n, int m, double* d) {
    if (d == NULL) return NULL;
    Matrix *matrix = newMatrix(n, m);
    if (matrix == NULL) return NULL;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (DEBUG) {printf("$%lf ", d[i * m + j]);}
            matrix->d[i * m + j] = d[i * m + j];
        }
    }
    return matrix;
}

void printMatrix(Matrix* matrix) {
    printf("\nMatrix %d x %d:\n", matrix->n, matrix->m);
    for (int i = 0; i < matrix->n; i++) {
        for (int j = 0; j < matrix->m; j++) {
            printf("% e ",matrix->d[i * matrix->m  + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void freeMatrix(Matrix* matrix) {
    if (matrix == NULL) return;
    free(matrix->d);
    free(matrix);
}

double normVector2(const Matrix* v) {
    if (v == NULL || v->m != 1 || v->d == NULL) return NAN;
    double sum = 0.0;
    for (int i = 0; i < v->n; ++i) sum += v->d[i] * v->d[i];
    return sqrt(sum);
}

double normMatrix(Matrix* T) {
    if (T == NULL || T->d == NULL) return NAN;
    double mx = 0.0;
    for (int j = 0; j < T->m; j++) {
        double sum = 0.0;
        for (int i = 0; i < T->n; i++) {
            sum += fabs(T->d[i*T->m + j]);
        }
        if (sum > mx) { mx = sum; }
    }
    return mx;
}

Matrix* multiplyMatrix(Matrix* A, Matrix* B) {
    if (A == NULL || B == NULL || A->d == NULL || B->d == NULL || A->m != B->n) return NULL;
    int n = A->n;
    int p = A->m;
    int m = B->m;
    Matrix *C = newMatrix(n, m);
    if (C == NULL) return NULL;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            for (int k = 0; k < p; k++) {
                sum += A->d[i * A->m + k] * B->d[k * B->m + j];
            }
            C->d[i * C->m + j] = sum;
        }
    }
    return C;
}

Matrix* createRandomMatrix(int n, int m, double a, double b) {
    if (n <= 0 || m <= 0 || !isfinite(a) || !isfinite(b) || b < a) return NULL;
    double *d = (double *)malloc(sizeof(double) * (size_t)n * (size_t)m);
    if (d == NULL) return NULL;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            d[i*m + j] = a + (b-a) * (double) rand() / (double)RAND_MAX;
        }
    }
    Matrix* T = createMatrix(n, m, d);
    free(d);
    return T;
}

Matrix* transposeMatrix(Matrix* M) {
    if (M == NULL || M->d == NULL) return NULL;
    double *d = (double *)malloc(sizeof(double) * (size_t)M->n * (size_t)M->m);
    if (d == NULL) return NULL;
    for (int i = 0; i < M->n; i++) {
        for (int j = 0; j < M->m; j++) {
            d[j*M->n + i] = M->d[i*M->m + j];
        }
    }
    Matrix* T = createMatrix(M->m, M->n, d);
    free(d);
    return T;
}

Matrix* addMatrix(Matrix* A, Matrix* B, char op) {
    if (A == NULL || B == NULL || A->d == NULL || B->d == NULL || A->m != B->m || A->n != B->n) {
        printf("Matrix dimensions do not match");
        return NULL;
    }
    double k = 0;
    switch (op) {
        case '+':
            k = 1;
            break;
        case '-':
            k = -1;
            break;
        default:
            printf("Invalid operation");
            return NULL;
            break;
    }

    Matrix* C = newMatrix(A->n, A->m);
    if (C == NULL) return NULL;

    for (int i = 0; i < A->n; i++) {
        for (int j = 0; j < A->m; j++) {
            C->d[i*C->m + j] = A->d[i*A->m + j] + k * B->d[i*B->m + j];
        }
    }
    return C;
}

double* maxInMatrix(Matrix* M) {
    if (M == NULL || M->d == NULL || M->n <= 0 || M->m <= 0) return NULL;
    double* max = &M->d[0];
    for (int i = 0; i < M->n * M->m; ++i) {
        if (*max < M->d[i]) max = &M->d[i];
    }
    return max;
}

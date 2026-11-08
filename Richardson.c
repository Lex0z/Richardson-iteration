//
// Created by Алексей on 08.11.2026.
//

#include "Richardson.h"
#include "Matrixs.h"

#include <math.h>
#include <stdio.h>

int CheckOrder(const int order[3])
{
    int used[3] = {0, 0, 0};

    for (int i = 0; i < 3; i++)
    {
        if (order[i] < 0 || order[i] >= 3)
            return 0;

        if (used[order[i]])
            return 0;

        used[order[i]] = 1;
    }

    return 1;
}

Matrix* RichardsonStep(const Matrix* A, const Matrix* b, const Matrix* x, double alpha) {
    if (A == NULL || b == NULL || x == NULL)
        return NULL;

    if (A->n != A->m ||
        b->m != 1 ||
        x->m != 1 ||
        A->n != b->n ||
        A->n != x->n)
        return NULL;

    int n = A->n;

    Matrix* x_new = newMatrix(n, n);
    if (x_new == NULL){ return  NULL; }

    Matrix* Ax = multiplyMatrix(A, x);
    if (Ax == NULL) {
        freeMatrix(x_new);
        return NULL;
    }
    for (int i = 0; i < n; i++) {
        x_new->d[i] = x->d[i] - alpha * (Ax->d[i] - b->d[i]);
    }
    return x_new;
}

Matrix* RichardsonMethod(const Matrix* A, const Matrix* b, const Matrix* x_fix,
    const double* alphas, const int* order, double M, double eps, int MAX_ITER, int* iterCount, double* error) {
    if (A == NULL || b == NULL || x_fix == NULL){ return NULL; }
    if (CheckOrder(order) != 1) { return NULL; }

    int n = A->n;
    Matrix* x = newMatrix(n, 1);
    if (x == NULL){ return  NULL; }
    for (int i = 0; i < n; i++) {
        x->d[i] = 0.0;
    }
    double C = eps*(1 - M)/M;

    int iter = 0;
    double err = 0.0;

    while (iter < MAX_ITER) {
        double alpha = alphas[iter%3];

        Matrix* x_new = RichardsonStep(A, b, x, alpha);

        Matrix* e = (x_new != NULL) ? addMatrix(x_new, x, '-') : NULL;
        err = normVector2(e);
        iter++;
        freeMatrix(e);
        freeMatrix(x);
        x = x_new;
        if (err < C) { break; }
    }

    *iterCount = iter;
    *error = err;
    return x;
}

Matrix* createQMatrix(int n) {
    if (n <= 0) return NULL;
    Matrix* w = createRandomMatrix(n, 1, 0 , 1);
    if (w == NULL) return NULL;

    double norm2 = normVector2(w);
    if (norm2 == 0.0) {
        freeMatrix(w);
        return NULL;
    }
    for (int i = 0; i < n; ++i) w->d[i] /= norm2;

    Matrix* Q = newMatrix(n, n);
    if (Q == NULL || Q->d == NULL) {
        if (Q != NULL) freeMatrix(Q);
        freeMatrix(w);
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const double identity = (i == j) ? 1.0 : 0.0;
            Q->d[i * n + j] = identity - 2.0 * w->d[i] * w->d[j];
        }
    }
    freeMatrix(w);
    return Q;
}

Matrix* builtDMatrix(int n, double det) {
    if (n <= 0) return NULL;
    Matrix* T = createRandomMatrix(n, 1, 0,10);
    Matrix* D = newMatrix(n, n);
    if (D == NULL || D->d == NULL || T == NULL) {
        if (D != NULL) freeMatrix(D);
        if (T != NULL) freeMatrix(T);
        return NULL;
    }
    double prod = 1;
    for (int i = 0; i < n; i++) {
        prod*= T->d[i];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) D->d[i * n + j] = 0.0;
        D->d[i * n + i] = T->d[i]/prod;
    }
    return D;
}

Matrix* decreaseDElement(Matrix* D, double value) {
    if (D == NULL || value == 0) return NULL;
    int n = D->n;
    Matrix* T = copyMatrix(D);
    double* max = (T!= NULL)? maxInMatrix(T) : NULL;
    if (max != NULL) {
        *max /= value;
    }
    return T;
}
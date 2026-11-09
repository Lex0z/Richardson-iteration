#include "Richardson.h"

#include <math.h>
#include <stdlib.h>

int CheckOrder(const int order[RICHARDSON_M])
{
    if (order == NULL) return 0;
    int used[RICHARDSON_M] = {0, 0, 0};
    for (int i = 0; i < RICHARDSON_M; ++i) {
        if (order[i] < 0 || order[i] >= RICHARDSON_M || used[order[i]]) return 0;
        used[order[i]] = 1;
    }
    return 1;
}

void RichardsonParams(double delta, double Delta, double alphas[RICHARDSON_M])
{
    if (alphas == NULL || !(delta > 0.0) || !(Delta >= delta)) return;

    const double c = sqrt(3.0) / 2.0;
    const double cosines[RICHARDSON_M] = {c, 0.0, -c};
    for (int i = 0; i < RICHARDSON_M; ++i) {
        double denominator = delta + Delta + (Delta - delta) * cosines[i];
        alphas[i] = 2.0 / denominator;
    }
}

double RichardsonM3(double delta, double Delta)
{
    if (!(delta > 0.0) || !(Delta > delta)) return NAN;
    double z = (Delta + delta) / (Delta - delta);
    double t3 = 4.0 * z * z * z - 3.0 * z;
    if (!isfinite(t3) || t3 == 0.0) return NAN;
    return 1.0 / fabs(t3);
}

Matrix* RichardsonStep(const Matrix* A, const Matrix* b, const Matrix* x, double alpha)
{
    if (A == NULL || b == NULL || x == NULL || A->d == NULL || b->d == NULL || x->d == NULL)
        return NULL;
    if (A->n != A->m || b->m != 1 || x->m != 1 || A->n != b->n || A->n != x->n)
        return NULL;

    const int n = A->n;
    Matrix* x_new = newMatrix(n, 1);
    if (x_new == NULL) return NULL;

    for (int i = 0; i < n; ++i) {
        double Ax_i = 0.0;
        for (int j = 0; j < n; ++j)
            Ax_i += A->d[i * A->m + j] * x->d[j];
        x_new->d[i] = x->d[i] + alpha * (b->d[i] - Ax_i);
    }
    return x_new;
}

static double vectorDifferenceNorm(const Matrix* a, const Matrix* b)
{
    if (a == NULL || b == NULL || a->m != 1 || b->m != 1 || a->n != b->n)
        return NAN;
    double sum = 0.0;
    for (int i = 0; i < a->n; ++i) {
        double d = a->d[i] - b->d[i];
        sum += d * d;
    }
    return sqrt(sum);
}

static double residualNorm(const Matrix* A, const Matrix* b, const Matrix* x)
{
    double sum = 0.0;
    for (int i = 0; i < A->n; ++i) {
        double ax = 0.0;
        for (int j = 0; j < A->n; ++j)
            ax += A->d[i * A->m + j] * x->d[j];
        double r = b->d[i] - ax;
        sum += r * r;
    }
    return sqrt(sum);
}

RichardsonResult RichardsonMethod(
    const Matrix* A, const Matrix* b, const Matrix* x0, const Matrix* x_star,
    const double alphas[RICHARDSON_M], const int order[RICHARDSON_M],
    double M3, double eps, int max_cycles)
{
    RichardsonResult result = {0};
    result.error_estimate = NAN;
    result.residual_norm = NAN;
    result.actual_error = NAN;
    result.first_cycle_error_step1 = NAN;
    result.first_cycle_error_step2 = NAN;
    result.first_cycle_error_step3 = NAN;

    if (A == NULL || b == NULL || x0 == NULL || x_star == NULL || alphas == NULL ||
        A->d == NULL || b->d == NULL || x0->d == NULL || x_star->d == NULL ||
        A->n != A->m || b->m != 1 || x0->m != 1 || x_star->m != 1 ||
        A->n != b->n || A->n != x0->n || A->n != x_star->n ||
        !CheckOrder(order) || !(M3 >= 0.0 && M3 < 1.0) || !(eps > 0.0) || max_cycles <= 0)
        return result;

    for (int i = 0; i < RICHARDSON_M; ++i)
        if (!isfinite(alphas[i])) return result;

    Matrix* x = copyMatrix(x0);
    Matrix* previous_cycle = NULL;
    if (x == NULL) return result;

    for (int p = 1; p <= max_cycles; ++p) {
        for (int s = 0; s < RICHARDSON_M; ++s) {
            Matrix* next = RichardsonStep(A, b, x, alphas[order[s]]);
            if (next == NULL) {
                freeMatrix(x);
                freeMatrix(previous_cycle);
                return result;
            }
            freeMatrix(x);
            x = next;

            if (p == 1) {
                double step_error = vectorDifferenceNorm(x, x_star);
                if (s == 0) result.first_cycle_error_step1 = step_error;
                if (s == 1) result.first_cycle_error_step2 = step_error;
                if (s == 2) result.first_cycle_error_step3 = step_error;
            }
        }

        result.cycles = p;
        result.iterations = p * RICHARDSON_M;

        if (previous_cycle != NULL) {
            double cycle_difference = vectorDifferenceNorm(x, previous_cycle);
            result.error_estimate = (M3 == 0.0) ? 0.0 : (M3 / (1.0 - M3)) * cycle_difference;
            if (result.error_estimate <= eps) {
                result.converged = 1;
                freeMatrix(previous_cycle);
                previous_cycle = NULL;
                break;
            }
        }

        freeMatrix(previous_cycle);
        previous_cycle = copyMatrix(x);
        if (previous_cycle == NULL) {
            freeMatrix(x);
            result.x = NULL;
            return result;
        }
    }

    freeMatrix(previous_cycle);
    result.x = x;
    result.residual_norm = residualNorm(A, b, x);
    result.actual_error = vectorDifferenceNorm(x, x_star);
    return result;
}

Matrix* createQMatrix(int n)
{
    if (n <= 0) return NULL;
    Matrix* w = createRandomMatrix(n, 1, 0.0, 1.0);
    if (w == NULL) return NULL;
    double norm = normVector2(w);
    if (!(norm > 0.0) || !isfinite(norm)) {
        freeMatrix(w);
        return NULL;
    }
    for (int i = 0; i < n; ++i) w->d[i] /= norm;

    Matrix* Q = newMatrix(n, n);
    if (Q == NULL) {
        freeMatrix(w);
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double identity = (i == j) ? 1.0 : 0.0;
            Q->d[i * n + j] = identity - 2.0 * w->d[i] * w->d[j];
        }
    }
    freeMatrix(w);
    return Q;
}

Matrix* builtDMatrix(int n, double det)
{
    if (n <= 0 || !(det > 0.0) || !isfinite(det)) return NULL;

    Matrix* D = newMatrix(n, n);
    if (D == NULL) return NULL;
    for (int i = 0; i < n * n; ++i) D->d[i] = 0.0;

    /* Фиксированный seed задаётся один раз в main; здесь генерируются положительные числа. */
    double product = 1.0;
    for (int i = 0; i < n - 1; ++i) {
        double value = 0.5 + (double)rand() / (double)RAND_MAX * 1.5;
        D->d[i * n + i] = value;
        product *= value;
    }
    D->d[(n - 1) * n + (n - 1)] = det / product;
    return D;
}

Matrix* decreaseDElement(const Matrix* D, double factor, int* changed_index)
{
    if (D == NULL || D->d == NULL || D->n != D->m || !(factor > 1.0) || !isfinite(factor))
        return NULL;

    int max_index = 0;
    double max_value = D->d[0];
    for (int i = 1; i < D->n; ++i) {
        double value = D->d[i * D->m + i];
        if (value > max_value) {
            max_value = value;
            max_index = i;
        }
    }
    if (!(max_value > 0.0)) return NULL;

    Matrix* result = copyMatrix(D);
    if (result == NULL) return NULL;
    result->d[max_index * result->m + max_index] /= factor;
    if (changed_index != NULL) *changed_index = max_index;
    return result;
}

double calculateDetD(const Matrix* D)
{
    if (D == NULL || D->d == NULL || D->n != D->m) return NAN;
    double det = 1.0;
    for (int i = 0; i < D->n; ++i) {
        double value = D->d[i * D->m + i];
        if (!(value > 0.0)) return NAN;
        det *= value;
    }
    return det;
}

int getSpectralBoundsD(const Matrix* D, double* delta, double* Delta)
{
    if (D == NULL || D->d == NULL || D->n != D->m || delta == NULL || Delta == NULL)
        return 0;
    double mn = D->d[0];
    double mx = D->d[0];
    for (int i = 0; i < D->n; ++i) {
        double value = D->d[i * D->m + i];
        if (!(value > 0.0) || !isfinite(value)) return 0;
        if (value < mn) mn = value;
        if (value > mx) mx = value;
    }
    *delta = mn;
    *Delta = mx;
    return 1;
}

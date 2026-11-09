#ifndef RICHARDSON_ITERATION_RICHARDSON_H
#define RICHARDSON_ITERATION_RICHARDSON_H

#include "Matrixs.h"

#define RICHARDSON_M 3

typedef struct {
    Matrix* x;                 /* Итоговое приближение; память принадлежит вызывающему коду */
    int cycles;                /* Число завершённых полных циклов */
    int iterations;            /* Число отдельных итераций (3 * cycles) */
    int converged;             /* 1 — точность достигнута, 0 — исчерпан лимит */
    double error_estimate;     /* Апостериорная оценка ошибки */
    double residual_norm;      /* ||b - A*x||_2 */
    double actual_error;       /* ||x - x_star||_2; доступна в тестах с известным решением */
    double first_cycle_error_step1;
    double first_cycle_error_step2;
    double first_cycle_error_step3;
} RichardsonResult;

int CheckOrder(const int order[RICHARDSON_M]);
void RichardsonParams(double delta, double Delta, double alphas[RICHARDSON_M]);
double RichardsonM3(double delta, double Delta);

Matrix* RichardsonStep(const Matrix* A, const Matrix* b, const Matrix* x, double alpha);
RichardsonResult RichardsonMethod(
    const Matrix* A,
    const Matrix* b,
    const Matrix* x0,
    const Matrix* x_star,
    const double alphas[RICHARDSON_M],
    const int order[RICHARDSON_M],
    double M3,
    double eps,
    int max_cycles
);

Matrix* createQMatrix(int n);
Matrix* builtDMatrix(int n, double det);
Matrix* decreaseDElement(const Matrix* D, double factor, int* changed_index);
double calculateDetD(const Matrix* D);
int getSpectralBoundsD(const Matrix* D, double* delta, double* Delta);

#endif

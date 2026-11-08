//
// Created by Алексей on 08.11.2026.
//

#ifndef RICHARDSON_ITERATION_RICHARDSON_H
#define RICHARDSON_ITERATION_RICHARDSON_H

#include "Matrixs.h"

int CheckOrder(const int order[3]);

Matrix* RichardsonStep(const Matrix* A, const Matrix* b, const Matrix* x, double alpha);

Matrix* createQMatrix(int n);
Matrix* builtDMatrix(int n, double det);
Matrix* decreaseDElement(Matrix* D, double value);

#endif //RICHARDSON_ITERATION_RICHARDSON_H
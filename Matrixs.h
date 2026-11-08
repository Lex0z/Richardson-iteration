#ifndef LAB2_MATRIXS_H
#define LAB2_MATRIXS_H

typedef struct {
    int n; // rows
    int m; // columns
    double* d;
} Matrix;

Matrix* newMatrix(int n, int m);
Matrix* createMatrix(int n, int m, double* d);
void printMatrix(Matrix* matrix);
void freeMatrix(Matrix* matrix);
double normMatrix(Matrix* T); // matrix 1-norm

double normVector2(const Matrix* v); // Euclidean norm of a column vector
Matrix* multiplyMatrix(Matrix* A, Matrix* B);
Matrix* createRandomMatrix(int n, int m);
Matrix* transposeMatrix(Matrix* M);
Matrix* addMatrix(Matrix* A, Matrix* B, char operation);

#endif

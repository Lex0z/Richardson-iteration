#include <stdio.h>
#include <stdlib.h>

#include "Matrixs.h"
#include "Richardson.h"

#ifndef PROJECT_ROOT_DIR
#define PROJECT_ROOT_DIR "."
#endif

static int openOutput(const char* relativePath, char* path, size_t pathSize, FILE** output) {
    snprintf(path, pathSize, "%s/%s", PROJECT_ROOT_DIR, relativePath);
    *output = fopen(path, "w");
    if (*output == NULL) {
        perror(path);
        return 0;
    }
    return 1;
}

int main(void) {
    printf("Hello, World!\n");
    char path[1024];
    FILE* output = NULL;

    int m = 3;
    int N = 10;

    srand(1);
    Matrix* x_fix = createRandomMatrix(N,1, -5.,5.);
    Matrix* Q = createQMatrix(N);
    Matrix* D = builtDMatrix(N, 1);




    return 0;
}

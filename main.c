#include <stdio.h>

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

    
    return 0;
}

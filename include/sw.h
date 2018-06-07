
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
/**/
#define MAX_LEN (pow(2,6))
#define MATCH (5)
#define MISMATCH (-5)
#define W (10)

#define GRIDX 8
#define GRIDY 8

/**/
int readFile(char*, char*, int);

void displayDPMatrix(int*, int, int);

void initDPMatrix(int*, int , int);

void align(int*, int, int, char*, char*);

void traceBack(int*, int, int, char*, char*);

void computeBlock(int *dpM, char *strA, char *strB, int nRow, int nCol, int ox, int oy, int blockx, int blocky);

void align_omp(int*, int, int, char*, char*);

void traceBack_omp(int*, int, int, char*, char*);


#include "sw.h"

int main(int argc, char** argv)
{	if(argc != 2) {
		printf("ussage: %s log2(max_len)\n", argv[0]);
		exit(0);
	}
	int max_len = pow(2,atoi(argv[1]));
	char *strA = (char *) malloc(max_len*sizeof(char)+1);
	char *strB = (char *) malloc(max_len*sizeof(char)+1);
	int nRow, nCol, lenA, lenB;
	struct timeval start, end;
	double execTime;
	lenA = readFile("/home/canhld/dcsWsp/term_project/Smith-Waterman-MPI/data/inputA.txt",strA, max_len);
	lenB = readFile("/home/canhld/dcsWsp/term_project/Smith-Waterman-MPI/data/inputB.txt",strB, max_len);
	printf("Read file done\n");	
	nRow = lenA + 1;
	nCol = lenB + 1;
	int* dpM = (int*) malloc(sizeof(int)*nRow*nCol);
	initDPMatrix(dpM, nRow, nCol);
	// displayDPMatrix(dpM, nRow, nCol);
	printf("init matrix done\n");
	/**/
    gettimeofday(&start, NULL);

	align(dpM, nRow, nCol, strA, strB);

//	displayDPMatrix(dpM, nRow, nCol);
	traceBack(dpM, nRow, nCol, strA, strB);

   	gettimeofday(&end, NULL);
	execTime = (end.tv_sec - start.tv_sec) * 1000.0;      // sec to ms
    execTime += (end.tv_usec - start.tv_usec) / 1000.0;   // us to ms
	printf("Sequential code execution time: %lf\n", execTime);

	// memset(dpM, 0, nRow*nCol*sizeof(int));

    gettimeofday(&start, NULL);

	align_omp(dpM, nRow, nCol, strA, strB);

	// displayDPMatrix(dpM, nRow, nCol);
	traceBack_omp(dpM, nRow, nCol, strA, strB);

   	gettimeofday(&end, NULL);
	/**/
	execTime = (end.tv_sec - start.tv_sec) * 1000.0;      // sec to ms
    execTime += (end.tv_usec - start.tv_usec) / 1000.0;   // us to ms
	printf("Parallel code execution time: %lf\n", execTime);
	free(dpM);
	free(strA);
	free(strB);
	return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
/**/
#define MATCH (5)
#define MISMATCH (-5)
#define W (10)

#define MPI01 (0)
#define MPI12 (1)
#define MPI23 (3)

#define LEN ((int) pow(2,15))
#define ROW (LEN+1)
#define COL (LEN+1)
// #define nProc (2)
#define RLOCAL (ROW/nProc+1)
#define RBLOCK RLOCAL
#define CBLOCK RLOCAL

/**/
int readFile(char*, char*, int);

void displayDPMatrix(int*, int, int);

void initDPMatrix(int*, int , int);

void align(int*, int, int, char*, char*);

void traceBack(int*, int, int, char*, char*);

void computeBlock(int*, char*, char*,int, int, int, int, int, int);

void align_omp(int*, int, int, char*, char*);

void traceBack_omp(int*, int, int, char*, char*);

int main(int argc, char** argv)
{
	char *strA = (char *) malloc(LEN*sizeof(char)+1);
	char *strB = (char *) malloc(LEN*sizeof(char)+1);
	int *dpM, *sub_dpM;
	int nProc, rank;
	double start, end;
	double execTime = 0;
	MPI_Status status;
    MPI_Request request;

//	MPI intit
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    // if(nP != nProc){
    //     printf("Number of nodes must be %d\n", nProc);
    //     exit(0);
    // }
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Type_vector(RLOCAL, RLOCAL, LEN+1, MPI_INT, &block);
	// MPI_Type_commit(&block);
	// MPI_Type_vector(RLOCAL, 1, LEN+1, MPI_INT, &column);
	// MPI_Type_commit(&column);
//	master node will read sequences and broadcast it to all node
	if(rank == 0){
		readFile("inputA.txt",strA, LEN);
		readFile("inputB.txt",strB, LEN);
		printf("Read file done\n");	
	}
	MPI_Bcast(strA, LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(strB, LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
	if(rank == 0){ // init matrix and compute the first block
		dpM = (int*) malloc(sizeof(int)*ROW*COL);
		sub_dpM = dpM;
		initDPMatrix(dpM, ROW, COL);
		start = MPI_Wtime();
		computeBlock(dpM, strA, strB, ROW, COL, 0, 0, CBLOCK, RBLOCK);	
	}
	else{
		sub_dpM = (int *) malloc(sizeof(int)*RLOCAL*COL);
		for(int i = 0; i < RLOCAL; ++i){
			sub_dpM[i*COL] = 0;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	/*Send row to 2nd node*/
	for(int k = 1; k < 2*nProc-1; ++k){
		if(rank == 0) {
			// printf("====== ROUND %d=======\n", k);
			MPI_Send(sub_dpM + (RLOCAL-1)*COL, COL, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
			// printf("Node %d: sent to node %d in round %d!\n", rank, rank + 1, k);
		}
		else if(rank == nProc-1){
			MPI_Recv(sub_dpM, COL, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
			// printf("Node %d: recv from node %d in round %d!\n", rank, rank -1, k);
		}
		else {
			MPI_Recv(sub_dpM, COL, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
			// printf("Node %d: recv from node %d in round %d!\n", rank, rank -1, k);
			MPI_Send(sub_dpM + (RLOCAL-1)*COL, COL, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
			// printf("Node %d: sent to node %d in round %d!\n", rank, rank + 1, k);
			// if (rank == 1) {
			// 	displayDPMatrix(sub_dpM, RLOCAL, COL);
			// }
		}
		if(k < nProc){
			if(rank <= k){
					computeBlock(sub_dpM, strA + rank*(RLOCAL-1), strB, RLOCAL, COL, (k-rank)*(CBLOCK-1), 0, RBLOCK, CBLOCK);
					// if (rank == 1) {
					// 	printf("Node %d: compute done in round %d!\n", rank, k);
					// 	// displayDPMatrix(sub_dpM, RLOCAL, COL);
					// }
			}
		}
		else {
			if(rank > k - nProc){
					computeBlock(sub_dpM, strA + rank*(RLOCAL-1), strB, RLOCAL, COL, (k-rank)*(CBLOCK-1), 0, RBLOCK, CBLOCK);
					// if (rank == 1) {
					// 	printf("Node %d: compute done in round %d!\n", rank, k);
					// 	// displayDPMatrix(sub_dpM, RLOCAL, COL);
					// }
			}
		}	
		// MPI_Barrier(MPI_COMM_WORLD);
	}
	// for(int k = 1; k < nProc; ++k){

	// }
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Gather(sub_dpM+COL, (RLOCAL-1)*COL, MPI_INT, dpM + COL, (RLOCAL-1)*COL, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank != 0)
		MPI_Send(sub_dpM, RLOCAL*COL, MPI_INT, 0, 1, MPI_COMM_WORLD);
	else{
		for(int id = 1; id < nProc; ++id){
			MPI_Recv(dpM+id*COL*(RLOCAL-1), RLOCAL*COL, MPI_INT, id, 1, MPI_COMM_WORLD, &status);
		}
		printf("Data gather at node %d!\n", rank);
		// computeBlock(dpM, strA, strB, (LEN+1), (LEN+1), RLOCAL - 1, RLOCAL - 1, RLOCAL, RLOCAL);
		// displayDPMatrix(dpM, ROW, COL);
		traceBack_omp(dpM, ROW, COL, strA, strB);

		end = MPI_Wtime();
		//execTime = (end.tv_sec - start.tv_sec) * 1000.0;      // sec to ms
		//execTime += (end.tv_usec - start.tv_usec) / 1000.0;   // us to ms
		printf("MPI code execution time: %lf\n", (end - start)*1000);
	}

    // gettimeofday(&start, NULL);

	// // displayDPMatrix(dpM, (LEN+1), (LEN+1));
	// traceBack(dpM, (LEN+1), (LEN+1), strA, strB);

   	// gettimeofday(&end, NULL);
	// execTime = (end.tv_sec - start.tv_sec) * 1000.0;      // sec to ms
    // execTime += (end.tv_usec - start.tv_usec) / 1000.0;   // us to ms
	// printf("Sequential code execution time: %lf\n", execTime);
	MPI_Finalize();
	if(rank == 0)
	free(dpM);
	else free(sub_dpM);
	free(strA);
	free(strB);
	return 0;
}

void initDPMatrix(int *dpM, int nRow, int nCol)
{
	int i;
	int j;
	memset(dpM, 0, 4*nRow*nCol);
	for(i=0;i<nRow;++i) dpM[i*nCol] = 0;
	for(j=0;j<nCol;++j) dpM[j] = 0;
	/**/	
}
/**/
int readFile(char* strInputFileNameA, char *strA, int max_len)
{
	int c;
	int i = 0, lenA;
	FILE* fIn;
	fIn = fopen(strInputFileNameA,"r");
	if(fIn == NULL)
	{
		printf("Can not open file %s\n",strInputFileNameA);
		exit(1);
	}
	/**/
	while((c = fgetc(fIn)) != EOF && i < max_len){  
		if(c == 'a' || c == 't' || c == 'g' || c == 'c') strA[i++] = c;
	}
	strA[i] = '\0';
	lenA = strlen(strA);
	fclose(fIn);	
	printf("len = %d\n",lenA);
	return lenA;
}
void displayDPMatrix(int* dpM, int nRow, int nCol)
{
	int i,j;
	for(i=0;i<nRow;++i)
	{
		for(j=0;j<nCol;++j)
			printf("%d\t",dpM[i*nCol+j]);
		printf("\n");
	}
}

void align(int *dpM, int nRow, int nCol, char *strA, char *strB)
{
	int tmpM, tmpI, tmpD,score;
	int i,j;
	for(i=1;i<nRow;++i)
		for(j=1;j<nCol;++j)
		{
			if(strA[i-1] == strB[j-1]) tmpM = dpM[(i-1)*nCol+(j-1)] + MATCH;
			else tmpM = dpM[(i-1)*nCol+(j-1)] + MISMATCH; 
			tmpD = dpM[(i-1)*nCol+j] - W;
			tmpI = dpM[i*nCol+(j-1)] - W;
			score = tmpM > tmpD ? tmpM:tmpD;
			score = score > tmpI ? score : tmpI;
			score = score > 0 ? score : 0;
			dpM[i*nCol+j] = score;
		}
}

void traceBack(int *dpM, int nRow, int nCol, char *strA, char *strB)
{
	int i,j,i1,j1;
	int maxCol, maxRow, maxScore;
	int traceType; //1: diagonal (match/mismatch), 2: left (insert to A), 3: up (insert to B)
	int seqLen;
	char *seqA, *seqB;
	seqA = (char *) malloc((nRow - 1)*sizeof(char));
	seqB = (char *) malloc((nCol - 1)*sizeof(char));
	maxCol = maxRow = maxScore = 0;
	/**/
	for(i=0;i<nRow;++i)
	{
		for(j=0;j<nCol;++j)
			if(dpM[i*nCol+j] > maxScore)
			{
				maxCol = j;
				maxRow = i;
				maxScore = dpM[i*nCol+j]; 
			} 
	}
	/**/
	printf("Max Score = %d, at row %d, col %d\n",maxScore,maxRow,maxCol);
	/**/
	i = maxRow;
	j = maxCol;
	seqLen = 0;
	while  (((i>0) || (j>0))  &&  (dpM[i*nCol+j]!=0))
	{
		i1 = i - 1;
		j1 = j - 1;
		if(dpM[i*nCol+j] == dpM[i1*nCol+j1] + MATCH) traceType = 1;
		if(dpM[i*nCol+j] == dpM[i1*nCol+j1] + MISMATCH) traceType = 1;
		if(dpM[i*nCol+j] == dpM[i*nCol+j1] - W) traceType = 2;
		if(dpM[i*nCol+j] == dpM[i1*nCol+j]  - W) traceType = 3;
		/**/
		switch (traceType)
		{
			case 1:
				seqA[seqLen] = strA[i-1];
				seqB[seqLen] = strB[j-1];
				i--; j--;
				break;
			case 2:
				seqA[seqLen] = '-';
				seqB[seqLen] = strB[j-1];
				j--;
				break;
			case 3:
				seqA[seqLen] = strA[i-1];
				seqB[seqLen] = '-';
				i--;
				break;
			default: exit(1);
		}
		seqLen++;
		/**/  
	}//end of while
	printf("Sequence Result:\n");
	for(i=seqLen-1;i>=0;--i) printf("%c",seqA[i]);
	printf("\n");
	for(i=seqLen-1;i>=0;--i) printf("%c",seqB[i]);
	printf("\n");
	free(seqA);
	free(seqB);
	return;
}


void computeBlock(int *dpM, char *strA, char *strB, int nRow, int nCol, int ox, int oy, int blockx, int blocky){
	int x, y, i, j;
	int tmpM, tmpI, tmpD, score;
	for(y=1;y<blocky;++y)
		for(x=1;x<blockx;++x){
			i = y + oy;
			j = x + ox;
			if(strA[i-1] == strB[j-1]) tmpM = dpM[(i-1)*nCol+(j-1)] + MATCH;
			else tmpM = dpM[(i-1)*nCol+(j-1)] + MISMATCH; 
			tmpD = dpM[(i-1)*nCol+j] - W;
			tmpI = dpM[i*nCol+(j-1)] - W;
			score = tmpM > tmpD ? tmpM:tmpD;
			score = score > tmpI ? score : tmpI;
			score = score > 0 ? score : 0;
			dpM[i*nCol+j] = score;
			// printf("i : %d j: %d\n", i, j);
		}

}

void traceBack_omp(int *dpM, int nRow, int nCol, char *strA, char *strB)
{
	int i,j,i1,j1;
	int maxCol, maxRow, maxScore =  -1, p;
	int traceType; //1: diagonal (match/mismatch), 2: left (insert to A), 3: up (insert to B)
	int seqLen;
	char *seqA, *seqB;
	seqA = (char *) malloc((nRow - 1)*sizeof(char));
	seqB = (char *) malloc((nCol - 1)*sizeof(char));
	maxCol = maxRow = maxScore = 0;
	/**/

	#pragma omp parallel 
	{
		int index_local = 0;
		int max_local = maxScore;
	#pragma omp for nowait
	for(i=0;i<nRow*nCol;++i)
			if(dpM[i] > max_local)
			{
				index_local = i;
				max_local = dpM[i]; 
			} 
	#pragma omp critical
	{
		if (max_local > maxScore){
			maxScore = max_local;
			p = index_local;
		}
	}
	}
	maxCol = p%nCol;
	maxRow = p/nCol;
	/**/
	printf("Max Score = %d, at row %d, col %d\n",maxScore,maxRow,maxCol);
	/**/
	i = maxRow;
	j = maxCol;
	seqLen = 0;
	while  (((i>0) || (j>0))  &&  (dpM[i*nCol+j]!=0))
	{
		i1 = i - 1;
		j1 = j - 1;
		if(dpM[i*nCol+j] == dpM[i1*nCol+j1] + MATCH) traceType = 1;
		if(dpM[i*nCol+j] == dpM[i1*nCol+j1] + MISMATCH) traceType = 1;
		if(dpM[i*nCol+j] == dpM[i*nCol+j1] - W) traceType = 2;
		if(dpM[i*nCol+j] == dpM[i1*nCol+j]  - W) traceType = 3;
		/**/
		switch (traceType)
		{
			case 1:
				seqA[seqLen] = strA[i-1];
				seqB[seqLen] = strB[j-1];
				i--; j--;
				break;
			case 2:
				seqA[seqLen] = '-';
				seqB[seqLen] = strB[j-1];
				j--;
				break;
			case 3:
				seqA[seqLen] = strA[i-1];
				seqB[seqLen] = '-';
				i--;
				break;
			default: exit(1);
		}
		seqLen++;
		/**/  
	}//end of while
	printf("Aligned Sub Sequence Result:\n");
	for(i=seqLen-1;i>=0;--i) printf("%c",seqA[i]);
	printf("\n");
	for(i=seqLen-1;i>=0;--i) printf("%c",seqB[i]);
	printf("\n");
	free(seqA);
	free(seqB);
	return;
}

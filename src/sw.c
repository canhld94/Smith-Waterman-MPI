#include "sw.h"

void initDPMatrix(int *dpM, int nRow, int nCol)
{
	int i;
	int j;
	// memset(dpM, -1, 4*(nRow*nCol-1));
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

void align_omp(int *dpM, int nRow, int nCol, char *strA, char *strB)
{
	int i, j, k;
	int blockx = nCol / GRIDX + 1, blocky = nRow / GRIDY + 1;
	for(k = 0; k < GRIDX; ++k){
	#pragma omp parallel for num_threads(2*GRIDX)  private(i, j) shared(dpM)
		for(j = 0; j <= k ; ++j){
			i = k - j;
			int ox = i*(blockx - 1);
			int oy = j*(blocky - 1);
			computeBlock(dpM, strA, strB,nRow, nCol,ox, oy, blockx, blocky);
		}
	}
	// displayDPMatrix(dpM, nRow, nCol);
	// printf("Region 1 done \n");
	for(k = 1 ; k < GRIDY; ++k){
	#pragma omp parallel for num_threads(2*GRIDX) private(i, j) shared(dpM)
		// #pragma omp for schedule(dynamic,1)
		for(j = k; j < GRIDY; ++j){
			i = GRIDY + k  - j - 1;
			int ox = i*(blockx - 1);
			int oy = j*(blocky - 1);
			computeBlock(dpM, strA, strB,nRow, nCol,ox, oy, blockx, blocky);
		}
	}
	// displayDPMatrix(dpM, nRow, nCol);
	// printf("Region 2 done \n");
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
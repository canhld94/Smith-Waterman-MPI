
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
/**/
#define MAX_LEN (400)
#define MATCH (5)
#define MISMATCH (-5)
#define W (5)
/**/
char *strA, *strB;
int lenA, lenB;
int nRow, nCol;
int* dpM;
/**/
void readFile(char* strInputFileNameA, char *strInputFilenameB);
/**/
void initDPMatrix();
/**/
void align();
/**/
void displayDPMatrix();
/**/
void traceBack();
int main(int argc, char** argv)
{
	strA = (char *) malloc(MAX_LEN*sizeof(char));
	strB = (char *) malloc(MAX_LEN*sizeof(char));
	readFile("inputA.txt","inputB.txt" );
	printf("Read file done\n");	
	nRow = lenA + 1;
	nCol = lenB + 1;
	dpM = (int*) malloc(sizeof(int)*nRow*nCol);
	initDPMatrix();
	printf("init matrix done\n");
	/**/
	align();
	/**/
//	displayDPMatrix();
	/**/
	traceBack();
	/**/
	free(dpM);
	free(strA);
	free(strB);
	return 0;
}
/**/
void initDPMatrix()
{
	int i;
	int j;
	for(i=0;i<nRow;++i) dpM[i*nCol] = 0;
	for(j=0;j<nCol;++j) dpM[j] = 0;
	/**/	
}
/**/
void align()
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
/**/
void readFile(char* strInputFileNameA, char *strInputFileNameB)
{
	int c;
	int i = 0, j = 0; 
	FILE* fIn;
	fIn = fopen(strInputFileNameA,"r");
	if(fIn == NULL)
	{
		printf("Can not open file %s\n",strInputFileNameA);
		exit(1);
	}
	/**/
	while((c = fgetc(fIn)) != EOF && i < MAX_LEN){
		if(c == 'a' || c == 't' || c == 'g' || c == 'c') strA[i++] = c;
	}
	strA[i] = '\0';
	lenA = strlen(strA);
	fclose(fIn);	
        fIn = fopen(strInputFileNameB,"r");
        if(fIn == NULL)
        {
                printf("Can not open file %s\n",strInputFileNameB);
                exit(1);
        }
        /**/
        while((c = fgetc(fIn)) != EOF && j < MAX_LEN){
		if(c == 'a' || c == 't' || c == 'g' || c == 'c') strB[j++] = c;
        }
	strB[j] = '\0';
	lenB = strlen(strB);
	fclose(fIn);
	printf("A: %s\n",strA);
	printf("B: %s\n",strB);
	printf("lenA = %d, lenB = %d\n",lenA,lenB);
	fflush(stdout);
}
void displayDPMatrix()
{
	int i,j;
	for(i=0;i<nRow;++i)
	{
		for(j=0;j<nCol;++j)
			printf("%d ",dpM[i*nCol+j]);
		printf("\n");
	}
}
/**/
void traceBack()
{
	int i,j,i1,j1;
	int maxCol, maxRow, maxScore;
	int traceType; //1: diagonal (match/mismatch), 2: left (insert to A), 3: up (insert to B)
	int seqLen;
	char *seqA, *seqB;
	seqA = (char *) malloc(lenA*sizeof(char));
	seqB = (char *) malloc(lenB*sizeof(char));
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
}

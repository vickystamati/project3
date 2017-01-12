#include <stdio.h>
#include <stdlib.h>
#include "rmsd.h"
#include "drmsd.h"
#include <time.h>
#define bufSize 2048

int main(int argc, char *argv[])
{

	FILE* fp;
	char buflen[bufSize];
	int i,j,tablescounter,rows,choice;
	double ***input,***inputbig;
	fp = fopen((argv[1]), "r");
	fscanf(fp, "%s", buflen);
	tablescounter=atoi(buflen);
	fscanf(fp, "%s", buflen);
	rows=atoi(buflen);
	printf("%d %d\n",tablescounter,rows);
	input=malloc(tablescounter*sizeof(double**));
	inputbig=malloc(tablescounter*sizeof(double**));
	for(i=0;i<tablescounter;i++)
	{
		input[i]=malloc(rows *sizeof(double*));
		inputbig[i]=malloc(rows *sizeof(double*));
		for(j=0;j<rows;j++)
		{
			input[i][j]=malloc(3*sizeof(double));
			inputbig[i][j]=malloc(3*sizeof(double));
		}
	}
	readfile(fp,tablescounter,rows,input,inputbig);
	fclose(fp);
 	//crmsd(tablescounter,rows,inputbig,0,3);

	//return 0;
	printf("Press 0 for c-RMSD or 1 for d-RMSD\n");
	scanf("%d",&choice);
	if(choice==0)
		crmsdmethod(input,inputbig,tablescounter,rows);
	else if(choice==1)
		drmsdmethod(input,tablescounter,rows);
	else
		printf("WRONG CHOICE!\n");
	for(i=0;i<tablescounter;i++)
	{
		for(j=0;j<rows;j++)
		{
			free(input[i][j]);
			free(inputbig[i][j]);
		}
		free(input[i]);
		free(inputbig[i]);
	}
	free(input);
	free(inputbig);
}

#include <stdio.h>
#include <stdlib.h>
#include "rmsd.h"
#include "drmsd.h"
#include <time.h>
#include "clrec.h"
#include <string.h>
#include "nnlsh.h"
#include <assert.h>
#define bufSize 2048

int main(int argc, char *argv[])
{
	if(argc!=2 && argc!=5 && argc!=6)
	{
		printf("Wrong arguments!\n");
		return -1;
	}
	FILE* fp;
	char buflen[bufSize];
	int i,j,tablescounter,rows,choice,val;
	double ***input,***inputbig;
	if(argc==2)
		choice=1;
	else 
		choice=0;
	if(argc==6)
		val=1;
	else
		val=0;
	if(choice==1)
	{
		printf("Press 0 for c-RMSD or 1 for d-RMSD\n");
		scanf("%d",&choice);
		if(choice==0)
		{
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
			assert(crmsdmethod(input,inputbig,tablescounter,rows));
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
		else if(choice==1)
		{
			fp = fopen((argv[1]), "r");
			fscanf(fp, "%s", buflen);
			tablescounter=atoi(buflen);
			fscanf(fp, "%s", buflen);
			rows=atoi(buflen);
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
			assert(drmsdmethod(input,tablescounter,rows));
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
		else
			printf("WRONG CHOICE!\n");
		
	}
	else if(choice==0)
	{
		printf("Press 0 for NN-LSH or 1 for Clustering\n");
		scanf("%d",&choice);
		if(choice==0)
			assert(nnlsh(argv[2],argv[4],val));
		else if(choice==1)
			assert(clusrec(argv[2],argv[4],val));
		else
			printf("WRONG CHOICE!\n");
	}
	else
		printf("WRONG CHOICE!\n");
}

#include <stdio.h>
#include <stdlib.h>
#include "list.h"
#include "rmsd.h"
#define bufSize 2048

int main(int argc, char *argv[])
{
	int z;
	struct node * temp;
	struct list * lista;
	int i,j;
	FILE* fp;
	int tablescounter;
	int rows;
	int less;
	double values[3];
	char buflen[bufSize],buflen1[bufSize],buflen2[bufSize],buflen3[bufSize];
	fp = fopen((argv[1]), "r");
	fscanf(fp, "%s", buflen);
	tablescounter=atoi(buflen);
	fscanf(fp, "%s", buflen);
	rows=atoi(buflen);
	printf("%d %d\n",tablescounter,rows);
	less=(rows*(rows-1))/2;
	lista=malloc(tablescounter*sizeof(struct list));
	for(i=0;i<tablescounter;i++)
	{
		lista[i].distance=malloc(rows*sizeof(double));
		lista[i].data=malloc(rows*sizeof(struct node));
		for(j=0;j<rows;j++)
		{
			lista[i].data[j].key1=malloc(3*sizeof(double));
		}
	}
	i=0;
	j=0;
	while(j<5)
	{
		i=0;
		while(i<rows)
		{
			fscanf(fp, "%s %s %s", buflen1,buflen2,buflen3);
			values[0]=atof(buflen1);
			values[1]=atof(buflen2);
			values[2]=atof(buflen3);
			inserteuclidean(lista,values,j,i);
			i++;
		}
		j++;
	}
	printf("Edwwww %f\n",lista[1].data[0].key1[1]);

	i=0;
	while(i<2)
	{
		j=0;
		
		while(j<rows)
		{
			temp=&lista[i].data[j];
			for(z=0;z<3;z++)
			{
				printf("%f ",temp->key1[z]);
			}
			printf(" \n ");
			j++;
		}
		printf("dsadas\n");
		i++;
	}
	sleep(4);
	/*
	distancefind(lista,tablescounter,rows);
		printf("Edwwww2\n");
	*/
	createnewc(lista,rows,tablescounter);
	fclose(fp);
	return 0;
}

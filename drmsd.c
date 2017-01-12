#include <stdio.h>
#include <stdlib.h>
#include "rmsd.h"
#include <math.h>
#include "drmsd.h"
#include <time.h>
#define bufSize 2048

int drmsdmethod(double ***input,int tablescounter,int rows)
{
	int i,j,r[3],loop,size,ktab[7],cent=5;
	double **firstm,**disttable,silh[7],time[7];
	r[0]=rows;
	r[1]=(int)pow(rows,1.5);
	r[2]=(rows*(rows-1))/2;
	firstm=malloc(rows*sizeof(double*));
	for(i=0;i<rows;i++)
		firstm[i]=malloc(rows*sizeof(double));
	findfirstmatr(input,firstm,tablescounter,rows);
	if(rows%2==0)
		size=((rows*rows)/2) - (rows/2) -1;
	else
		size=((rows*rows)/2) - (rows/2) ;
	disttable=malloc(size*sizeof(double*));
	printf("size= %d\n",size);
	for(i=0;i<size;i++)
		disttable[i]=malloc(3*sizeof(double));
	distancesort(firstm,disttable,size,rows);
//loop gia ta diaforetika k
//desmeysi kentroieidous analoga me k clist[k]
	for(i=0;i<3;i++)
	{
		clustering(input,disttable,tablescounter,rows,r[i],ktab,silh,time,cent);
	}
	

	for(i=0;i<rows;i++)
		free(firstm[i]);
	free(firstm);
	return 1;
}

int findfirstmatr(double *** input,double ** firstm,int tablescounter,int rows)
{
	int i,j;
	for(i=0;i<rows;i++)
	{
		for(j=i;j<rows;j++)
		{
			if(i==j)
				firstm[i][j]=0;
			else
				firstm[i][j]=eucldistance(input,0,i,j);
		}
	}
	for(i=0;i<rows;i++)
	{
		for(j=i;j<rows;j++)
		{
			if(i==j)
				firstm[i][j]=0;
			else
				firstm[j][i]=firstm[i][j];
		}
	}
	return 1;
}

double eucldistance(double *** input,int d,int pos1,int pos2)
{
	double distance,sum,sum1,sum2;
	int j,z;
	distance=0;
	sum=0;
	sum2=0;
	for(z=0;z<3;z++)//typos euclidean
	{
		sum=input[d][pos1][z] - input[d][pos2][z];
		sum=sum*sum;
		sum2+=sum;
		sum=0;
	}
	distance=sqrt(sum2);//upologismos apostasis
	return distance;
}
int distancesort(double **firstm,double**disttable,int size,int rows)
{
	int i,j,v;
	v=0;
	for(i=0;i<rows;i++)
	{
		for(j=i;j<rows;j++)
		{
			if(i<j)
			{
				disttable[v][0]=firstm[i][j];
				disttable[v][1]=i;
				disttable[v][2]=j;		
				v++;
			}
		}
	}
	bubble_sort2d_3cell(disttable, v);
	return 1;

}


void bubble_sort2d_3cell(double** matr, int counter)
{
	double t1,t2,t3;
	int i,j;
	for (i=0;i<(counter-1);i++)
	{
		for (j=0;j<counter-i-1;j++)
		{
			if (matr[j][0] > matr[j+1][0])
			{
				t1 = matr[j][1];
				t2= matr[j][0];
				t3= matr[j][2];
				matr[j][1] = matr[j+1][1];
				matr[j][0] = matr[j+1][0];
				matr[j][2] =matr[j+1][2];
				matr[j+1][0] = t2;
				matr[j+1][1] = t1;
				matr[j+1][2] = t3;
			}
		}
	}
}



int clustering(double ***input,double ** disttable,int tablescounter,int rows,int num,int *ktab,double *silh,double *time,int cent)//num=r
{
	int i,j,z;
	struct clustlist * clist;
	double ** matrix,**drmsdmatr,sum,sum2;
	clist=malloc(cent*sizeof(struct clustlist));
	for(i=0;i<cent;i++)
	{
		clist[i].key=malloc(rows*sizeof(double*));
		for(j=0;j<rows;j++)
			clist[i].key[j]=malloc(3*sizeof(double*));
		clist[i].head=NULL;
	}
	matrix=malloc(tablescounter*sizeof(double*));
	for(i=0;i<tablescounter;i++)
		matrix[i]=malloc(num*sizeof(double));

	drmsdmatr=malloc(tablescounter*sizeof(double*));
	for(i=0;i<tablescounter;i++)
		drmsdmatr[i]=malloc(tablescounter*sizeof(double));
	
	for(i=0;i<tablescounter;i++)
	{
		for(j=0;j<num;j++)
		{
			if(i==0)
				matrix[i][j]=disttable[j][0];
			else
				matrix[i][j]=eucldistance(input,i,disttable[j][1],disttable[j][2]);
		}
	
	}
	//k loops tou means
	kmeansstart(clist,input,matrix,tablescounter,rows,cent);
	
	for(i=0;i<cent;i++)
	{
		for(j=0;j<rows;j++)
			free(clist[i].key[j]);
		free(clist[i].key);
	}
	free(clist);
}

int meanspam(double *** input,struct clustlist * clist,int tablescounter,int rows,int cent,int *matr,int nnum)
{
	int i,j,k,flag,loop;
	double **pammatr;
	if(nnum==0)
		loop=tablescounter;
	else
		loop=nnum;
	pammatr=malloc(cent*sizeof(double*));
	for(i=0;i<cent;i++){
		pammatr[i]=malloc(2*sizeof(double));
	}
	for(k=0;k<loop;k++)
	{
		flag=1;
		for(i=0;i<cent;i++)
		{	
			pammatr[i][0]=i;
			if(nnum==0)
			{
				if(k==clist[i].id)
				{
					flag=0;
					break;
				}
			}
			else
			{
				if(matr[k]==clist[i].id)
				{
					flag=0;
					break;
				}
			}
		}
		if(flag==1)
		{
			for(j=0;j<cent;j++)
			{
				if(nnum==0)
					pammatr[j][1]=crmsd(tablescounter,rows,inputbig,k,clist[j].id);
				else
					pammatr[j][1]=crmsd(tablescounter,rows,inputbig,matr[k],clist[j].id);
			}
			bubble_sort2d(pammatr, cent);
			if(nnum==0)
				insertassig(pammatr,clist,k,input,rows);	
			else
				insertassig(pammatr,clist,matr[k],input,rows);	
		}
	}
	for(i=0;i<cent;i++)
		free(pammatr[i]);
	free(pammatr);
	return 1;
}*/


int kmeansstart(struct clustlist * clist,double ***input,double ** matrix,int tablescounter,int rows,int cent)
{
	int i,j,z,flag,random,listmatr[cent];
	i=0;
	while(i<cent)
	{
		flag=0;
		random=1+ (rand() /( RAND_MAX + 1.0))*(tablescounter-1);
		if(i==0)
		{
			listmatr[i]=random;
		}
		else
		{
			flag=0;
			for(j=0;j<i;j++)
				if(listmatr[j]==random)
					flag=1;
			if(flag==1)
				continue;
			else
				listmatr[j]=random;
		}
		i++;
	}
	for(i=0;i<cent;i++)
	{
		for(j=0;<rows;j++)
			for(z=0;z<3;z++)
				clist[i].key[j][z]=input[listmatr[i]][j][z];
	}
	return 1;
}





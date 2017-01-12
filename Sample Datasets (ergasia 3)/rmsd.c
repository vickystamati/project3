#include <stdio.h>
#include <stdlib.h>
#include "list.h"
#include "rmsd.h"
//#include <gsl/gsl_blas.h>
#include "cblas.h"
#include "lapacke.h"
//#include <lapacke.h>

void createnewc(struct list * lista,int rows,int tablescounter)
{
	double xc[3];
	double yc[3];
	int i,j,z,v,t;
	double ** reverse;
	double *rev,*y,*last,*q;
	double *utable,*vt,*sbig,*temputable;
	double detq;
	lapack_int out;
	utable=LAPACKE_malloc(sizeof(double)*9);
	temputable=LAPACKE_malloc(sizeof(double)*9);
	vt=LAPACKE_malloc(sizeof(double)*9);
	sbig=LAPACKE_malloc(sizeof(double)*9);
	rev=malloc(3*rows*sizeof(double));
	y=malloc(3*rows*sizeof(double));
	last=malloc(3*3*sizeof(double));
	q=malloc(3*3*sizeof(double));
	reverse=malloc(3*sizeof(double*));
	for(z=0;z<3;z++)
	{
		reverse[z]=malloc(rows*sizeof(double));
	}
	i=0;
	struct node* temp;
	while(i<tablescounter)
	{
		for(z=0;z<3;z++)
		{
			xc[z]=0;
			yc[z]=0;	
		}
		j=0;
		printf("MOUNAKIAAAA %d\n",i);
		while(j<rows)
		{
			temp=&lista[i].data[j];
			for(z=0;z<3;z++)
			{
				printf("%f\n",temp->key1[z]);
				if(i%2==0)
					xc[z]+=temp->key1[z];
				else
					yc[z]+=temp->key1[z];
			}
				printf("MOUNAKIAAAA %d\n",j);
			j++;
		}
		printf("MOUNAKIAAAA\n");
		for(z=0;z<3;z++)
		{
			xc[z]=xc[z]/rows;
			yc[z]=yc[z]/rows;	
		}
				printf("MOUNAKIAAAA2\n");
		j=0;
		while(j<rows)
		{
			temp=&lista[i].data[j];
			for(z=0;z<3;z++)
			{
				if(tablescounter%2==0)
					temp->key1[z]-=xc[z];
				else
					temp->key1[z]-=yc[z];
			}
			j++;
		}
		printf("MOUNAKIAAAA3\n");
		for(z=0;z<rows;z++)
		{
			printf("TIMES %d\n",z);
			temp=&lista[i].data[z];
			for(v=0;v<3;v++)
			{
				reverse[v][z]=temp->key1[v];
				printf("%f \n",reverse[v][z]);
			}
		}
		sleep(2);
		printf("MOUNAKIAAAA4\n");
		for(v=0;v<3;v++)
		{
			for(z=0;z<rows;z++)
			{
				printf("%f ",reverse[v][z]);
			}
			printf(" \n ");
		}
		sleep(2);
		printf("MOUNAKIAAAA5\n");
		if(i%2==0)
		{
			t=0;
			for(v=0;v<3;v++)
			{
				for(z=0;z<rows;z++)
				{
					rev[t]=reverse[v][z];
					t++;
				}
			}
		}
		else
		{
			t=0;
			for(v=0;v<rows;v++)
			{
				temp=&lista[i].data[v];
				for(z=0;z<3;z++)
				{
					y[t]=temp->key1[z];
					t++;
				}
			}
			cblas_dgemm (CblasRowMajor, CblasNoTrans,CblasNoTrans,3,3,rows,1.0,rev,rows,y,3,0.0,last,3);
			for(z=0;z<9;z++){
			printf("skata %f\n",last[z]);
			
			}
			sleep(3);
			out=LAPACKE_dgesvj(LAPACK_ROW_MAJOR,'G','U','V',3,3,last,3,sbig,0,vt,3,utable);
			for(z=0;z<9;z++)
				printf("utable %f\n",utable[z]);
			for(z=0;z<9;z++)
				printf("sbig %f\n",sbig[z]);
			for(z=0;z<9;z++)
				printf("vt %f\n",vt[z]);
			sleep(2);
			if(sbig[2]<=0)//EDW THELW SBIG[8]
			{
				i++;
				continue;//termatise
			}
			cblas_dgemm (CblasRowMajor, CblasNoTrans,CblasNoTrans,3,3,3,1.0,vt,3,utable,3,0.0,q,3);
			detq=(q[0]*q[4]*q[8])+(q[1]*q[5]*q[6])+(q[2]*q[3]*q[7])-(q[2]*q[4]*q[6])-(q[1]*q[3]*q[8])-(q[0]*q[5]*q[7]);
			if(detq<0)
			{
				for(z=0;z<9;z++)
				{
					if(z==2 || z==5 || z==8)
						temputable[z]=-utable[z];
					else
						temputable[z]=utable[z];
				}
				cblas_dgemm (CblasRowMajor, CblasNoTrans,CblasNoTrans,3,3,3,1.0,vt,3,temputable,3,0.0,q,3);
			}
			
		}
		i++;	
	}
	

}

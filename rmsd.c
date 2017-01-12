#include <stdio.h>
#include <stdlib.h>
#include "rmsd.h"
#include <math.h>
#include "cblas.h"
#include "lapacke.h"
#include "drmsd.h"
#include <time.h>
#define bufSize 2048

int readfile(FILE *fp,int tablescounter,int rows,double ***input,double ***inputbig)
{
	int i,j,z;
	double cmatr[3];
	char buflen1[bufSize],buflen2[bufSize],buflen3[bufSize];
	for(i=0;i<tablescounter;i++)
	{
		for(j=0;j<3;j++)
			cmatr[j]=0;
		for(j=0;j<rows;j++)
		{
			fscanf(fp, "%s %s %s", buflen1,buflen2,buflen3);
			input[i][j][0]=atof(buflen1);
			input[i][j][1]=atof(buflen2);
			input[i][j][2]=atof(buflen3);
			for(z=0;z<3;z++)
				cmatr[z]+=input[i][j][z];
		}
		for(z=0;z<3;z++)
			cmatr[z]=cmatr[z]/rows;
		for(j=0;j<rows;j++)
			for(z=0;z<3;z++)
				inputbig[i][j][z]=input[i][j][z]-cmatr[z];
	}
	return 1;
}

double crmsd(int tablescounter,int rows,double ***inputbig,int pos1,int pos2)
{
	long double trace;
	double **matrix,**reverse,*x,*y,*xt,*utable,*vt,*temputable,*sbig,*q,detq,*mulxty,*mulxq,*final,dist,*res,*rest,*superb;
	int i,j,z,v;
	lapack_int out;
	matrix=malloc(rows *sizeof(double*));
	for(i=0;i<rows;i++)
		matrix[i]=malloc(3*sizeof(double));
	reverse=malloc(3*sizeof(double*));
	for(i=0;i<3;i++)
		reverse[i]=malloc(rows*sizeof(double));

	x=LAPACKE_malloc(3*rows*sizeof(double));
	xt=LAPACKE_malloc(3*rows*sizeof(double));
	y=LAPACKE_malloc(3*rows*sizeof(double));
	mulxty=LAPACKE_malloc(3*3*sizeof(double));
	utable=LAPACKE_malloc(3*3*sizeof(double));
	temputable=LAPACKE_malloc(sizeof(double)*9);
	vt=LAPACKE_malloc(sizeof(double)*9);
	sbig=LAPACKE_malloc(sizeof(double)*3);
	q=LAPACKE_malloc(3*3*sizeof(double));
	mulxq=LAPACKE_malloc(3*rows*sizeof(double));
	res=LAPACKE_malloc(3*rows*sizeof(double));
	rest=LAPACKE_malloc(3*rows*sizeof(double));
	final=LAPACKE_malloc(3*3*sizeof(double));
	superb=LAPACKE_malloc(9*sizeof(double));


	for(j=0;j<rows;j++)
		for(z=0;z<3;z++)
			reverse[z][j]=inputbig[pos1][j][z];
	v=0;
	for(j=0;j<rows;j++)
	{
		for(z=0;z<3;z++)
		{
			x[v]=inputbig[pos1][j][z];
			//printf("x %f\n",x[v]);
			v++;
		}
	}
	v=0;
	for(z=0;z<3;z++)
	{

		for(j=0;j<rows;j++)
		{
			xt[v]=reverse[z][j];
			//printf("xt %f\n",xt[v]);
			v++;
		}
	}	
	//printf("\n");
	v=0;
	for(j=0;j<rows;j++)
	{
		for(z=0;z<3;z++)
		{
			y[v]=inputbig[pos2][j][z];
			//printf("y %f\n",y[v]);
			v++;
		}
	}
	cblas_dgemm(CblasRowMajor, CblasNoTrans,CblasNoTrans,3,3,rows,1.0,xt,rows,y,3,0.0,mulxty,3);
	//out=LAPACKE_dgesvj(LAPACK_ROW_MAJOR,'G','U','V',3,3,mulxty,3,sbig,0,vreal,3,utable);

	out = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', 3, 3, mulxty, 3,sbig, utable, 3, vt, 3, superb );
	if(sbig[2]<=0)
	{
		printf("error with S\n");
		//sleep(2);
		return -1;//termatise
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
			q[z]=0;
		}
		cblas_dgemm (CblasRowMajor, CblasNoTrans,CblasNoTrans,3,3,3,1.0,temputable,3,vt,3,0.0,q,3);
	}
	cblas_dgemm (CblasRowMajor, CblasNoTrans,CblasNoTrans,rows,3,3,1.0,x,3,q,3,0.0,mulxq,3);
	for(j=0;j<rows*3;j++)
		mulxq[j]-=y[j];
	v=0;
	for(j=0;j<rows;j++)
	{
		for(z=0;z<3;z++)
		{
			matrix[j][z]=mulxq[v];
			v++;				
		}
	}
	for(j=0;j<rows;j++)
		for(z=0;z<3;z++)
			reverse[z][j]=matrix[j][z];

	v=0;
	for(j=0;j<rows;j++)
	{
		for(z=0;z<3;z++)
		{
			res[v]=matrix[j][z];
			v++;
		}
	}
	v=0;
	for(z=0;z<3;z++)
	{
		for(j=0;j<rows;j++)
		{
			rest[v]=reverse[z][j];
			v++;
		}
	}	
	cblas_dgemm (CblasRowMajor, CblasNoTrans,CblasNoTrans,3,3,rows,1.0,rest,rows,res,3,0.0,final,3);
	trace=final[0]+final[4]+final[8];
	if(trace<0)
		trace=-1*trace;
	trace=sqrt(trace);
	trace=trace/sqrt(rows);
	dist=trace;
	//printf("dist= %f\n",dist);
	LAPACKE_free(superb);
	LAPACKE_free(final);
	LAPACKE_free(rest);
	LAPACKE_free(res);
	LAPACKE_free(mulxq);
	LAPACKE_free(q);
	LAPACKE_free(sbig);
	LAPACKE_free(vt);
	LAPACKE_free(temputable);
	LAPACKE_free(utable);
	LAPACKE_free(mulxty);
	LAPACKE_free(y);
	LAPACKE_free(xt);
	LAPACKE_free(x);

	for(i=0;i<3;i++)
		free(reverse[i]);
	free(reverse);
	for(i=0;i<rows;i++)
		free(matrix[i]);
	free(matrix);

	return dist;
}


int claraf(struct clustlist * clist , double *** input,double *** inputbig,int rows,int tablescounter,int cent)
{
	int i,j,z,m,r,counter;
	int random,nnum,flag;
	int randmatr[cent];
	int *listmatr;
	double jval,jtemp;
	struct clustlist * clisttemp,*clistfinal;
	struct centlist * centtemp;
	srand(time(NULL));
	//struct node * temp;
	//struct list * lista=malloc(sizeof(struct list));
	//lista->head=NULL;
	nnum=20 + 2*cent;
	listmatr=malloc(nnum*sizeof(int));
	i=0;
	while(i<nnum)
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
	i=0;
	while(i<cent)
	{
		random=1+ (rand() /( RAND_MAX + 1.0))*(nnum-1);
		if(i==0)
		{
			randmatr[i]=random;
		}
		else
		{
			flag=0;
			for(j=0;j<i;j++)
				if(randmatr[j]==random)
					flag=1;
			if(flag==1)
				continue;
			else
				randmatr[j]=random;
		}
		clist[i].id=randmatr[i];
		for(z=0;z<rows;z++)
			for(m=0;m<3;m++)
				clist[i].key[z][m]=input[randmatr[i]][z][m];
		i++;
	}
	medpam(input,inputbig,clist,tablescounter,rows,cent,listmatr,nnum);
	jval=calcj(clist ,cent);
	clistfinal=malloc(cent*sizeof(struct clustlist));
	clisttemp=malloc(cent*sizeof(struct clustlist));
	for(i=0;i<cent;i++)
	{
		clistfinal[i].key=malloc(rows*sizeof(double*));
		clisttemp[i].key=malloc(rows*sizeof(double*));
		for(j=0;j<rows;j++)
		{
			clisttemp[i].key[j]=malloc(3*sizeof(double));
			clistfinal[i].key[j]=malloc(3*sizeof(double));
		}
		clisttemp[i].head=NULL;
		clistfinal[i].head=NULL;
	}
	counter=0;
	for(i=0;i<cent;i++)
	{
		for(j=0;j<cent;j++)
		{
			centtemp=clist[j].head;
			while(centtemp!=NULL)
			{	
				for(m=0;m<cent;m++)//pername ta kentroidi sto temp struct
				{
					if(i==m)
					{
						clisttemp[m].id=centtemp->id;
						for(z=0;z<rows;z++)
							for(r=0;r<3;r++)
								clisttemp[m].key[z][r]=centtemp->key[z][r];	
					}
					else
					{
						clisttemp[m].id=clist[m].id;
						for(z=0;z<rows;z++)
							for(r=0;r<3;r++)
								clisttemp[m].key[z][r]=clist[m].key[z][r];
					}
				}
				medpam(input,inputbig,clisttemp,tablescounter,rows,cent,listmatr,nnum);
				jtemp=calcj(clisttemp ,cent);
				if(jtemp<jval)
				{
					freeclustlist(clistfinal,cent,rows);
					swapclist(clistfinal,clisttemp,cent,rows);
					freeclustlist(clisttemp,cent,rows);
					jval=jtemp;
				}
				else
					freeclustlist(clisttemp,cent,rows);
				centtemp=centtemp->next;
				counter++;

			}
		}
	}
	if(clistfinal[0].head!=NULL)
	{
		freeclustlist(clist,cent,rows);
		swapclist(clist,clistfinal,cent,rows);
	}
	freeclustlist(clisttemp,cent,rows);
	freeclustlist(clistfinal,cent,rows);
	freeclustlist(clist,cent,rows);
		
	medpam(input,inputbig,clist,tablescounter,rows,cent,listmatr,0);

	for(i=0;i<cent;i++)
	{
		for(j=0;j<rows;j++)
		{
			free(clisttemp[i].key[j]);
			free(clistfinal[i].key[j]);
		}
		free(clisttemp[i].key);
		free(clistfinal[i].key);
	}

	free(clistfinal);
	free(clisttemp);
	free(listmatr);
	return 1;
}


int clara(struct clustlist * clist,double *** input,double *** inputbig,int rows,int tablescounter,int cent)
{
	struct clustlist*clistptr;
	int i,j,s=3;//apo theoria
	double jsum,jsumclara;
	clistptr=malloc(cent*sizeof(struct clustlist));
	for(i=0;i<cent;i++)
	{
		clistptr[i].key=malloc(rows*sizeof(double*));
		for(j=0;j<rows;j++)
			clistptr[i].key[j]=malloc(3*sizeof(double*));
		clistptr[i].head=NULL;
	}
	for(i=0;i<s;i++)
	{
		claraf(clistptr,input,inputbig,rows,tablescounter,cent);
		if(i==0)
		{
			jsumclara=calcj(clistptr ,cent );
			swapclist(clist,clistptr,cent,rows);
			freeclustlist(clistptr,cent,rows);
		}
		else
		{
			jsum=calcj(clistptr ,cent );
			if(jsum<jsumclara)
			{
				jsumclara=jsum;
				freeclustlist(clist,cent,rows);
				swapclist(clist,clistptr,cent,rows);
			}
			freeclustlist(clistptr,cent,rows);
		}
	}
	for(i=0;i<cent;i++)
	{

		for(j=0;j<rows;j++)
			free(clistptr[i].key[j]);
		free(clistptr[i].key);
	}
	free(clistptr);
	return 1;
}

double calcj(struct clustlist * clist,int cent)
{
	int i;
	double sum=0;
	struct centlist * temp;
	for(i=0;i<cent;i++)
	{
		temp=clist[i].head;
		while(temp!=NULL)
		{
			sum+=temp->distance;
			temp=temp->next;		
		}
	}
	return sum;
}


int freeclustlist(struct clustlist * clist,int cent,int rows)
{
	struct centlist * temp;
	int i,j;
	for(i=0;i<cent;i++)
	{
		temp=clist[i].head;	
		while (clist[i].head != NULL)
		{
			temp=clist[i].head;	
			clist[i].head=clist[i].head->next;
			for(j=0;j<rows;j++)		
				free(temp->key[j]);
			free(temp->key);
			free(temp);
		}
		clist[i].head=NULL;
	}
	return 1;
}

int swapclist(struct clustlist * clist,struct clustlist * tempclist,int cent,int rows)
{
	int i,j,z;
	struct centlist * temp, *temp2;
	for(i=0;i<cent;i++)
	{
		clist[i].id=tempclist[i].id;
		if(clist[i].key==NULL)
		{
			clist[i].key=malloc(rows*sizeof(double*));
			for(j=0;j<rows;j++)
				clist[i].key[j]=malloc(3*sizeof(double));
		}
		for(j=0;j<rows;j++)
			for(z=0;z<3;z++)
				clist[i].key[j][z]=tempclist[i].key[j][z];
		temp=tempclist[i].head;
		while(temp!=NULL)
		{
			if(clist[i].head==NULL)
			{
				clist[i].head=malloc(sizeof(struct centlist));
				clist[i].head->id=temp->id;
				clist[i].head->key=malloc(rows*sizeof(double*));
				for(j=0;j<rows;j++)
					clist[i].head->key[j]=malloc(3*sizeof(double));
				for(j=0;j<rows;j++)
					for(z=0;z<3;z++)
						clist[i].head->key[j][z]=temp->key[j][z];
				clist[i].head->distance=temp->distance;
				clist[i].head->listnearid=temp->listnearid;
				clist[i].head->distancenear=temp->distancenear;
				clist[i].head->next=NULL;
			}
			else
			{	
				temp2=clist[i].head;
				while(temp2->next!=NULL)
				{
					temp2=temp2->next;
				}
				temp2->next=malloc(sizeof(struct centlist));
				temp2->next->id=temp->id;
				temp2->next->key=malloc(rows*sizeof(double*));
				for(j=0;j<rows;j++)
					temp2->next->key[j]=malloc(3*sizeof(double));
				for(j=0;j<rows;j++)
					for(z=0;z<3;z++)
						temp2->next->key[j][z]=temp->key[j][z];
				temp2->next->distance=temp->distance;
				temp2->next->listnearid=temp->listnearid;
				temp2->next->distancenear=temp->distancenear;
				temp2->next->next=NULL;
			}			
			temp=temp->next;
		}
	}
	return 1;
}

int medpam(double *** input,double *** inputbig,struct clustlist * clist,int tablescounter,int rows,int cent,int *matr,int nnum)//medoids
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
}

int insertassig(double ** pammatr,struct clustlist * clist,int k,double ***input,int rows)
{
	int i,j,z,m;
	i=pammatr[0][0];
	if (clist[i].head == NULL)//an head null kane eisagwgi
	{
		clist[i].head=malloc(sizeof(struct centlist));
		clist[i].head->id=k;
		clist[i].head->key=malloc(rows * sizeof(double*));
		for(z=0;z<rows;z++)
			clist[i].head->key[z]=malloc(3*sizeof(double));
		for(z=0;z<rows;z++)
			for(m=0;m<3;m++)
				clist[i].head->key[z][m]=input[k][z][m];
		clist[i].head->distance=pammatr[0][1];
		j=pammatr[1][0];
		clist[i].head->listnearid=j;
		clist[i].head->distancenear=pammatr[1][1];
		clist[i].head->next=NULL;
	}
	else
	{
		struct centlist *current = clist[i].head;
		while (current->next != NULL) 
			current = current->next;
		current->next= malloc(sizeof(struct centlist));
		current->next->id=k;
		current->next->key=malloc(rows * sizeof(double*));
		for(z=0;z<rows;z++)
			current->next->key[z]=malloc(3*sizeof(double));
		for(z=0;z<rows;z++)
			for(m=0;m<3;m++)
				current->next->key[z][m]=input[k][z][m];
		current->next->distance=pammatr[0][1];
		j=pammatr[1][0];
		current->next->listnearid=j;
		current->next->distancenear=pammatr[1][1];
		current->next->next=NULL;
	}
	return 1;
}

void bubble_sort2d(double** matr, int counter)
{
	double t1,t2;
	int i,j;
	for (i=0;i<(counter-1);i++)
	{
		for (j=0;j<counter-i-1;j++)
		{
			if (matr[j][1] > matr[j+1][1])
			{
				t1 = matr[j][1];
				t2= matr[j][0];
				matr[j][1] = matr[j+1][1];
				matr[j][0] = matr[j+1][0];
				matr[j+1][0] = t2;
				matr[j+1][1] = t1;
			}
		}
	}
}

double silhouette(struct clustlist * clist,double ***inputbig,int cent,int rows,int tablescounter)//FILE * fp)
{
	int i;
	double suma,sumb,sumaver,sumall=0,counta,countb,countaver;
	struct centlist * temp;
	struct centlist *temp2;
	int flag;
	printf("Silhuette per cluster {");
	for(i=0;i<cent;i++)
	{
		temp=clist[i].head;
		sumaver=0;
		countaver=0;
		while(temp!=NULL)
		{
			flag=0;
			counta=0;
			suma=0;
			temp2=clist[i].head;
			while(temp2!=NULL)
			{
				if(temp2->id!=temp->id)
				{
					suma+=crmsd(tablescounter,rows,inputbig,temp->id,temp2->id);
					counta++;
				}
				temp2=temp2->next;
			}
			countb=0;
			sumb=0;
			temp2=clist[temp->listnearid].head;
			while(temp2!=NULL)
			{
				sumb+=crmsd(tablescounter,rows,inputbig,temp->id,temp2->id);
				countb++;
				temp2=temp2->next;
			}
			if(counta!=0 && countb!=0)
			{		
				suma=suma/counta;
				sumb=sumb/countb;
			}
			else
				flag=1;
			if(suma>sumb && suma!=0)
				sumaver+=(sumb/suma)-1;
			else if(suma<sumb && sumb!=0)
				sumaver+=1-(suma/sumb);
			else if(suma==sumb)
				sumaver+=0;
			else
				flag=1;
			temp=temp->next;
			if(flag==0)
				countaver++;
		}
		if(countaver!=0)
			sumaver=sumaver/countaver;
		else
			sumaver=0;
		if(i==cent-1)
			printf("%f}\n",sumaver);
		else
			printf("%f, ",sumaver);
		sumall+=sumaver;
	}
	if(sumall<0) sumall*=-1;
	sumall/=cent;
	return sumall;

}

void bubble_sort(int* matr, int counter)
{
	int i,j,t;
	for (i=0;i<(counter-1);i++)
	{
		for (j=0;j<counter-i-1;j++)
		{
			if (matr[j] > matr[j+1])
			{
				t = matr[j];
				matr[j] = matr[j+1];
				matr[j+1] = t;
			}
		}
	}
}

int crmsdmethod(double *** input,double *** inputbig,int tablescounter,int rows)
{
	struct clustlist * clist, *clistfinal;
	struct centlist *temp;
	int knum=20,cent,loop,bestcent,*matr,counter,i,j;
	double silh,bestsilh;
	for(loop=0;loop<1;loop++)
	{
		cent=1+ (rand() /( RAND_MAX + 1.0))*(knum-1);
		clist=malloc(cent*sizeof(struct clustlist));
		for(i=0;i<cent;i++)
		{
			clist[i].key=malloc(rows*sizeof(double*));
			for(j=0;j<rows;j++)
				clist[i].key[j]=malloc(3*sizeof(double*));
			clist[i].head=NULL;
		}
		clara(clist,input,inputbig,rows,tablescounter,cent);
		silh=silhouette(clist,inputbig,cent,rows,tablescounter);//FILE * fp)
		if(loop==0)
		{
			bestsilh=silh;
			clistfinal=malloc(cent*sizeof(struct clustlist));
			for(i=0;i<cent;i++)
			{
				clistfinal[i].key=malloc(rows*sizeof(double*));
				for(j=0;j<rows;j++)
					clistfinal[i].key[j]=malloc(3*sizeof(double*));
				clistfinal[i].head=NULL;
			}
			swapclist(clistfinal,clist,cent,rows);
			bestcent=cent;
		}
		else
		{
			if(silh>bestsilh)
			{
				bestsilh=silh;
				for(i=0;i<bestcent;i++)
				{
					for(j=0;j<rows;j++)
						free(clistfinal[i].key[j]);
					free(clistfinal[i].key);
				}
				free(clistfinal);
				clistfinal=malloc(cent*sizeof(struct clustlist));
				for(i=0;i<cent;i++)
				{
					clistfinal[i].key=malloc(rows*sizeof(double*));
					for(j=0;j<rows;j++)
						clistfinal[i].key[j]=malloc(3*sizeof(double*));
					clistfinal[i].head=NULL;
				}
				swapclist(clistfinal,clist,cent,rows);
				bestcent=cent;
			}
		}
		for(i=0;i<cent;i++)
		{
			for(j=0;j<rows;j++)
				free(clist[i].key[j]);
			free(clist[i].key);
		}
		free(clist);
	}
	printf("k : %d\n",bestcent);
	printf("s : %f\n",bestsilh);
	for(i=0;i<bestcent;i++)
	{
		counter=0;
		temp=clistfinal[i].head;
		while(temp!=NULL)
		{
			temp=temp->next;
			counter++;
		}
		matr=malloc(counter*sizeof(int));
		temp=clistfinal[i].head;
		counter=0;
		while(temp!=NULL)
		{
			matr[counter]=temp->id;
			temp=temp->next;
			counter++;
		}
		bubble_sort(matr,counter);
		if(counter==0)
		{
			printf("Empty cluster!\n");
			continue;
		}
		else
			for(j=0;j<counter;j++)
				printf("%d ",matr[j]);
		printf("\n");
		free(matr);
	}
	for(i=0;i<bestcent;i++)
	{
		for(j=0;j<rows;j++)
			free(clistfinal[i].key[j]);
		free(clistfinal[i].key);
	}
	free(clistfinal);
	return 1;
}

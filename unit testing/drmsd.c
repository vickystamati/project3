#include <stdio.h>
#include <stdlib.h>
#include "rmsd.h"
#include <math.h>
#include "drmsd.h"
#include <time.h>
#include <assert.h>
#define bufSize 2048

int drmsdmethod(double ***input,int tablescounter,int rows)
{
	int i,j,r[3],loop,size;
	double **firstm,**disttable;
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
	for(i=0;i<size;i++)
		disttable[i]=malloc(3*sizeof(double));
	assert(distancesort(firstm,disttable,size,rows));
	assert(clustering(input,disttable,tablescounter,rows,r));
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



int clustering(double ***input,double ** disttable,int tablescounter,int rows,int* num)//num=r
{
	FILE * fp;
	clock_t start,end;
	int i,j,z,rec,cent,choice,count,count2,loop,rtab[7],ktab[7],choicetab[7],point;
	struct clustlist * clist;
	double ** matrix,**centmatrix,sum,sum2,silh[7],timetab[7],tempsilh;
	for(i=0;i<7;i++)//paradoxi arxikopoiisi
		silh[i]=100;
	srand(time(NULL));

	for(rec=0;rec<2;rec++)//posa diaforetika k
	{
		point=0;
		cent=rand() % 20 + 5;
		for(loop=0;loop<3;loop++)//gia na ginoun oles oi 7 ekteleseis
		{
			for(choice=0;choice<3;choice++)//gia na ginoun oles oi 7 ekteleseis
			{
				clist=malloc(cent*sizeof(struct clustlist));
				for(i=0;i<cent;i++)
				{
					clist[i].key=malloc(rows*sizeof(double*));
					for(j=0;j<rows;j++)
						clist[i].key[j]=malloc(3*sizeof(double));
					clist[i].head=NULL;
				}
				centmatrix=malloc(cent*sizeof(double*));
				for(i=0;i<cent;i++)
					centmatrix[i]=malloc(num[loop]*sizeof(double));
				matrix=malloc(tablescounter*sizeof(double*));
				for(i=0;i<tablescounter;i++)
					matrix[i]=malloc(num[loop]*sizeof(double));
				for(i=0;i<tablescounter;i++)
				{
					count=0;
					count2=0;
					for(j=0;j<num[loop];j++)
					{
						if(choice==1)//ayksousa
						{
							if(i==0)
								matrix[i][j]=disttable[j][0];
							else
								matrix[i][j]=eucldistance(input,i,disttable[j][1],disttable[j][2]);
						}
						else if(choice==2)//fthinousa
						{
							if(i==0)
								matrix[i][j]=disttable[num[loop]-j][0];
							else
								matrix[i][j]=eucldistance(input,i,disttable[num[loop]-j][1],disttable[num[loop]-j][2]);
						}
						else //random
						{
								count2++;
								if(count>=rows)
									count=0;
								if(count2>=rows)
									count2=0;
								matrix[i][j]=eucldistance(input,i,count,count2);
								count++;
						}
					}
				}
				start=clock();
				assert(kmeansstart(clist,input,matrix,centmatrix,tablescounter,rows,cent,num[loop]));
				assert(kmeansmethod(clist, input,matrix,centmatrix,tablescounter,rows,num[loop],cent,disttable,choice));
				end=clock();
				tempsilh=drmsdsilhouette(clist,matrix,centmatrix,cent,num[loop]);
				if(tempsilh<silh[point])//ama vrei kalitero
				{
					choicetab[point]=choice;
					ktab[point]=cent;
					silh[point]=tempsilh;
					timetab[point]=(double)(end - start) / CLOCKS_PER_SEC;
					rtab[point]=num[loop];
				}
				point++;//ayksanei i thesi tou pinaka
				for(i=0;i<cent;i++)
					free(centmatrix[i]);
				free(centmatrix);
				for(i=0;i<tablescounter;i++)
					free(matrix[i]);
				free(matrix);
				for(i=0;i<cent;i++)
				{
					for(j=0;j<rows;j++)
						free(clist[i].key[j]);
					free(clist[i].key);
				}
				free(clist);
				if(loop==2) break;//an einai to r=n an 2 tote kane mono mia klisi tis sinartisis
			}
		}
	}
	fp=fopen("experim.dat","w");
	for(i=0;i<7;i++)
		fprintf(fp,"r=%d , T=%d , k=%d , s=%f , time=%f\n",rtab[i],choicetab[i],ktab[i],silh[i],timetab[i]);
	fclose(fp);
}

int meanspam(double *** input,struct clustlist * clist,int tablescounter,int rows,int cent,double ** matrix,double ** centmatrix,int num)
{
	int i,j,k,flag;
	double **pammatr;
	pammatr=malloc(cent*sizeof(double*));
	for(i=0;i<cent;i++){
		pammatr[i]=malloc(2*sizeof(double));
	}
	for(k=0;k<tablescounter;k++)
	{
		flag=1;
		for(i=0;i<cent;i++)
		{	
			pammatr[i][0]=i;
			if(k==clist[i].id)
			{
				flag=0;
				break;
			}
		}
		if(flag==1)
		{
			for(j=0;j<cent;j++)
				pammatr[j][1]=drmsd(matrix,centmatrix,k,j,num);
			bubble_sort2d(pammatr, cent);
			assert(insertassig(pammatr,clist,k,input,rows));
		}
	}
	for(i=0;i<cent;i++)
		free(pammatr[i]);
	free(pammatr);
	return 1;
}


int kmeansstart(struct clustlist * clist,double ***input,double ** matrix,double ** centmatrix,int tablescounter,int rows,int cent,int num)
{
	int i,j,z,flag,random,listmatr[cent];
	i=0;
	srand (time(NULL));
	while(i<cent)//vriskei cent tuxaia stoixeia kai ta kanei kentroidi
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
		clist[i].id=listmatr[i];
		for(j=0;j<rows;j++)
			for(z=0;z<3;z++)
				clist[i].key[j][z]=input[listmatr[i]][j][z];
	}
	for(i=0;i<cent;i++)
		for(j=0;j<num;j++)
			centmatrix[i][j]=matrix[listmatr[i]][j];
	return 1;
}

double drmsd(double ** matrix,double ** centmatrix,int k,int j,int counter)
{
	int i;
	double sum,sum2;
	sum=0;
	for(i=0;i<counter;i++)
	{
		sum2=centmatrix[j][i] - matrix[k][i];
		sum2=sum2*sum2;
		sum+=sum2;
	}
	sum=sqrt(sum);
	return sum;
}

int kmeansmethod(struct clustlist * clist,double *** input,double ** matrix,double ** centmatrix,int tablescounter,int rows,int num,int cent,double ** disttable,int choice)
{
	struct clustlist * clistfinal;
	struct centlist * temp;
	int i,j,z,k,centid[cent],end=0,idn,count,count1,count2;
	double jsum1,jsum2,**finalcentmatrix;
	finalcentmatrix=malloc(cent*sizeof(double*));
	for(i=0;i<cent;i++)
		finalcentmatrix[i]=malloc(num*sizeof(double));
	idn=tablescounter+1;
	assert(meanspam(input,clist,tablescounter,rows,cent,matrix,centmatrix,num));
	jsum1=calcj(clist,cent);
	clistfinal=malloc(cent*sizeof(struct clustlist));
	for(i=0;i<cent;i++)//desmevei ta kentroidi
	{
		clistfinal[i].key=malloc(rows*sizeof(double*));
		for(j=0;j<rows;j++)
			clistfinal[i].key[j]=malloc(3*sizeof(double));
		clistfinal[i].head=NULL;
	}
	while(end==0)//epanalipsi mexri na vrei ta kalitera
	{
		for(i=0;i<cent;i++)
			centid[i]=clist[i].id;
		for(i=0;i<cent;i++)
		{		
			clistfinal[i].id=idn;
			idn++;
			for(j=0;j<rows;j++)
			{
				for(z=0;z<3;z++)
				{
					temp=clist[i].head;
					count=0;
					while(temp!=NULL)
					{
						clistfinal[i].key[j][z]+=temp->key[j][z];//dokimazei ws neo kentroides
						temp=temp->next;
						count++;
					}
					clistfinal[i].key[j][z]/=count;
				}
			}
			for(j=0;j<cent;j++)
			{		
				if(j!=i)
				{
					clistfinal[j].id=clist[j].id;
					for(z=0;z<rows;z++)
						for(k=0;k<3;k++)
							clistfinal[j].key[z][k]=clist[j].key[z][k];
				}
			}
			for(j=0;j<cent;j++)
			{		
				count1=0;
				count2=0;
				for(z=0;z<num;z++)
				{	
					if(i==j)
					{
						if(choice==0)
					 		finalcentmatrix[j][z]=eucldistancecluster(clistfinal[j].key,disttable[z][1],disttable[z][2]);
						else if(choice==1)
					 		finalcentmatrix[j][z]=eucldistancecluster(clistfinal[j].key,disttable[num-z][1],disttable[num-z][2]);
						else
						{
							count2++;
							if(count1>=rows)
								count1=0;
							if(count2>=rows)
								count2=0;
							finalcentmatrix[j][z]=eucldistancecluster(clistfinal[j].key,count1,count2);
							count1++;
						}
					}
					else
						finalcentmatrix[j][z]=centmatrix[j][z];
				}
			}
			assert(meanspam(input,clistfinal,tablescounter,rows,cent,matrix,finalcentmatrix,num));
			jsum2=calcj(clistfinal,cent);
			if(jsum2<jsum1)
			{

				swapclist(clist,clistfinal,cent,rows);
				for(j=0;j<num;j++)
					centmatrix[i][j]=finalcentmatrix[i][j];
				jsum1=jsum2;

			}
			freeclustlist(clistfinal,cent,rows);
		}	
		for(i=0;i<cent;i++)
			if(centid[i]==clist[i].id)
				end++;
		if(end!=cent)
			end=0;
			
	}
	for(i=0;i<cent;i++)
	{
		for(j=0;j<rows;j++)
			free(clistfinal[i].key[j]);
		free(clistfinal[i].key);
	}
	free(clistfinal);
	for(i=0;i<cent;i++)
		free(finalcentmatrix[i]);
	free(finalcentmatrix);
	return 1;
}


double drmsdsilhouette(struct clustlist * clist,double **matrix,double** centmatrix,int cent,int num)//FILE * fp)
{

	int i;
	double suma,sumb,sumaver,sumall=0,counta,countb,countaver;
	struct centlist * temp;
	struct centlist *temp2;
	int flag;
	//printf("Silhuette per cluster {");
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
					suma+=drmsd(matrix,matrix,temp->id,temp2->id,num);
					counta++;
				}
				temp2=temp2->next;
			}
			countb=0;
			sumb=0;
			temp2=clist[temp->listnearid].head;
			while(temp2!=NULL)
			{
				sumb+=drmsd(matrix,matrix,temp->id,temp2->id,num);
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
		//if(i==cent-1)
		//	printf("%f}\n",sumaver);
		//else
		//	printf("%f, ",sumaver);
		sumall+=sumaver;
	}
	if(sumall<0) sumall*=-1;
	sumall/=cent;
	//printf("s = %f\n",sumall);
	return sumall;

}

double eucldistancecluster(double ** matr,int pos1,int pos2)
{
	double distance,sum,sum1,sum2;
	int j,z;
	distance=0;
	sum=0;
	sum2=0;
	for(z=0;z<3;z++)//typos euclidean
	{
		sum=matr[pos1][z] - matr[pos2][z];
		sum=sum*sum;
		sum2+=sum;
		sum=0;
	}
	distance=sqrt(sum2);//upologismos apostasis
	return distance;
}



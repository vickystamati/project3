#include <stdio.h>
#include <stdlib.h>
#include "rmsd.h"
#include "drmsd.h"
#include <time.h>
#include "clrec.h"
#include <string.h>
#include <math.h>
#define bufSize 2048



int clusrec(char* filename,char*output,int val)
{
	FILE *fp;
	int i,j,idnum,countid,prodnum,**matrix,prodpos,rating,cpro,choice;
	char buf1[bufSize],buf2[bufSize],buf3[bufSize],last[10];
	double sum,sumcount,**finalmatrix,error;
	struct user * suser;
	printf("Choose algorithm: Hamming->0  |   Cosine->1   |   Euclidean->2\n");
	scanf("%d",&choice);

	idnum=5400;
	prodnum=1000;
	suser= malloc(idnum*sizeof(struct user));
	for(i=0;i<idnum;i++)
	{
		suser[i].id=i+1;
		suser[i].rateitems=0;
	}
	matrix=malloc(idnum*sizeof(int*));
	finalmatrix=malloc(idnum*sizeof(double*));
	for(i=0;i<idnum;i++)
	{
		matrix[i]=malloc(prodnum*sizeof(int));
		finalmatrix[i]=malloc(prodnum*sizeof(double));
	}
	for(i=0;i<idnum;i++)
	{
		for(j=0;j<prodnum;j++)
		{
			matrix[i][j]=0;
			finalmatrix[i][j]=0;
		}
	}	
	countid=0;
	fp=fopen(filename,"r");
	while(fscanf(fp, "%s %s %s", buf1,buf2,buf3)!=EOF)
	{
		if((strcmp(last,buf1))!=0)
			countid++;
		suser[countid-1].rateitems++;
		strcpy(last,buf1);
		prodpos=atoi(buf2);	
		rating=atoi(buf3);
		matrix[countid-1][prodpos-1]=rating;
	}
	fclose(fp);
	memset(buf1,0,sizeof(buf1));
	memset(buf2,0,sizeof(buf2));
	memset(buf3,0,sizeof(buf3));
	fp=fopen(filename,"r");
	for(i=0;i<idnum;i++)
		suser[i].ratekey=malloc(suser[i].rateitems*sizeof(int));
	countid=0;
	while(fscanf(fp, "%s %s %s", buf1,buf2,buf3)!=EOF)
	{
		if((strcmp(last,buf1))!=0)
		{
			countid++;
			cpro=0;
		}

		strcpy(last,buf1);
		suser[countid-1].ratekey[cpro]=atoi(buf2);	
		cpro++;
	}
	fclose(fp);
	if(choice==2)
	{
		for(i=0;i<idnum;i++)
		{
			sum=0;
			sumcount=0;
			for(j=0;j<prodnum;j++)
			{
				if(matrix[i][j]!=0)//iparxei vathmologia
				{
					sum+=matrix[i][j];
					sumcount++;
				}
			}
			sum/=sumcount;
			for(j=0;j<prodnum;j++)
				if(matrix[i][j]!=0)
					finalmatrix[i][j]=matrix[i][j]-sum;
		}
	}
	else if(choice==0)
	{
		for(i=0;i<idnum;i++)
		{
			for(j=0;j<prodnum;j++)
			{
				if(matrix[i][j]<=2)//pairnei 0
					finalmatrix[i][j]=0;
				else
					finalmatrix[i][j]=1;
		
			}
		}
	}
	else
		for(i=0;i<idnum;i++)
			for(j=0;j<prodnum;j++)
				finalmatrix[i][j]=matrix[i][j];

	if(val==0)
		clusteringmethod(matrix,finalmatrix,prodnum,idnum,suser,output,val,&error,choice);
	else
	{
		error=0;
		for(i=0;i<1;i++)
			clusteringmethod(matrix,finalmatrix,prodnum,idnum,suser,output,val,&error,choice);
		error=error/10;
		fp=fopen(output,"w");
		if(choice==0)
			fprintf(fp,"Hamming NN-LSH  ");
		else if(choice==1)
			fprintf(fp,"Cosine NN-LSH  ");
		else if(choice==2)
			fprintf(fp,"Euclidean NN-LSH  ");
		fprintf(fp,"MAE: %f\n",error);
		fclose(fp);
	}

}

int clusteringmethod(int ** matrix,double ** finalmatrix,int productnum,int idnum,struct user * suser,char * output,int val,double*error,int choice)
{
	struct clustlist * clist, *clistfinal;
	struct centlist *temp;
	int knum=20,cent,loop,bestcent,*matr,counter,i,j;
	double silh,bestsilh;
	srand(time(NULL));
	for(loop=0;loop<10;loop++)
	{
		cent=rand() %knum+5;
		clist=malloc(cent*sizeof(struct clustlist));
		for(i=0;i<cent;i++)
			clist[i].head=NULL;
		claramethod(clist,matrix,finalmatrix,productnum,idnum,cent,suser);
		silh=clustsilhouette(clist,matrix,finalmatrix,cent,productnum,choice,suser);
		if(loop==0)
		{
			bestsilh=silh;
			clistfinal=malloc(cent*sizeof(struct clustlist));
			for(i=0;i<cent;i++)
				clistfinal[i].head=NULL;
			swapclistnokey(clistfinal,clist,cent);
			bestcent=cent;
		}
		else
		{
			if(silh>bestsilh)
			{
				bestsilh=silh;
				free(clistfinal);
				clistfinal=malloc(cent*sizeof(struct clustlist));
				for(i=0;i<cent;i++)
					clistfinal[i].head=NULL;
				swapclistnokey(clistfinal,clist,cent);
				bestcent=cent;
			}
		}
		free(clist);
	}
	if(val==0)
		findbest(clistfinal,finalmatrix,suser,productnum,idnum,choice,bestcent,matrix,output);
	else
		*error+=clustvalidate(clistfinal,finalmatrix,suser,productnum,idnum,choice,matrix);
	free(clistfinal);
	return 1;
}

int claramethod(struct clustlist * clist,int ** matrix,double ** finalmatrix,int productnum,int idnum,int cent,struct user * suser)
{
	struct clustlist*clistptr;
	int i,j,s=2;//apo theoria
	double jsum,jsumclara;
	clistptr=malloc(cent*sizeof(struct clustlist));
	for(i=0;i<cent;i++)
		clistptr[i].head=NULL;
	for(i=0;i<s;i++)
	{
		clarafclu(clistptr,matrix,finalmatrix,productnum,idnum,cent,suser);
		if(i==0)
		{
			jsumclara=calcj(clistptr ,cent );
			swapclistnokey(clist,clistptr,cent);
			freeclustlistnokey(clistptr,cent);
		}
		else
		{
			jsum=calcj(clistptr ,cent );
			if(jsum<jsumclara)
			{
				jsumclara=jsum;
				freeclustlistnokey(clist,cent);
				swapclistnokey(clist,clistptr,cent);
			}
			freeclustlistnokey(clistptr,cent);
		}
	}

	free(clistptr);
	return 1;
}


int freeclustlistnokey(struct clustlist * clist,int cent)
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
			free(temp);
		}
		clist[i].head=NULL;
	}
	return 1;
}

int swapclistnokey(struct clustlist * clist,struct clustlist * tempclist,int cent)
{
	int i,j,z;
	struct centlist * temp, *temp2;
	for(i=0;i<cent;i++)
	{
		clist[i].id=tempclist[i].id;
		temp=tempclist[i].head;
		while(temp!=NULL)
		{
			if(clist[i].head==NULL)
			{
				clist[i].head=malloc(sizeof(struct centlist));
				clist[i].head->id=temp->id;
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


int clarafclu(struct clustlist * clist , int ** matrix,double **finalmatrix,int productnum,int idnum,int cent,struct user * suser)
{
	int i,j,z,m,r,counter;
	int random,nnum,flag;
	int randmatr[cent];
	int *listmatr;
	double jval,jtemp;
	struct clustlist * clisttemp,*clistfinal;
	struct centlist * centtemp;
	srand(time(NULL));
	nnum=20 + 2*cent;
	listmatr=malloc(nnum*sizeof(int));
	i=0;
	while(i<nnum)
	{
		flag=0;
		random=rand()%idnum+1;
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
		random=rand()%nnum+1;
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
		i++;
	}
	medpamclu(matrix,finalmatrix,clist,productnum,idnum,cent,listmatr,nnum,suser);
	jval=calcj(clist ,cent);
	clistfinal=malloc(cent*sizeof(struct clustlist));
	clisttemp=malloc(cent*sizeof(struct clustlist));
	for(i=0;i<cent;i++)
	{
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
						clisttemp[m].id=centtemp->id;
					else
						clisttemp[m].id=clist[m].id;
				}
				medpamclu(matrix,finalmatrix,clisttemp,productnum,idnum,cent,listmatr,nnum,suser);
				jtemp=calcj(clisttemp ,cent);
				if(jtemp<jval)
				{
					freeclustlistnokey(clistfinal,cent);
					swapclistnokey(clistfinal,clisttemp,cent);
					freeclustlistnokey(clisttemp,cent);
					jval=jtemp;
				}
				else
					freeclustlistnokey(clisttemp,cent);
				centtemp=centtemp->next;
				counter++;

			}
		}
	}
	if(clistfinal[0].head!=NULL)
	{
		freeclustlistnokey(clist,cent);
		swapclistnokey(clist,clistfinal,cent);
	}
	freeclustlistnokey(clisttemp,cent);
	freeclustlistnokey(clistfinal,cent);
	freeclustlistnokey(clist,cent);
	medpamclu(matrix,finalmatrix,clist,productnum,idnum,cent,listmatr,0,suser);
	free(clistfinal);
	free(clisttemp);
	free(listmatr);
	return 1;
}

int medpamclu(int **matrix,double ** finalmatrix,struct clustlist * clist,int productnum,int idnum,int cent,int *matr,int nnum,struct user * suser)
{
	int i,j,k,flag,loop;
	double **pammatr;
	if(nnum==0)
		loop=idnum;
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
				if(k+1==clist[i].id)
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
				{	
					pammatr[j][1]=eucldistanceclu(finalmatrix,k+1,clist[j].id,productnum,suser);//euclidean
					//pammatr[j][1]=cosdistanceclu(matrix,k+1,clist[j].id,productnum,suser);//cosine
					//pammatr[j][1]=hamdistanceclu(finalmatrix,k+1,clist[j].id,productnum,suser);//hamm

				}
				else
				{	
					pammatr[j][1]=eucldistanceclu(finalmatrix,matr[k],clist[j].id,productnum,suser);
					//pammatr[j][1]=cosdistanceclu(matrix,matr[k],clist[j].id,productnum,suser);//cosine
					//pammatr[j][1]=hamdistanceclu(finalmatrix,matr[k],clist[j].id,productnum,suser);//hamm
				}
			}
			bubble_sort2d(pammatr, cent);
			if(nnum==0)
				insertassigclu(pammatr,clist,k+1,suser);	
			else
				insertassigclu(pammatr,clist,matr[k],suser);	
		}
	}
	for(i=0;i<cent;i++)
		free(pammatr[i]);
	free(pammatr);
	return 1;
}

int insertassigclu(double ** pammatr,struct clustlist * clist,int k,struct user * suser)
{
	int i,j,z,m;
	i=pammatr[0][0];
	if (clist[i].head == NULL)//an head null kane eisagwgi
	{
		clist[i].head=malloc(sizeof(struct centlist));
		clist[i].head->id=k;
		clist[i].head->distance=pammatr[0][1];
		j=pammatr[1][0];
		clist[i].head->listnearid=j;
		clist[i].head->distancenear=pammatr[1][1];
		suser[k-1].centroid=i;
		clist[i].head->next=NULL;
	}
	else
	{
		struct centlist *current = clist[i].head;
		while (current->next != NULL) 
			current = current->next;
		current->next= malloc(sizeof(struct centlist));
		current->next->id=k;
		current->next->distance=pammatr[0][1];
		j=pammatr[1][0];
		current->next->listnearid=j;
		current->next->distancenear=pammatr[1][1];
		suser[k-1].centroid=i;
		current->next->next=NULL;
	}
	return 1;
}

double eucldistanceclu(double ** matr,int pos1,int pos2,int counter,struct user * suser)
{
	double distance,sum,sum1,sum2,flag;
	int i,j,z,diff,*diffmatr;
	distance=0;
	diff=0;
	for(i=0;i<suser[pos2-1].rateitems;i++)//metraei ta diaforetika
	{
		flag=0;
		for(j=0;j<suser[pos1-1].rateitems;j++)
		{	
			if(suser[pos2-1].ratekey[i]==suser[pos1-1].ratekey[j])
			{
				flag=1;
				break;
			}
		}
		if(flag==0)
			diff++;
	}
	if(diff!=0)
	{
		diffmatr=malloc(diff*sizeof(int));
		diff=0;
		for(i=0;i<suser[pos2-1].rateitems;i++)
		{
			flag=0;
			for(j=0;j<suser[pos1-1].rateitems;j++)
			{	
				if(suser[pos2-1].ratekey[i]==suser[pos1-1].ratekey[j])
				{
					flag=1;
					break;
				}
			}
			if(flag==0)
			{
				diffmatr[diff]=suser[pos2-1].ratekey[i];
				diff++;
			}
		}
		diff=1;
	}
	sum=0;
	sum2=0;
	if(diff!=0)
	{		
		for(i=0;i<diff;i++)
		{
			sum=matr[pos1-1][diffmatr[i]-1] - matr[pos2-1][diffmatr[i]-1];
			sum=sum*sum;
			sum2+=sum;
			sum=0;
		}
		free(diffmatr);
	}
	else
	{
		for(i=0;i<suser[pos1-1].rateitems;i++)
		{
			sum=matr[pos1-1][suser[pos1-1].ratekey[i]-1] - matr[pos2-1][suser[pos1-1].ratekey[i]-1];
			sum=sum*sum;
			sum2+=sum;
			sum=0;
		}
	}

	distance=sqrt(sum2);//upologismos apostasis
	return distance;
}

double cosdistanceclu(int ** matr,int pos1,int pos2,int counter,struct user * suser)
{
	double distance,sum,sum1,sum2;
	int i,j,z,diff,*diffmatr,flag;
	distance=0;
	diff=0;
	for(i=0;i<suser[pos2-1].rateitems;i++)//metraei ta diaforetika
	{
		flag=0;
		for(j=0;j<suser[pos1-1].rateitems;j++)
		{	
			if(suser[pos2-1].ratekey[i]==suser[pos1-1].ratekey[j])
			{
				flag=1;
				break;
			}
		}
		if(flag==0)
			diff++;
	}
	if(diff!=0)
	{
		diffmatr=malloc(diff*sizeof(int));
		diff=0;
		for(i=0;i<suser[pos2-1].rateitems;i++)
		{
			flag=0;
			for(j=0;j<suser[pos1-1].rateitems;j++)
			{	
				if(suser[pos2-1].ratekey[i]==suser[pos1-1].ratekey[j])
				{
					flag=1;
					break;
				}
			}
			if(flag==0)
			{
				diffmatr[diff]=suser[pos2-1].ratekey[i];
				diff++;
			}
		}
		diff=1;
	}
	if(diff!=0)
	{		
		sum=0;
		sum1=0;
		sum2=0;
		for(i=0;i<diff;i++)
		{
			sum+=matr[pos1-1][suser[pos1-1].ratekey[i]-1] * matr[pos2-1][suser[pos1-1].ratekey[i]-1];
			sum1+=matr[pos1-1][suser[pos1-1].ratekey[i]-1] * matr[pos1-1][suser[pos1-1].ratekey[i]-1];
			sum2+=matr[pos2-1][suser[pos1-1].ratekey[i]-1] * matr[pos2-1][suser[pos1-1].ratekey[i]-1];
		}
		free(diffmatr);
	}
	else
	{
		sum=0;
		sum1=0;
		sum2=0;
		for(i=0;i<suser[pos1-1].rateitems;i++)
		{
			sum+=matr[pos1-1][suser[pos1-1].ratekey[i]-1] * matr[pos2-1][suser[pos1-1].ratekey[i]-1];
			sum1+=matr[pos1-1][suser[pos1-1].ratekey[i]-1] * matr[pos1-1][suser[pos1-1].ratekey[i]-1];
			sum2+=matr[pos2-1][suser[pos1-1].ratekey[i]-1] * matr[pos2-1][suser[pos1-1].ratekey[i]-1];
		}
	}
	if(sum1>0 && sum2>0)//epeidi einai riza, elegxos na einai diaforetiko ap to 0
	{
		distance=1-(sum/(sqrt(sum1) * sqrt(sum2)));
	}
	else 
	{
		distance=0;
	}
	return distance;
}


double hamdistanceclu(double ** matr,int pos1,int pos2,int counter,struct user * suser)
{
	int i,j;
	double distance=0;
	for(i=0;i<counter;i++)
		if(matr[pos1-1][i]!=matr[pos2-1][i])
			distance++;
	return distance;
}




double clustsilhouette(struct clustlist * clist,int **matrix,double** finalmatrix,int cent,int productnum,int choice,struct user * suser)//FILE * fp)
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
					suma+=eucldistanceclu(finalmatrix,temp->id,temp2->id,productnum,suser);
					//suma+=cosdistanceclu(matrix,temp->id,temp2->id,productnum,suser);//cosine
					//suma+=hamdistanceclu(finalmatrix,temp->id,temp2->id,productnum,suser);//hamm
					counta++;
				}
				temp2=temp2->next;
			}
			countb=0;
			sumb=0;
			temp2=clist[temp->listnearid].head;
			while(temp2!=NULL)
			{
				sumb+=eucldistanceclu(finalmatrix,temp->id,temp2->id,productnum,suser);
				//sumb+=cosdistanceclu(matrix,temp->id,temp2->id,productnum,suser);
				//sumb+=hamdistanceclu(finalmatrix,temp->id,temp2->id,productnum,suser);
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

int findbest(struct clustlist *clist,double **finalmatrix,struct user *suser,int productnum,int idnum,int choice,int cent,int **matrix,char * output)
{
	FILE * fp;
	int i,j,cnt;
	clock_t start,end;
	double **bestmatrix,dista,zsum,sumdist,**bubblematr,ftime;
	struct centlist  *temp,*temp2;
	bestmatrix=malloc(idnum*sizeof(double*));
	for(i=0;i<idnum;i++)
		bestmatrix[i]=malloc(productnum*sizeof(double));
	start=clock();
	for(i=0;i<cent;i++)
	{
		temp=clist[i].head;
		while(temp!=NULL)
		{
			for(j=0;j<productnum;j++)
			{
				zsum=0;
				sumdist=0;
				if(finalmatrix[temp->id-1][productnum]==0)
				{
					temp2=clist[i].head;
					while(temp2!=NULL)
					{
						if(temp->id!=temp2->id)
						{
							if(finalmatrix[temp2->id-1][j]!=0)
							{
								if(choice==0)
									dista=hamdistanceclu(finalmatrix,temp->id,temp2->id,productnum,suser);
								else if(choice==1)
									dista=cosdistanceclu(matrix,temp->id,temp2->id,productnum,suser);
								else//eucl
									dista=eucldistanceclu(finalmatrix,temp->id,temp2->id,productnum,suser);
								zsum+=dista;
								sumdist+=dista * finalmatrix[temp2->id-1][j];

							}
	
						}
						temp2=temp2->next;
					}
				}
				bestmatrix[temp->id-1][j]=(1/(1+zsum)) * sumdist;
			}
			temp=temp->next;
		}
	}
	end=clock();
	fp=fopen(output,"w");
	if(choice==2)
		fprintf(fp,"Euclidean Clustering\n");
	else if(choice==0)
		fprintf(fp,"Hamming Clustering\n");
	else
		fprintf(fp,"Cosine Clustering\n");
	
	for(i=0;i<idnum;i++)
	{
		cnt=0;
		for(j=0;j<productnum;j++)
			if(bestmatrix[i][j]!=0)
				cnt++;
		bubblematr=malloc(cnt*sizeof(double*));
		for(j=0;j<cnt;j++)
			bubblematr[j]=malloc(2*sizeof(double));
		cnt=0;
		for(j=0;j<productnum;j++)
		{
			if(bestmatrix[i][j]!=0)
			{
				bubblematr[cnt][0]=bestmatrix[i][j];
				bubblematr[cnt][1]=j+1;
				cnt++;
			}
		}	
		bubble_sort2darr(bubblematr, cnt);
		if(cnt>5)
			cnt=5;
		fprintf(fp,"u%d ",i+1);
		for(j=0;j<cnt;j++)
		{		

			fprintf(fp,"item%d ",(int)bubblematr[j][1]);
		}
		fprintf(fp,"\n");
	}
	ftime=(double)(end - start) / CLOCKS_PER_SEC;
	fprintf(fp,"Execution Time : %f\n",ftime);
	fclose(fp);
	return 1;
}

void bubble_sort2darr(double** matr, int counter)
{
	double t1,t2;
	int i,j;
	for (i=0;i<(counter-1);i++)
	{
		for (j=0;j<counter-i-1;j++)
		{
			if (matr[j][0] < matr[j+1][0])
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


double clustvalidate(struct clustlist*clist,double **finalmatrix,struct user *suser,int prodnum,int idnum,int choice,int**matrix)
{
	struct centlist * temp;
	int i,j,z,size,*randidtab,idr,flag,ratedid,pos;
	double zsum,sumdist,dista,**bestmatrix;
	size=idnum/10;
	randidtab=malloc(size*sizeof(int));
	bestmatrix=malloc(idnum*sizeof(double*));
	for(i=0;i<idnum;i++)
		bestmatrix[i]=malloc(prodnum*sizeof(double));
	i=0;
	srand(time(NULL));
	while(i<size)
	{
		idr=rand()%idnum+1;
		if(i!=0)
		{
			flag=0;
			for(j=0;j<i;j++)
			{
				if(randidtab[j]==idr)
					flag=1;
			}
			if(flag!=1)
			{
				randidtab[i]=idr;
				i++;
			}		
		}
		else
		{
			randidtab[i]=idr;
			i++;
		}
	}
	for(i=0;i<size;i++)
	{
		pos=randidtab[i];
		for(j=0;j<10;j++)
		{	
			ratedid=suser[pos-1].ratekey[j];			
			zsum=0;
			sumdist=0;	
			temp=clist[suser[pos-1].centroid].head;
			while(temp!=NULL)
			{		
				if(temp->id!=pos)	
				{		
					if(finalmatrix[temp->id-1][ratedid-1]!=0)
					{				
						if(choice==0)
							dista=hamdistanceclu(finalmatrix,pos,temp->id,prodnum,suser);
						else if(choice==1)
							dista=cosdistanceclu(matrix,pos,temp->id,prodnum,suser);
						else//eucl
							dista=eucldistanceclu(finalmatrix,pos,temp->id,prodnum,suser);
						zsum+=dista;
						sumdist+=dista * matrix[temp->id-1][ratedid-1];
						
					}
				}
				temp=temp->next;
			}
			
			bestmatrix[pos-1][ratedid-1]=(1/(1+zsum)) * sumdist;
			if(choice==1)
				if(bestmatrix[pos-1][ratedid-1]==0)
					bestmatrix[pos-1][ratedid-1]=(double)rand()/RAND_MAX*2.0-1.0;
		}
	}
	sumdist=0;
	for(i=0;i<size;i++)
	{
		pos=randidtab[i];
		for(j=0;j<10;j++)
		{
			ratedid=suser[pos-1].ratekey[j];
			zsum=matrix[pos-1][ratedid-1]-bestmatrix[pos-1][ratedid-1];
			if(zsum<0) zsum*=-1;
			sumdist+=zsum;
		}
	}
	sumdist/=(size * 10);

	free(randidtab);
	return sumdist;

}


/*
void readfile(char * filename,int idnum,int prodnum)
{
	FILE *fp;
	int i,j,v,countid,prodnum,temp,y,flag,kma[1000],* prod;
	char buf1[bufSize],buf2[bufSize],buf3[bufSize], last[10];
	fp=fopen(filename,"r");
	countid=0;
	prodnum=0;
	while(fscanf(fp, "%s %s %s", buf1,buf2,buf3)!=EOF)
	{
		if(strcmp(last,buf1)==0)
			
			i=0;
		else
			countid++;
		prodnum++;
		strcpy(last,buf1);

	}
	printf("prodnum %d , id  %d\n",prodnum,countid);
	
	fclose(fp);
	prod=malloc(prodnum*sizeof(int));
	fp=fopen(filename,"r");
	i=0;
	while(fscanf(fp, "%s %s %s", buf1,buf2,buf3)!=EOF)
	{
		prod[i]=atoi(buf2);
		i++;


	}
	bubble_sort1(prod, prodnum);

	v=1;
	kma[0]=1;
	for(i=0;i<prodnum;i++)
	{
		for(j=0;j<prodnum;j++)
		{
			if(prod[i]==prod[j])
			{
				flag=0;
				for(y=0;y<v;y++)
				{
					if(prod[i]==kma[y])
					{
						flag=1;
					}
				}
			if(flag==0){
					kma[v]=prod[i];
						v++;}
		
			}
		}

	}
}

*/

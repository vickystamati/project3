#include <stdio.h>
#include <stdlib.h>
#include "rmsd.h"
#include "drmsd.h"
#include <time.h>
#include "clrec.h"
#include <string.h>
#include "nnlsh.h"
#include <math.h>
#include <assert.h>
#define bufSize 2048

int nnlsh(char * input, char * output,int val)
{
	long hashsize=2;
	FILE *fp;
	char buf1[bufSize],buf2[bufSize],buf3[bufSize],last[10];
	struct hashtable *hasht;
	struct user * suser;
	int i,j,idnum,cpro,prodnum,countid,k=4,l=3,choice,**matrix,prodpos,rating;//k kai l by default
	double **finalmatrix,sum,sumcount,error;
	printf("Choose algorithm: Hamming->0  |   Cosine->1   |   Euclidean->2\n");
	scanf("%d",&choice);
	for(i=1;i<k;i++)
	{
		hashsize*=2;
	}
	hasht=malloc(l*sizeof(struct hashtable));
	for(i=0;i<l;i++)
		hasht[i].lista=malloc(hashsize*sizeof(struct list));
	for(i=0;i<l;i++)
		for(j=0;j<hashsize;j++)
			hasht[i].lista[j].head=NULL;
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
	fp=fopen(input,"r");
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
	fp=fopen(input,"r");
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
	if(choice==0)
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
	if(choice==1)
		for(i=0;i<idnum;i++)
			for(j=0;j<prodnum;j++)
				finalmatrix[i][j]=matrix[i][j];
	if(val==0)
	{
		if(choice==0)
			assert(hamminglsh(hasht,finalmatrix,idnum,prodnum,suser,output,l,k,matrix,val,&error));
		else if(choice==1)
			assert(cosinelsh(hasht,finalmatrix,idnum,prodnum,suser,output,l,k,matrix,val,&error));
		else if(choice==2)
			assert(euclideanlsh(hasht,finalmatrix,idnum,prodnum,suser,output,l,k,matrix,val,&error));
	}
	else
	{
		error=0;
		for(i=0;i<10;i++)
		{
			if(choice==0)
				assert(hamminglsh(hasht,finalmatrix,idnum,prodnum,suser,output,l,k,matrix,val,&error));
			else if(choice==1)
				assert(cosinelsh(hasht,finalmatrix,idnum,prodnum,suser,output,l,k,matrix,val,&error));
			else if(choice==2)
				assert(euclideanlsh(hasht,finalmatrix,idnum,prodnum,suser,output,l,k,matrix,val,&error));
		}
		error=error/10;
		fp=fopen(output,"w");
		if(choice==0)
			fprintf(fp,"Hamming NN-LSH\n");
		else if(choice==1)
			fprintf(fp,"Cosine NN-LSH\n");
		else if(choice==2)
			fprintf(fp,"Euclidean NN-LSH\n");
		fprintf(fp,"MAE: %f\n",error);
		fclose(fp);
	}
	return 1;

}



int hamminglsh(struct hashtable * hasht,double **matrix,int idnum,int prodnum,struct user * suser,char * output,int L,int k,int **matrix2,int val,double *error)
{
	long size=2;
	int i,j,m,y,z,gfun[L][k],binary,fbin,P,end,***nearest,found,flag,**noth,loop=5*idnum,**bestid;
	double radius,more;
	for(i=1;i<k;i++)
	{
		size*=2;
	}
	srand(time(NULL));
	//P=rand()%31+20;//diastima [20,50]
	P=30;
	bestid=malloc(idnum*sizeof(int*));
	for(i=0;i<idnum;i++)
		bestid[i]=malloc(P*sizeof(int));
	nearest=malloc(L*sizeof(int**));
	for(i=0;i<L;i++)
	{
		nearest[i]=malloc(idnum*sizeof(int*));
		for(j=0;j<idnum;j++)
			nearest[i][j]=malloc(P*sizeof(int));
	}
	//gfun=malloc(L*sizeof(int*));	
	//for(i=0;i<L;i++)//kataskevazoyme ti sinartisi g
	//	gfun[i]=malloc(k*sizeof(int));
	for(i=0;i<L;i++)
		for(j=0;j<k;j++)
			gfun[i][j]=rand()%(prodnum-5)+1;
	for(z=0;z<idnum;z++)
	{	
		for(i=0;i<L;i++)
		{
			fbin=0;
			for(j=0;j<k;j++)
			{
				if(matrix[z][gfun[i][j]]==1)
				{
					binary=1;
					for(m=0;m<j;m++)
						binary*=2;
				}
				else
					binary=0;
				fbin+=binary;
			}
			assert(insertht(&hasht[i].lista[fbin],z+1));

		}
	}
	radius=firstradham(matrix,prodnum);
	loop=0;
	for(y=0;y<loop;y++)//paradoxi 5 epanalipsewn,se periptwsi pou den to vrei,tha kanei mia ekti kai tha valei ta P prwta
	{	
		for(z=0;z<idnum;z++)
		{
			end=0;
			for(i=0;i<L;i++)
			{
				fbin=0;
				for(j=0;j<k;j++)
				{			
					if(matrix[z][gfun[i][j]]==1)
					{
						binary=1;
						for(m=0;m<j;m++)
							binary*=2;
					}
					else
						binary=0;
					fbin+=binary;
				}
				found=searchht(&hasht[i].lista[fbin],z+1,nearest[i][z],radius,P,noth,matrix,0,suser,prodnum);
				if(found==P)
					end++;
			}
		}
		if(found>P)
		{
			more=radius;//kratame tin aktina pou vrike ta parapanw gia tin extra epanalipsi(an xreiastei)
			radius=radius/2;
		}
		else
			radius=radius*2;
		if(end==idnum-1)//an vrethei P tote kanei to flag 1 
		{
			flag=1;
			break;
		}
		else//alliws kanei to flag 0 kai kanei mia extra epanalipsi gia na valei ta prwta P pou tha sinantisei
			flag=0;
	}
	radius=radius*2;
	if(flag==0)
		for(z=0;z<idnum;z++)
			for(i=0;i<L;i++)
				for(j=0;j<size;j++)
					found=searchht(&hasht[i].lista[j],z+1,nearest[i][z],radius,P,noth,matrix,0,suser,prodnum);
	for(i=0;i<idnum;i++)
		for(j=0;j<P;j++)
			if(nearest[0][i][j]>0 && nearest[0][i][j]<=idnum )
				bestid[i][j]=nearest[0][i][j];
			else
				bestid[i][j]=rand()%idnum+1;
	if(val==0)
		assert(findbestlsh(bestid,matrix, L,idnum,prodnum, P,suser, 0,matrix2,output));
	else
		*error+=validate(bestid,matrix,L,idnum,prodnum,P,suser,0,matrix2,output);


	for(i=0;i<idnum;i++)
		free(bestid[i]);
	free(bestid);
	return 1;
}


int cosinelsh(struct hashtable * hasht,double **matrix,int idnum,int prodnum,struct user * suser,char * output,int L,int k,int **matrix2,int val,double *error)
{
	long size=2;
	char token[k];
	int i,j,m,y,z,gfun[L][k],binary,fbin,P,end,***nearest,found,flag,**noth,loop=5*idnum,**bestid;
	double radius,more,hfix[k*L][prodnum],ran,y1,y2,htable[L][k][prodnum],sum;
	for(i=1;i<k;i++)
	{
		size*=2;
	}
	srand(time(NULL));
	//P=rand()%31+20;//diastima [20,50]
	P=30;
	bestid=malloc(idnum*sizeof(int*));
	for(i=0;i<idnum;i++)
		bestid[i]=malloc(P*sizeof(int));
	nearest=malloc(L*sizeof(int**));
	for(i=0;i<L;i++)
	{
		nearest[i]=malloc(idnum*sizeof(int*));
		for(j=0;j<idnum;j++)
			nearest[i][j]=malloc(P*sizeof(int));
	}
	for(i=0;i<(k*L);i++)
	{
		for(j=0;j<prodnum;j++)
		{
			ran=2;
			while(ran>=1)//ftiaxnw y1<1 kai y2<1, kai ftiaxnw to rand<1 (typos marsaglia)
			{
				y1=-1+(rand()/(RAND_MAX +1.0))*(1+1);	
				y2=-1+(rand()/(RAND_MAX +1.0))*(1+1);	
				ran=(y1*y1) + (y2*y2);	
			}
			hfix[i][j]=sqrt((-2*log10(ran))/(ran)) * y1;//ftiaxno ta h
		}
	}
	for(i=0;i<L;i++)//kataskevazoyme ti sinartisi g pou periexei ton arithmo tou h poy tha mpei se kathe g(px h1 h2 h3 h5)
		for(j=0;j<k;j++)
			gfun[i][j]=1+ (rand() /( RAND_MAX + 1.0))*(L*k-1);
	for(i=0;i<L;i++)
		for(j=0;j<k;j++)
			for(z=0;z<prodnum;z++)		
				htable[i][j][z]=hfix[gfun[i][j]][z];
	for(z=0;z<idnum;z++)
	{	
		for(i=0;i<L;i++)
		{
			fbin=0;
			for(j=0;j<k;j++)
			{
				sum=0;
				for(y=0;y<suser[z].rateitems;y++)
				{
					sum+=htable[i][j][suser[z].ratekey[y]-1] * matrix[z][suser[z].ratekey[y]-1];
				}
				if(sum>=0)
					token[j]='1';
				else 
					token[j]='0';
				sum=0;
			}
			binary=atoi(token);
			fbin=turnintodecimal(binary);
			assert(insertht(&hasht[i].lista[fbin],z+1));
		}
	}
	radius=1;
	loop=0;
	for(m=0;m<loop;m++)//paradoxi 5 epanalipsewn,se periptwsi pou den to vrei,tha kanei mia ekti kai tha valei ta P prwta
	{	
		for(z=0;z<idnum;z++)
		{
			end=0;
			for(i=0;i<L;i++)
			{
				fbin=0;
				for(j=0;j<k;j++)
				{			
					sum=0;
					for(y=0;y<suser[z].rateitems;y++)
					{
						sum+=htable[i][j][suser[z].ratekey[y]-1] * matrix[z][suser[z].ratekey[y]-1];
					}
					if(sum>=0)
						token[j]='1';
					else 
						token[j]='0';
					sum=0;
				}
				binary=atoi(token);
				fbin=turnintodecimal(binary);
				found=searchht(&hasht[i].lista[fbin],z+1,nearest[i][y],radius,P,matrix2,matrix,1,suser,prodnum);
				if(found==P)
					end++;
			}
		}
		if(found>P)
		{
			more=radius;//kratame tin aktina pou vrike ta parapanw gia tin extra epanalipsi(an xreiastei)
			radius=radius/2;
		}
		else
			radius=radius*2;
		if(end==idnum-1)//an vrethei P tote kanei to flag 1 
		{
			flag=1;
			break;
		}
		else//alliws kanei to flag 0 kai kanei mia extra epanalipsi gia na valei ta prwta P pou tha sinantisei
			flag=0;
	}
	radius=40;
	for(z=0;z<idnum;z++)
	{
		for(i=0;i<L;i++)
		{
			fbin=0;
			for(j=0;j<k;j++)
			{			
				sum=0;
				for(y=0;y<suser[z].rateitems;y++)
				{
					sum+=htable[i][j][suser[z].ratekey[y]-1] * matrix[z][suser[z].ratekey[y]-1];
				}
				if(sum>=0)
					token[j]='1';
				else 
					token[j]='0';
				sum=0;
			}
			binary=atoi(token);
			fbin=turnintodecimal(binary);
			found=searchht(&hasht[i].lista[fbin],z+1,nearest[i][y],radius,P,matrix2,matrix,1,suser,prodnum);
		}
	}
	for(i=0;i<idnum;i++)
		for(j=0;j<P;j++)
			if(nearest[0][i][j]>0 && nearest[0][i][j]<=idnum )
				bestid[i][j]=nearest[0][i][j];
			else
				bestid[i][j]=rand()%idnum+1;
	if(val==0)
		assert(findbestlsh(bestid,matrix, L,idnum,prodnum, P,suser, 1,matrix2,output));
	else
		*error+=validate(bestid,matrix, L,idnum,prodnum, P,suser, 1,matrix2,output);
	for(i=0;i<idnum;i++)
		free(bestid[i]);
	free(bestid);

	return 1;
}

int euclideanlsh(struct hashtable * hasht,double **matrix,int idnum,int prodnum,struct user * suser,char * output,int L,int k,int **matrix2,int val,double *error)
{
	int hashsize=2;
	int **gfun,rfix[k],i,j,z,m,y,pos,end,P,flag,***nearest,loop=5*idnum,**bestid;
	double **vfix,*tfix,w=4,y1,y2,ran,sum,fsum,radius,found,more;
	//P=rand()%31+20;//diastima [20,50]
	P=20;
	vfix=malloc(k*L*sizeof(double*));
	for(i=0;i<k*L;i++)
		vfix[i]=malloc(prodnum*sizeof(double));
	gfun=malloc(L*sizeof(int*));
	for(i=0;i<L;i++)
		gfun[i]=malloc(k*sizeof(int));
	tfix=malloc(k*L*sizeof(int));
	bestid=malloc(idnum*sizeof(int*));
	for(i=0;i<idnum;i++)
		bestid[i]=malloc(P*sizeof(int));
	for(i=1;i<k;i++)
	{
		hashsize*=2;
	}
	nearest=malloc(L*sizeof(int**));
	for(i=0;i<L;i++)
	{
		nearest[i]=malloc(idnum*sizeof(int*));
		for(j=0;j<idnum;j++)
			nearest[i][j]=malloc(P*sizeof(int));
	}
	for(i=0;i<L;i++)//kataskevazoyme ti sinartisi g pou periexei ton arithmo tou h poy tha mpei se kathe g(px h1 h2 h3 h5)
		for(j=0;j<k;j++)
			gfun[i][j]=(rand() /( RAND_MAX + 1.0))*(L*k);
	for(i=0;i<k*L;i++)
	{
		tfix[i]=(rand() /( RAND_MAX + 1.0))*w;	//Kataskevazoume to t toy kathe h
		for(j=0;j<prodnum;j++)
		{
			ran=2;
			while(ran>=1)
			{
				y1=-1+(rand()/(RAND_MAX +1.0))*(1+1);	
				y2=-1+(rand()/(RAND_MAX +1.0))*(1+1);	
				ran=(y1*y1) + (y2*y2);	
			}
			vfix[i][j]=sqrt((-2*log10(ran))/(ran)) * y1;//to v tou kathe h
		}
	}
	for(i=0;i<k;i++)
		rfix[i]=1+(rand() /( RAND_MAX + 1.0))*128;//kataskevi r(i) gia ton upologismo tis Φ
	for(m=0;m<idnum;m++)
	{
		for(i=0;i<L;i++)
		{	
			fsum=0;
			for(j=0;j<k;j++)
			{	
				sum=0;
				for(z=0;z<10;z++)
					sum+=matrix[m][suser[m].ratekey[z]-1] * vfix[gfun[i][j]][suser[m].ratekey[z]-1];	
				sum=(sum+tfix[gfun[i][j]])/w;
				sum=floor(sum);//ypologismos h
				fsum+=rfix[j]*sum;//athroisma twn h me ta katallila r(i) gia tin kataskevi tis Φ

			}
			if(fsum<0)
				fsum*=-1;
			pos=(int)fsum%hashsize;
			assert(insertht(&hasht[i].lista[pos],m+1));
		}
	}
	radius=eucldistanceclu(matrix,10,1000,prodnum,suser);
	if(radius==0)
		radius+=5;
	//loop=0;
	for(y=0;y<loop;y++)//paradoxi 5 epanalipsewn,se periptwsi pou den to vrei,tha kanei mia ekti kai tha valei ta P prwta
	{	
		for(m=0;m<idnum;m++)
		{
			end=0;
			for(i=0;i<L;i++)
			{
				for(j=0;j<k;j++)
				{			
					sum=0;
					for(z=0;z<suser[m].rateitems;z++)
						sum+=matrix[m][suser[m].ratekey[z]-1] * vfix[gfun[i][j]][suser[m].ratekey[z]-1];
			
					sum=(sum+tfix[gfun[i][j]])/w;
					sum=floor(sum);//ypologismos h
					fsum+=rfix[j]*sum;//athroisma twn h me ta katallila r(i) gia tin kataskevi tis Φ
				}
				if(fsum<0)
					fsum*=-1;
				pos=(int)fsum%hashsize;
				found=searchht(&hasht[i].lista[pos],m+1,nearest[i][z],radius,P,matrix2,matrix,2,suser,prodnum);
				if(found==P)
					end++;
			}
		}
		if(found>P)
		{
			more=radius;//kratame tin aktina pou vrike ta parapanw gia tin extra epanalipsi(an xreiastei)
			radius=radius/2;
		}
		else
			radius=radius*2;
		if(end==idnum-1)//an vrethei P tote kanei to flag 1 
		{
			flag=1;
			break;
		}
		else//alliws kanei to flag 0 kai kanei mia extra epanalipsi gia na valei ta prwta P pou tha sinantisei
			flag=0;
	}
	flag=0;
	if(flag==0)
	{
		radius=10;
		for(m=0;m<idnum;m++)
		{
			end=0;
			for(i=0;i<L;i++)
			{
				for(j=0;j<k;j++)
				{			
					sum=0;
					for(z=0;z<suser[m].rateitems;z++)
						sum+=matrix[m][suser[m].ratekey[z]-1] * vfix[gfun[i][j]][suser[m].ratekey[z]-1];
					sum=(sum+tfix[gfun[i][j]])/w;
					sum=floor(sum);//ypologismos h
					fsum+=rfix[j]*sum;//athroisma twn h me ta katallila r(i) gia tin kataskevi tis Φ
				}
				if(fsum<0)
					fsum*=-1;
				pos=(int)fsum%hashsize;
				found=searchht(&hasht[i].lista[pos],m+1,nearest[i][z],radius,P,matrix2,matrix,2,suser,prodnum);
			}
		}
	}
	for(i=0;i<idnum;i++)
		for(j=0;j<P;j++)
			if(nearest[0][i][j]>0 && nearest[0][i][j]<=idnum )
				bestid[i][j]=nearest[0][i][j];
			else
				bestid[i][j]=rand()%idnum+1;
	if(val==0)
		assert(findbestlsh(bestid,matrix, L,idnum,prodnum, P,suser, 2,matrix2,output));
	else
		*error+=validate(bestid,matrix, L,idnum,prodnum, P,suser, 2,matrix2,output);
	//for(i=0;i<idnum;i++)
	//	free(bestid[i]);
	//free(bestid);

}

int insertht(struct list * lista,int id)
{
	struct user * temp;

	if (lista->head == NULL)
	{
		lista->head=malloc(sizeof(struct user));
		lista->head->id=id;
		lista->head->next=NULL;	
	}
	else
	{
		temp=lista->head;
		while (temp->next != NULL) 
		{
			temp = temp->next;
		}
		temp->next=malloc(sizeof(struct user));
		temp->next->id=id;
		temp->next->next=NULL;
	}	
	return 1;
}

double firstradham(double **matrix,int prodnum)
{
	int i;
	double rad=0;
	for(i=0;i<prodnum;i++)
		if(matrix[1][i]!=matrix[5][i])//random simeia gia tin prwti aktina
			rad++;
	return rad;

}



int searchht(struct list *lista,int id,int * nearest,double radius,int P,int ** matrix,double ** finalmatrix,int choice,struct user*suser,int prodnum)
{
	int i,getcount,count=0;
	double dist;
	int fin=P;
	struct user * temp;
	temp=lista->head;
	getcount=0;
	while(temp!=NULL)
	{
		if(temp->id!=id)
		{		
			if(choice==0)
				dist=hamdistanceclu(finalmatrix,id,temp->id,prodnum,suser);
			else if(choice==1)
				dist=cosdistanceclu(matrix,id,temp->id,prodnum,suser);
			else
				dist=eucldistanceclu(finalmatrix,id,temp->id,prodnum,suser);
			count = getcount-1;
			if((dist<=radius) && (count<P))
			{		
				nearest[getcount]=temp->id;
				getcount++;
			}
		}
		count++;
		temp=temp->next;
		if(choice==0 && count==fin)
			break;
	}
	return getcount;
}


int findbestlsh(int **bestid,double **finalmatrix,int L,int idnum,int prodnum,int P,struct user *suser,int choice,int ** matrix,char * output)
{
	FILE * fp;
	int i,j,z,m,f,y,flag,cnt;
	clock_t start,end;
	double **bestmatrix,dista,zsum,sumdist,**bubblematr,ftime;
	bestmatrix=malloc(idnum*sizeof(double*));
	for(i=0;i<idnum;i++)
		bestmatrix[i]=malloc(prodnum*sizeof(double));
	/*for(j=0;j<idnum;j++)
	{
		for(i=0;i<L;i++)
		{
			f=0;
			for(m=i;m<L;m++)
			{
				if(m==i)
					continue;
				for(z=0;z<P;z++)
				{
					flag=0;
					for(y=0;y<P;y++)
					{
					
						if(nearest[i][j][z]==nearest[m][j][y])
						{
							flag=1;
							break;
						}	
					}
					if(flag==0 && f<P)
					{
						best[j][f]=nearest[i][j][z];
						f++;
					}
				}
			}
		}

	}*/
	start=clock();
	for(i=0;i<idnum;i++)
	{
		for(j=0;j<prodnum;j++)
		{				
			zsum=0;
			sumdist=0;
			if(finalmatrix[i][j]==0)
			{				
				for(z=0;z<P;z++)
				{				
					if(finalmatrix[bestid[i][z]-1][j]!=0)
					{				
						if(choice==0)
							dista=hamdistanceclu(finalmatrix,i+1,bestid[i][z],prodnum,suser);
						else if(choice==1)
							dista=cosdistanceclu(matrix,i+1,bestid[i][z],prodnum,suser);
						else//eucl
							dista=eucldistanceclu(finalmatrix,i+1,bestid[i][z],prodnum,suser);
						zsum+=dista;
						sumdist+=dista * matrix[bestid[i][z]-1][j];
					}
				}
			}
			bestmatrix[i][j]=(1/(1+zsum)) * sumdist;
			if(choice==1)
				if(bestmatrix[i][j]==0)
					bestmatrix[i][j]=(double)rand()/RAND_MAX*2.0-1.0;
		}
	}
	end=clock();
	fp=fopen(output,"w");
	if(choice==0)
		fprintf(fp,"Hamming NN-LSH\n");
	else if(choice==1)
		fprintf(fp,"Cosine NN-LSH\n");
	else
	fprintf(fp,"Euclidean NN-LSH\n");
	for(i=0;i<idnum;i++)
	{
		cnt=0;
		for(j=0;j<prodnum;j++)
			if(bestmatrix[i][j]!=0)
				cnt++;
		bubblematr=malloc(cnt*sizeof(double*));
		for(j=0;j<cnt;j++)
			bubblematr[j]=malloc(2*sizeof(double));
		cnt=0;
		for(j=0;j<prodnum;j++)
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
double validate(int **bestid,double **finalmatrix,int L,int idnum,int prodnum,int P,struct user *suser,int choice,int ** matrix,char * output)
{
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
			for(z=0;z<P;z++)
			{				
				if(finalmatrix[bestid[pos-1][z]-1][ratedid-1]!=0)
				{				
					if(choice==0)
						dista=hamdistanceclu(finalmatrix,pos,bestid[pos-1][z],prodnum,suser);
					else if(choice==1)
						dista=cosdistanceclu(matrix,pos,bestid[pos-1][z],prodnum,suser);
					else//eucl
						dista=eucldistanceclu(finalmatrix,pos,bestid[pos-1][z],prodnum,suser);
					zsum+=dista;
					sumdist+=dista * matrix[bestid[i][z]-1][ratedid-1];
				}
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






int turnintodecimal(int bn)
{
	int binaryNumber,decimalNumber=0,j=1,remainder;
	binaryNumber=bn;
	while(binaryNumber!=0)
	{
		remainder=binaryNumber%10;
		decimalNumber=decimalNumber+remainder*j;
		j=j*2;
		binaryNumber=binaryNumber/10;
    }
	return decimalNumber;
}

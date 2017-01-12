#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include "list.h"
#include "rmsd.h"

/*void inserthamming(struct list * lista,long long binary,long itemid)
{
	if (lista->head == NULL)
	{
		lista->head=malloc(sizeof(struct node));
		lista->head->key=binary;
		lista->head->next=NULL;	
	}
	else
	{
		struct node *current = lista->head;
		while (current->next != NULL) 
		{
			current = current->next;
		}
		current->next=malloc(sizeof(struct node));
		current->next->key=binary;
		current->next->next=NULL;	
	}	
}


void insertcosine(struct list *lista,double * vector,long itemid,int counter)
{
	int i;
	struct node * new;
	new = malloc(sizeof(struct node));
	if(new == NULL)
	{
		printf("Failed to allocate memory!\n");
	}
	new -> key1 =malloc(counter*sizeof(double));
	for(i=0;i<counter;i++)
	{	
		new->key1[i]=vector[i];//vazei sto key1 ta stoixeia tou vector
	}
	new -> next = NULL;
	if (lista->head == NULL)//an to head einai null kanei eisagwgi
	{
		lista->head=new;
		lista->head->next=NULL;	
	}
	else
	{
		struct node *current = lista->head;
		while (current->next != NULL) //allwis vriskei ton adeio komvo
		{
			current = current->next;
		}
		current->next = new;
	}	
}


void inserteuclidean(struct list * lista,double *vector,int place)
{
	int i;
	struct node * new;
	new = malloc(sizeof(struct node));
	if(new == NULL)
	{
		printf("Failed to allocate memory!\n");
	}
	new -> key1 =malloc(3*sizeof(double));
	for(i=0;i<3;i++)
	{	
		new->key1[i]=vector[i];
	}
	new -> next = NULL;
	if (lista[place].head == NULL)
	{
		printf("mptika head sto %d\n",place);
		sleep(2);
		lista[place].head=new;
		lista[place].head->next=NULL;	
	}
	else
	{
				printf("mptika else sto %d\n",place);
		struct node *current = lista->head;
		while (current->next != NULL) 
		{
			current = current->next;
		}
		current->next = new;
	}	
}
*/

void inserteuclidean(struct list * lista,double *vector,int placelist,int placenode)
{
	int i;
	for(i=0;i<3;i++)
	{	
		lista[placelist].data[placenode].key1[i]=vector[i];	
	}

}




/*
double hamdistance(struct node *temp,struct node *temp2)//sinartisi ipologismou hamming
{
	char token[65],token2[65];
	double distance;
	int z;
	turnintobinary(temp->key ,64 ,token);
	distance=0;
	turnintobinary(temp2->key ,64 ,token2);
	for(z=0;z<64;z++)
	{
		if(token[z]!=token2[z])
			distance++;
	}
	return distance;
}


double cosdistance(struct node *temp,struct node *temp2)
{
	double distance,sum,sum1,sum2;
	int z;
	distance=0;
	for(z=0;z<3;z++)
	{
		sum=0,sum1=0,sum2=0;
		for(z=0;z<3;z++)//upologismos typou cosine
		{
			sum+=temp->key1[z] * temp2->key1[z];
			sum1+=temp->key1[z] * temp->key1[z];
			sum2+=temp2->key1[z] * temp2->key1[z];
		}
		if(sum1>0 && sum2>0)//epeidi einai riza, elegxos na einai diaforetiko ap to 0
		{
			distance=1-(sum/(sqrt(sum1) * sqrt(sum2)));
		}
		else 
		{
			distance=0;
		}
	}
	return distance;
}


double eucldistance(struct node *temp,struct node *temp2)
{
	double distance,sum,sum2;
	int z;
	distance=0;
	sum=0;
	sum2=0;
	for(z=0;z<3;z++)//typos euclidean
	{
		sum=temp->key1[z] - temp2->key1[z];
		sum=sum*sum;
		sum2+=sum;
		sum=0;
	}
	distance=sqrt(sum2);//upologismos apostasis
	return distance;
}



void distancefind(struct list * lista,int tablescounter,int rows)
{
	struct node * temp,*temp2;
	int i,j,pos;
	i=0;
	while(i<1)
	{
		pos=0;
		temp=lista[i].head;
		temp2=lista[i].head->next;
		lista[i].distance[pos]= eucldistance(temp,temp2);
		//printf("apost  %f\n",lista[i].distance[pos]);
		pos++;
		while(pos<rows-1)
		{
			temp=temp2->next;
			temp2=temp->next;
			lista[i].distance[pos]= eucldistance(temp,temp2);
			//		printf("apost  %f\n",lista[i].distance[pos]);
			pos++;
		}
		lista[i].distance[pos]= eucldistance(temp2,lista[i].head);
		//printf("apost  %f %d\n",lista[i].distance[pos],pos);
		i++;
	}
	
}



long long turnintodecimal(long long bn)
{
	long long binaryNumber,decimalNumber=0,j=1,remainder;
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
void turnintobinary(long long dn,int length,char* token)
{
	int i,k,count=0;
	char binarray[length+1];
	for (i = length-1; i >= 0; i--)
	{
		k = dn >> i;
 
		if (k & 1)
			binarray[count]='1';
		else
			binarray[count]='0';
		count++;
	}
	binarray[length]='\0';
	strcpy(token,binarray);
}

*/

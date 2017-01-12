struct node{
	long long  key;
	double * key1;
};

struct list{
	double * distance;
	struct node * data;
};

void inserthamming(struct list * ,long long,long);
void insertcosine(struct list *,double * , long,int);
void inserteuclidean(struct list *,double * ,int,int);
double hamdistance(struct node * ,struct node *);
double cosdistance(struct node *,struct node *);
double eucldistance(struct node *,struct node *);
void distancefind(struct list * ,int ,int);
void turnintobinary(long long dn,int length,char* token);
long long turnintodecimal(long long bn);

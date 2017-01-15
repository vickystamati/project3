struct list
{
	struct user * head;
};

struct hashtable
{
	struct list * lista;
};


int hamminglsh(struct hashtable * hasht,double **finalmatrix,int idnum,int prodnum,struct user* suser,char* output,int,int,int **,int,double*);
int cosinelsh(struct hashtable * hasht,double **matrix,int idnum,int prodnum,struct user * suser,char * output,int L,int k,int **,int,double*);
int nnlsh(char * input, char * output,int val);
int insertht(struct list * lista,int id);
double firstradham(double **matrix,int prodnum);
int searchht(struct list *lista,int id,int * nearest,double radius,int P,int ** matrix,double ** finalmatrix,int choice,struct user*suser,int);
int findbestlsh(int **nearest,double **,int L,int idnum,int prodnum,int P,struct user *suser,int choice,int ** matrix,char * output);	
int turnintodecimal(int bn);
int euclideanlsh(struct hashtable * hasht,double **matrix,int idnum,int prodnum,struct user * suser,char * output,int L,int k,int **matrix2,int,double*);
double validate(int **nearest,double **,int L,int idnum,int prodnum,int P,struct user *suser,int choice,int ** matrix,char * output);



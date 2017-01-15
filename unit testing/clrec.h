struct user
{
	int id;
	int *ratekey;
	int rateitems;
	int centroid;
	struct user * next;
};



int clusrec(char *,char *,int);
//void bubble_sort1(int* matr, int counter);
//void readfile(char * filename);
int clusteringmethod(int ** matrix,double ** finalmatrix,int productnum,int idnum,struct user * suser,char *,int val,double*error,int);
int freeclustlistnokey(struct clustlist * clist,int cent);
int swapclistnokey(struct clustlist * clist,struct clustlist * tempclist,int cent);
int claramethod(struct clustlist * clist,int ** matrix,double ** finalmatrix,int productnum,int idnum,int cent,struct user * suser);
int medpamclu(int **matrix,double ** finalmatrix,struct clustlist * clist,int productnum,int idnum,int cent,int *matr,int nnum,struct user * suser);
int clarafclu(struct clustlist * clist , int ** matrix,double **finalmatrix,int productnum,int idnum,int cent,struct user * suser);
int insertassigclu(double ** pammatr,struct clustlist * clist,int k,struct user * suser);
double eucldistanceclu(double ** matr,int pos1,int pos2,int counter,struct user * suser);
double clustsilhouette(struct clustlist * clist,int **matrix,double** finalmatrix,int cent,int productnum,int choice,struct user * suser);
int findbest(struct clustlist *clist,double **finalmatrix,struct user *suser,int productnum,int idnum,int choice,int cent,int **,char*);
void bubble_sort2darr(double** matr, int counter);
double cosdistanceclu(int ** matr,int pos1,int pos2,int counter,struct user * suser);
double hamdistanceclu(double ** matr,int pos1,int pos2,int counter,struct user * suser);
double clustvalidate(struct clustlist*clistfinal,double **finalmatrix,struct user *suser,int productnum,int idnum,int choice,int  **matrix);

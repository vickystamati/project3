int drmsdmethod(double ***input,int tablescounter,int rows);
int findfirstmatr(double *** input,double ** firstm,int tablescounter,int rows);
double eucldistance(double *** input,int d,int pos1,int pos2);
void bubble_sort2d_3cell(double** matr, int counter);
int distancesort(double **firstm,double**disttable,int size,int rows);
int clustering(double ***input,double ** distttable,int tablescounter,int rows,int* num);//num=r
int kmeansstart(struct clustlist * clist,double ***input,double ** matrix,double ** centmatrix,int tablescounter,int rows,int cent,int num);
double drmsd(double ** matrix,double ** centmatrix,int k,int j,int counter);
int meanspam(double *** input,struct clustlist * clist,int tablescounter,int rows,int cent,double ** matrix,double ** centmatrix,int num);
int kmeansmethod(struct clustlist * clist,double *** input,double ** matrix,double ** centmatrix,int tablescounter,int rows,int num,int cent,double ** disttable,int choice);
double drmsdsilhouette(struct clustlist * clist,double **matrix,double** centmatrix,int cent,int num);
double eucldistancecluster(double ** matr,int pos1,int pos2);

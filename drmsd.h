int drmsdmethod(double ***input,int tablescounter,int rows);
int findfirstmatr(double *** input,double ** firstm,int tablescounter,int rows);
double eucldistance(double *** input,int d,int pos1,int pos2);
void bubble_sort2d_3cell(double** matr, int counter);
int distancesort(double **firstm,double**disttable,int size,int rows);
int clustering(double ***input,double ** distttable,int tablescounter,int rows,int num,int *ktab,double *silh,double *time,int cent);//num=r
int kmeansstart(struct clustlist * clist,double ***input,double ** matrix,int tablescounter,int rows,int cent);

struct centlist
{
	int id;
	double **key;
	double distance;
	int listnearid;
	double distancenear;
	//double sav;
	struct centlist * next;
};

struct clustlist
{
	long id;
	double **key;
	struct centlist * head;
};

double crmsd(int tablescounter,int rows,double ***inputbig,int pos1,int pos2);
int readfile(FILE *fp,int tablescounter,int rows,double ***input,double***);
int clara(struct clustlist * clist,double *** input,double ***inputbig,int rows,int tablescounter,int cent);
int claraf(struct clustlist * clist , double *** input,double ***inputbig,int rows,int tablescounter,int cent);
double calcj(struct clustlist * clist,int cent);
int freeclustlist(struct clustlist * clist,int cent,int rows);
int swapclist(struct clustlist * clist,struct clustlist * tempclist,int cent,int rows);
int medpam(double *** input,double *** inputbig,struct clustlist * clist,int tablescounter,int rows,int cent,int*,int);
void bubble_sort2d(double** matr, int cent);
int insertassig(double ** pammatr,struct clustlist * clist,int k,double ***input,int rows);
double silhouette(struct clustlist * clist,double ***inputbig,int cent,int rows,int tablescounter);
void bubble_sort(int* matr, int counter);
int crmsdmethod(double *** input,double *** inputbig,int tablescounter,int rows);

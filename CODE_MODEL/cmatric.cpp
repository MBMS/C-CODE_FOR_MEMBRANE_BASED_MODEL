//This is the process for allocating memory by Zhiming Gao (gaoz@ornl.gov)

#include <iostream>
#include "cmatrix.h"

#define NR_END 0
#define FREE_ARG char*

using namespace std;

void nrerror(char error_text[])
{
	cout<<"Numerical Recipes run-time error..."<<endl;
	cout<<error_text<<endl;
	cout<<"...now exiting to system..."<<endl;
	exit(1);
}

double * *dmatrix(int nrl, int nrh, int ncl, int nch)
{	
	int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	double * *matrixs;

	/* allocate pointers to rows */
	matrixs=(double * *)malloc((nrow+NR_END)*sizeof(double *));
	if(!matrixs) nrerror("allocation failure 1 in matrix()");
	matrixs +=NR_END;
	matrixs -=nrl;

	/* allocate rows and set pointers to them */
	matrixs[nrl]=(double *)malloc((nrow*ncol+NR_END)*sizeof(double));
	if(!matrixs[nrl]) nrerror("allocation failure 2 in matrix()");
	matrixs[nrl] +=NR_END;
	matrixs[nrl] -=ncl;

	for(i=nrl+1;i<=nrh;i++) matrixs[i]=matrixs[i-1]+ncol;
	/* return pointer to array of pointers to rows */
	return matrixs;
}

void free_dmatrix(double **matrixs, int nrl, int nrh, int ncl, int nch)
{
	free((FREE_ARG) (matrixs[nrl]+ncl-NR_END));
	free((FREE_ARG) (matrixs+nrl-NR_END));
}

double max(double a, double b)
{
	if (a > b)
	{return a;}
	else
	{return b;}
}//end of max function

double min(double a, double b)
{
	if (a > b)
	{return b;}
	else
	{return a;}
}//end of min function
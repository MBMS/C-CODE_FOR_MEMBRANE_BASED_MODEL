#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

void nrerror(char error_text[]);

double **dmatrix(int nrl, int nrh, int ncl, int nch);

void free_dmatrix(double **matrixs, int nrl, int nrh, int ncl, int nch);

double max(double a, double b);

double min(double a, double b);

#endif

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#define inf 1e19;
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

double const_UBCases(double a, double b, char c, int ra, int rb, int w, int gap)
{
	double v = 0;
	if (ra == 0 && rb == 0 && gap <= w)
		v = (a - b)*(a - b);

	else if (ra == 1 && rb == 1)
	{
		v = 0;
	}
	else {
		if (c == 'd') {
			if (ra == 0 && rb == 1) {
				v = b*a*a;
			}
			else if (ra == 1 && rb == 0) {
				v = a*b*b;

			}
			else {
				v = inf;
			}
		}
		else if (c == 'l') {
			if (ra == 0 && rb == 1 && gap <= w) {
				v = b*a*a;
			}
			else if (ra == 1 && rb == 0) {
				v = b*b;
			}
			else
				v = inf;

		}
		else if (c == 't') {
			if (ra == 0 && rb == 1) {
				v = a*a;
			}
			else if (ra == 1 && rb == 0 && gap <= w) {
				v = a*b*b;
			}
			else
				v = inf;
		}
	}

	return v;
}


double const_AWarp(std::vector<double> s, std::vector<double> t, int* rs, int* rt, int w)
{
	int ns = s.size();
	int nt = t.size();

	double *x = new double[ns + 1];
	double *y = new double[nt + 1];

	int *tx = new int[ns + 1];
	int *ty = new int[nt + 1];

	double **D = new double*[ns + 1];
	for (int i = 0; i <= ns; i++)
		D[i] = new double[nt + 1];

	if (w < abs(ns - nt)) {
		printf("w need to be greater than %d.", abs(ns - nt));
		return(-1);
	}


	//transfer s to x and append one to them
	for (int i = 0; i<ns; i++)
		x[i] = s[i];
	
	x[ns] = 1;

	//transfer t to y and append one to them
	for (int i = 0; i<nt; i++)
		y[i] = t[i];
	
	y[nt] = 1;


	//calculate the timestamp of the event in s
	int iit = 0;
	for (int i = 0; i<ns + 1; i++) {
		if (x[i] > 0) {
			iit = iit + 1;
			tx[i] = iit;
		}
		else {
			iit = iit + abs(x[i]);
			tx[i] = iit;
		}
	}

	iit = 0;
	for (int i = 0; i<nt + 1; i++) {
		if (y[i] > 0) {
			iit = iit + 1;
			ty[i] = iit;
		}
		else {
			iit = iit + abs(y[i]);
			ty[i] = iit;
		}
	}

	for (int i = 0; i<ns + 1; i++)
		D[i][0] = inf;


	for (int i = 0; i<nt + 1; i++)
		D[0][i] = inf;

	D[0][0] = 0;

	int gap = 0;
	for (int i = 0; i<ns; i++)	{
		for (int j = 0; j<nt; j++)	{
			gap = abs(tx[i] - ty[j]);
			if (gap > w && ((j > 0 && ty[j - 1] - tx[i] > w) || (i > 0 && tx[i - 1] - ty[j] > w))) {
				D[i + 1][j + 1] = inf;
			}
			else {
				double a1 = D[i][j] + (s[i] - t[j])*(s[i] - t[j]);

				if (i > 0 && j > 0)
					a1 = D[i][j] + const_UBCases(s[i], t[j], 'd', rs[i], rt[j], w, gap);
				double a2 = D[i + 1][j] + const_UBCases(s[i], t[j], 'l', rs[i], rt[j], w, gap);
				double a3 = D[i][j + 1] + const_UBCases(s[i], t[j], 't', rs[i], rt[j], w, gap);

				D[i + 1][j + 1] = MIN(a3, MIN(a1, a2));
			}
		}

	}
	//print the D matrix  
	/*for( i = 0 ; i <= ns ; i++ ){
	for (j = 0 ; j <= nt ; j++){


	if (D[i][j] >= 10000000000000000000.000000){
	printf("inf ");
	}else{
	printf("%f ", D[i][j] );
	}
	}
	printf("\n");
	}*/
	double d = sqrt(D[ns][nt]);

	for (int i = 0; i <= ns; i++) 
		delete[] D[i];

	delete[] D;
	delete[] x;
	delete[] y;
	delete[] tx;
	delete[] ty;

	return d;//(double)matchLength;
}

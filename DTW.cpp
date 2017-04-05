#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>
#define inf 1e19;
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

double DTW(std::vector<double> s, std::vector<double> t)
{
	int ns = s.size();
	int nt = t.size();

	double **D = new double*[ns + 1];
	for (int i = 0; i <= ns; i++)
		D[i] = new double[nt + 1];


	for (int i = 0; i<ns + 1; i++)
		D[i][0] = inf;

	for (int i = 0; i<nt + 1; i++)
		D[0][i] = inf;

	D[0][0] = 0;


	for (int i = 0; i<ns; i++)
	{

		for (int j = 0; j<nt; j++)
		{
			double cost = (s[i] - t[j])*(s[i] - t[j]);
			D[i + 1][j + 1] = cost + MIN(D[i][j], MIN(D[i + 1][j], D[i][j + 1]));

		}
	}

	double d = sqrt(D[ns][nt]);

	for (int i = 0; i <= ns; i++)
		delete[] D[i];
	delete[] D;

	return d;//(double)matchLength;
}
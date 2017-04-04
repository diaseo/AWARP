#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include "Stack.h"
#define inf 1e19;

double UBCases(double a, double b, char c, bool ra, bool rb);
double LBCases(double a, double b, char c, bool ra, bool rb);
double Alter_path(double a, double b, bool ra, bool rb) {
	if (!ra && rb)
		return (b-1)*a*a;
	
	else if(ra && !rb)
		return (a-1)*b*b;
	
	else if(ra && rb)
		return 0;

	else return inf;
}
double Ncost(double a, double b, bool ra, bool rb) {
	if (!ra && !rb)
		return (a - b)*(a - b);

	else if (ra && !rb)
		return b*b;

	else if (!ra && rb)
		return a*a;
	else return 0.0;
}

double MIN(double a1, double a2, double a3) {
	double temp = a1;

	if (a2 < temp) temp = a2;
	if (a3 < temp) temp = a3;

	return temp;
}


double Mine_UB(std::vector<double> s, std::vector<double> t, bool* run_s, bool* run_t)
{
	int ns = s.size();
	int nt = t.size();

	double **D = new double*[ns + 1];
	char **align = new char*[ns];
	for (int i = 0; i <= ns; i++) {
		D[i] = new double[nt + 1];	
		D[i][0] = inf;
	}

	for (int i = 0; i < nt + 1; i++) 
		D[0][i] = inf;
	
	for (int i = 0; i < ns; i++) {
		align[i] = new char[nt];
	}

	D[0][0] = 0;

	for (int i = 0; i<ns; i++)
	{

		for (int j = 0; j<nt; j++)
		{

			double a1 = D[i][j] + (s[i] - t[j])*(s[i] - t[j]);

			if (i > 0 && j > 0) {
				a1 = D[i][j] + UBCases(s[i], t[j], 'd', run_s[i], run_t[j]);
				if (!run_s[i] && run_t[j]) {
					double a4 = D[i][j] + Alter_path(s[i-1], t[j], run_s[i-1], run_t[j])
					   	+ (s[i] - t[j])*(s[i] - t[j]);
					if (a4 < a1)
						a1 = a4;
				}
				
				else if(run_s[i] && !run_t[j]) {
					double a4 = D[i][j] + Alter_path(s[i], t[j-1], run_s[i], run_t[j-1])
						+ (s[i] - t[j])*(s[i] - t[j]);
					if (a4 < a1)
						a1 = a4;
				}
			}
			double a2 = D[i + 1][j] + UBCases(s[i], t[j], 't', run_s[i], run_t[j]);
			double a3 = D[i][j + 1] + UBCases(s[i], t[j], 'l', run_s[i], run_t[j]);

			double cost = MIN(a1, a2, a3);
			D[i + 1][j + 1] = cost;

			if (cost == a1) align[i][j] = 'd';
			else if (cost == a2) align[i][j] = 't';
			else align[i][j] = 'l';

		}
	}
	align[0][0] = 's';

	for (int i = 0; i <= ns; i++)
		delete[] D[i];
	delete[] D;

	// Track alignment and recompute cost
	double cost = 0;

	int s_ind = ns-1, t_ind = nt-1;
	char prev_d = 'e';
	double min = inf;
	int run_count = 0;

	// top : s_ind is fixed
	// left: t_ind is fixed
	
	while(true) {
		char dir = align[s_ind][t_ind];

		// add cost
		// s[s_ind] is non-zero and moving on the matrix while s_ind is fixed
		if (dir == 't' && !run_s[s_ind] && run_t[t_ind])
			cost += t[t_ind]*s[s_ind]*s[s_ind];
		
		// t[t_ind] is non-zero and moving on the matrix while t_ind is fixed
		else if (dir == 'l' && run_s[s_ind] && !run_t[t_ind])
			cost += s[s_ind]*t[t_ind]*t[t_ind];

		else cost += Ncost(s[s_ind], t[t_ind], run_s[s_ind], run_t[t_ind]);
	
		// moving to diagonal direction
		if (dir == 'd') {
			// previous left direction
			if (prev_d == 'l') {
				if (run_s[s_ind]) min = 0;
				else if (std::abs(s[s_ind]) < min) min = std::abs(s[s_ind]);
				run_count++;
				
				if(run_t[t_ind] && (double)run_count < t[t_ind])
					cost += (t[t_ind] - (double)run_count) * min * min;
			}
			// previous top direction
			else if(prev_d == 't') {
				if (run_t[t_ind]) min = 0;
				else if (std::abs(t[t_ind]) < min) min = std::abs(s[s_ind]);
				run_count++;

				if(run_s[s_ind] && (double)run_count < t[t_ind])
					cost += (s[s_ind] - (double)run_count) * min * min;
			}

			s_ind--; t_ind--;	
		}
		// moving to top direction
		else if (dir == 't') {
			// still in the same direction
			if(prev_d == 't') {
				run_count++;
				if(run_t[t_ind]) min = 0;
				else if(std::abs(t[t_ind]) < min) min = std::abs(t[t_ind]);
			}

			// direction changed to top
			else {
				// previous left direction
				if (prev_d == 'l') {
					if (run_s[s_ind]) min = 0;
					else if (std::abs(s[s_ind]) < min) min = std::abs(s[s_ind]);
					run_count++;
				
					if(run_t[t_ind] && (double)run_count < t[t_ind])
						cost += (t[t_ind] - run_count) * min * min;
				}
				
				run_count = 1;
				if(run_t[t_ind]) min = 0;
				else min = std::abs(t[t_ind]);

			}
			t_ind--;
		}
		// moving to left direction
		else if (dir == 'l') {
			if(prev_d == 'l') {
				run_count++;
				if(run_s[s_ind]) min = 0;
				else if (std::abs(s[s_ind]) < min) min = std::abs(s[s_ind]);
			}
			else {
				if(prev_d == 't') {
					run_count++;
					if(run_t[t_ind]) min = 0;
					else if(std::abs(t[t_ind]) < min) min = std::abs(t[t_ind]);

					if(run_s[s_ind] && (double)run_count < s[s_ind])
						cost += (s[s_ind]-run_count)*min*min;

				}

				run_count = 1;
				if(run_s[s_ind]) min = 0;
				else min = std::abs(s[s_ind]);
			}
		
			s_ind--;
		}

		// dir = 's'
		else {
			if(prev_d == 'l') {
				run_count++;
				if (run_s[s_ind]) min = 0;
				else if (std::abs(s[s_ind]) < min) min = std::abs(s[s_ind]);
				if (run_t[t_ind] && (double)run_count < t[t_ind])
					cost += t[t_ind]*min*min;
					
			}
			else if(prev_d == 't') {
				run_count++;
				if(run_t[t_ind]) min = 0;
				else if(std::abs(t[t_ind]) < min) min = std::abs(t[t_ind]);

				if (run_s[s_ind] && (double)run_count < s[s_ind])
					cost += s[s_ind]*min*min;
			}
		
			break;
		}

		prev_d = dir;
	}


	return sqrt(cost);
}


double Mine_LB(std::vector<double> s, std::vector<double> t, bool* rs, bool* rt) {
	int ns = s.size();
	int nt = t.size();
	double **D = new double*[ns+1];
	for(int i = 0; i<=ns; i++)
		D[i] = new double[nt+1];
	

	for(int i = 0; i<=ns; i++)
		D[i][0] = inf;

	for(int i = 0; i<=nt; i++)
		D[0][i] = inf;
	D[0][0] = 0;

	for(int i = 0; i < ns; i++) {
		for(int j = 0; j < nt; j++) {
			double a1 = D[i][j] + (s[i]-t[j])*(s[i]-t[j]);
			if(i > 0 && j > 0) {
				a1 = D[i][j] + LBCases(s[i], t[j], 'd', rs[i], rt[j]);

				if(rs[i-1] && !rt[j-1] && s[i-1] > j) {
					if(!rt[j])
						a1 += (s[i-1] - j)*t[j-1]*t[j-1];
					
				}

				else if(!rs[i] && rt[j] && t[j-1] > i) {
					if(!rs[i]) 
						a1 += (t[j-1] - i)*s[i-1]*s[i-1];
					
				}
			}
			
			double a2 = D[i][j+1] + LBCases(s[i], t[j], 'l', rs[i], rt[j]);
			if(rs[i-1] && !rt[j] && s[i-1] > j) {
				a2 += (s[i-1] - (j+1))*t[j]*t[j];
			}
			
			double a3 = D[i+1][j] + LBCases(s[i], t[j], 't', rs[i], rt[j]);
			if(!rs[i] && rt[j-1] && t[j-1] > i) {
				a3 += (t[j-1] - (i+1))*s[i]*s[i];
			}


			D[i+1][j+1] = MIN(a1,a2,a3);

		}
	}

	double cost = sqrt(D[ns][nt]);
	
	for(int i = 0; i<=ns; i++)
		delete[] D[i];

	delete[] D;
	return cost;
}

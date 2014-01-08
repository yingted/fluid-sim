#include <vector>
#include <utility>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "util.hpp"

void set_boundary(grid& a, double val){
	for (int i = 0; i < a.size(); ++i)
		a[i][0] = a[i][a[i].size()-1] = 0;
	for (int j = 0; j < a[0].size(); ++j)
		a[0][j] = a[a.size()-1][j]= 0;
}

grid diffuse(const grid& a0, double mu, double boundary){
	grid a = a0;
	for (int t = 0; t < 20; ++t){
		for (int i = 1; i+1 < a0.size(); ++i)
			for (int j = 1; j+1 < a0[i].size(); ++j)
				a[i][j] = (a0[i][j] + mu*(
							+ a[i-1][j  ]
							+ a[i+1][j  ]
							+ a[i  ][j-1]
							+ a[i  ][j+1]
							))/(1+4*mu);
		set_boundary(a, boundary);
	}
	return a;
}

grid advect(const grid& a0, const grid& dx, const grid& dy){
	grid a = a0;

	return a;
}

int main(){
	int N = 20, M = 20;
	double mu = .1;
	// coordinates: math-style
	grid dx = make_grid(N, M), dy = dx, p = dx, fx = dx, fy = dx, s = dx;
	for (int i = N/4; i < 3*N/4; ++i)
		for (int j = M/4; j <= 3*M/4; ++j)
			s[i][j] = .1; // smoke in middle quarter
	for (int i = 0; i < M/2; ++i)
		fx[1][i] = .2; // wind on the left
	for (int t = 0; t < 100; ++t){
		// 1. add forces
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j){
				dx[i][j] += fx[i][j];
				dy[i][j] += fy[i][j];
				p [i][j] += s [i][j];
			}

		// 2. diffuse, 3. advect
		grid np  = advect(diffuse(p , mu, 0), dx, dy);
		grid ndx = advect(diffuse(dx, mu, 0), dx, dy);
		grid ndy = advect(diffuse(dy, mu, 0), dx, dy);
		set_boundary(np, 0);

		// copy frame
		p = std::move(np); dx = std::move(ndx); dy = std::move(ndy);

		// write output
		char *name;
		assert(asprintf(&name, "smoke-grid-%04d.ppm", t+1) >= 0);
		std::ofstream f(name);
		free(name);
		double dpeak = 1;
		unsigned short peak = ~0;
		f << "P3\n" << N << ' ' << M << "\n" << peak << "\n";
		for (int j = M-1; j >=0; --j){
			for (int i = 0; i < N; ++i)
				for (const double& val : {dx[i][j], p[i][j]*2-dpeak, dy[i][j]}) // rgb
					f << (unsigned short)(std::max(0., std::min((double)peak, round((val+1)/2/dpeak*peak)))) << ' ';
			f << '\n';
		}
	}
	return 0;
} // vim: set ts=4 sw=4 noet:

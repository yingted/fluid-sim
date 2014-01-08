#include <vector>
#include <utility>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "util.hpp"

grid diffuse(const grid& a0, double mu){
	grid a = a0;

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
	for (int i = 0; i < M; ++i)
		fx[0][i] = .2; // wind on the left
	for (int t = 0; t < 100; ++t){
		// 1. add forces
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j){
				dx[i][j] += fx[i][j];
				dy[i][j] += fy[i][j];
				p [i][j] += s [i][j];
			}

		// 2. diffuse, 3. advect
		grid np  = advect(diffuse(p , mu), dx, dy);
		grid ndx = advect(diffuse(dx, mu), dx, dy);
		grid ndy = advect(diffuse(dy, mu), dx, dy);
		for (int i = 0; i < N; ++i) // clear boundaries
			np[i][0] = 0;
		for (int j = 0; j < M; ++j)
			np[0][j] = 0;

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

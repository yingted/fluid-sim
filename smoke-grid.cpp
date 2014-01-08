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
		a[i][0] = a[i].back() = val;
	fill(a[0]    .begin(), a[0]    .end(), val);
	fill(a.back().begin(), a.back().end(), val);
}

grid diffusion(const grid& a0, double mu, double boundary){
	grid a = a0;
	for (int t = 0; t < 20; ++t){
		for (int i = 1; i+1 < a0.size(); ++i)
			for (int j = 1; j+1 < a0[i].size(); ++j)
				a[i][j] = (a0[i][j]+mu*(
					+a[i-1][j  ]
					+a[i+1][j  ]
					+a[i  ][j-1]
					+a[i  ][j+1]
					))/(1+4*mu);
		set_boundary(a, boundary);
	}
	return std::move(a);
}

grid advection(const grid& a0, const grid& dx, const grid& dy){
	grid a = a0;
	for (int i = 1; i+1 < a0.size(); ++i)
		for (int j = 1; j+1 < a0[i].size(); ++j){
			double pi = std::max(0., std::min(nextafter(a0   .size()-1, 0), i-dx[i][j])), // coords
			       pj = std::max(0., std::min(nextafter(a0[0].size()-1, 0), j-dy[i][j]));
			int ii = floor(pi), ij = floor(pj);
			double s = pi-ii, t = pj-ij;
			a[i][j] =
				(1-s)*((1-t)*a0[ii  ][ij]+t*a0[ii  ][ij+1])
				+  s *((1-t)*a0[ii+1][ij]+t*a0[ii+1][ij+1]);
		}
	set_boundary(a, 0);
	return std::move(a);
}

void project(grid& dx, grid& dy){
	grid div = dx; // divergence = del dot v
	for (int i = 1; i+1 < dx.size(); ++i)
		for (int j = 1; j+1 < dx[i].size(); ++j)
			div[i][j] = (dx[i+1][j]-dx[i-1][j]+dy[i][j+1]-dy[i][j-1])/2;
	set_boundary(div, 0);

	grid grad = make_grid(dx.size(), dx.empty() ? 0 : dx[0].size());
	for (int t = 0; t < 20; ++t){
		for (int i = 1; i+1 < dx.size(); ++i)
			for (int j = 1; j+1 < dx[i].size(); ++j)
				grad[i][j] = (div[i][j]+
					+grad[i-1][j  ]
					+grad[i+1][j  ]
					+grad[i  ][j-1]
					+grad[i  ][j+1]
					)/4;
		set_boundary(grad, 0);
	}

	for (int i = 1; i+1 < dx.size(); ++i)
		for (int j = 1; j+1 < dx[i].size(); ++j){
			dx[i][j] += (grad[i+1][j  ]-grad[i-1][j  ])/2;
			dy[i][j] += (grad[i  ][j+1]-grad[i  ][j-1])/2;
		}

	for (int i = 0; i < dx.size(); ++i)
		dy[i][0] = dy[i].back() = 0;
	std::fill(dx[0]    .begin(), dx[0]    .end(), 0);
	std::fill(dx.back().begin(), dx.back().end(), 0);
}

int main(){
	int N = 50, M = 50;
	double mu = .1;
	// coordinates: math-style
	grid dx = make_grid(N, M), dy = dx, p = dx, fx = dx, fy = dx, s = dx;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			s[N/3+i][M/4+j] = .1;
	for (int j = 0; j < M/2; ++j)
		fx[1][j] = fx[2][j] = .2; // wind on the left
	for (int t = 0; t < 100; ++t){
		// 1. add forces
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j){
				dx[i][j] += fx[i][j];
				dy[i][j] += fy[i][j];
				p [i][j] += s [i][j];
			}

		// 2. diffusion
		dx = diffusion(dx, mu, 0);
		dy = diffusion(dy, mu, 0);
		p  = diffusion(p , mu, 0);
		set_boundary(p, 0);
		project(dx, dy);

		//3. advection
		dx = advection(dx, dx, dy);
		dy = advection(dy, dx, dy);
		p  = advection(p , dx, dy);
		set_boundary(p, 0);
		project(dx, dy);

		// write output
		char *name;
		assert(asprintf(&name, "smoke-grid-%04d.ppm", t+1) >= 0);
		std::ofstream f(name);
		free(name);
		double dpeak = 1;
		//unsigned short peak = ~0;
		unsigned short peak = 255;
		f << "P3\n" << N << ' ' << M << "\n" << peak << "\n";
		for (int j = M-1; j >=0; --j){
			for (int i = 0; i < N; ++i)
				for (const double& val : {dx[i][j], p[i][j]*2-dpeak, dy[i][j]}) // rgb
				//for (const double& val : {p[i][j]*2-dpeak, p[i][j]*2-dpeak, p[i][j]*2-dpeak}) // grey
					f << (unsigned short)(std::max(0., std::min((double)peak, round((val+1)/2/dpeak*peak)))) << ' ';
			f << '\n';
		}
	}
	return 0;
} // vim: set ts=4 sw=4 noet:

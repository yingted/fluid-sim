#include <vector>
#include <utility>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "util.hpp"

#define BOUNDARY_VERTICAL (1)
#define BOUNDARY_HORIZONTAL (2)
#define BOUNDARY_BORDER (BOUNDARY_VERTICAL|BOUNDARY_HORIZONTAL)

void set_boundary(grid& a, double val, int type){
	assert (!a.empty() && !a[0].empty());
	if (type & BOUNDARY_HORIZONTAL)
		for (int i = 0; i < a.size(); ++i)
			a[i][0] = a[i].back() = val; // y = 0
	if (type & BOUNDARY_VERTICAL){
		fill(a[0]    .begin(), a[0]    .end(), val); // x = 0
		fill(a.back().begin(), a.back().end(), val);
	}
}

grid diffusion(const grid& a0, double mu, double boundary, int type){
	assert (!a0.empty() && !a0[0].empty());
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
		set_boundary(a, boundary, type);
	}
	return std::move(a);
}

grid advection(const grid& a0, const grid& dx, const grid& dy, double ox, double oy, int type){
	grid a = a0;
	for (int i = 1; i+1 < a0.size(); ++i)
		for (int j = 1; j+1 < a0[i].size(); ++j){
			double cur_dx = (1-ox)*dx[i][j]+ox*dx[i+1][j  ],
			       cur_dy = (1-oy)*dy[i][j]+oy*dy[i  ][j+1],
			       pi = std::max(0., std::min(nextafter(a0   .size()-1, 0), i-cur_dx)), // coords
			       pj = std::max(0., std::min(nextafter(a0[0].size()-1, 0), j-cur_dy));
			int ii = floor(pi), ij = floor(pj);
			double s = pi-ii, t = pj-ij;
			a[i][j] = (1-s)*((1-t)*a[ii  ][ij]+t*a[ii  ][ij+1])
			          +  s *((1-t)*a[ii+1][ij]+t*a[ii+1][ij+1]);
		}
	set_boundary(a, 0, type);
	return std::move(a);
}

void project(grid& dx, grid& dy){
	grid div = make_grid(dy.size(), dx[0].size()), p = div; // divergence = del dot v
	for (int i = 0; i < div.size(); ++i)
		for (int j = 0; j < div[i].size(); ++j)
			div[i][j] = (dx[i+1][j]-dx[i][j]+dy[i][j+1]-dy[i][j])/2;
	set_boundary(div, 0, BOUNDARY_BORDER);

	for (int t = 0; t < 20; ++t){
		for (int i = 1; i+1 < p.size(); ++i)
			for (int j = 1; j+1 < p.size(); ++j)
				p[i][j] = (div[i][j]+
					+p[i-1][j  ]
					+p[i+1][j  ]
					+p[i  ][j-1]
					+p[i  ][j+1]
					)/4;
		set_boundary(p, 0, BOUNDARY_BORDER);
	}

	for (int i = 1; i+1 < dx.size(); ++i)
		for (int j = 0; j < dx[i].size(); ++j)
			dx[i][j] += (p[i][j]-p[i-1][j])/2;
	for (int i = 0; i < dy.size(); ++i)
		for (int j = 1; j+1 < dy[i].size(); ++j)
			dy[i][j] += (p[i][j]-p[i][j-1])/2;

	set_boundary(dx, 0, BOUNDARY_VERTICAL);
	set_boundary(dy, 0, BOUNDARY_HORIZONTAL);
}

int main(){
	int N = 50, M = 50;
	double mu = .1;
	// coordinates: math-style
	grid dx = make_grid(N+1, M), dy = make_grid(N, M+1), p = make_grid(N, M), fx = dx, fy = dy, s = p;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			s[N/3+i][M/4+j] = .1;
	for (int j = 1; j < M/2; ++j)
		fx[3][j] = .2; // wind on the left
	for (int t = 0; t < 100; ++t){
		// 1. add forces
		dx += fx;
		dy += fy;
		p  += s ;

		// 2. diffusion
		dx = diffusion(dx, mu, 0, BOUNDARY_VERTICAL);
		dy = diffusion(dy, mu, 0, BOUNDARY_HORIZONTAL);
		p  = diffusion(p , mu, 0, BOUNDARY_BORDER);
		set_boundary(p, 0, BOUNDARY_BORDER);
		project(dx, dy);

		//3. advection
		dx = advection(dx, dx, dy,  0, .5, BOUNDARY_VERTICAL);
		dy = advection(dy, dx, dy, .5,  0, BOUNDARY_HORIZONTAL);
		p  = advection(p , dx, dy, .5, .5, BOUNDARY_BORDER);
		set_boundary(p, 0, BOUNDARY_BORDER);
		project(dx, dy);

		// write output
		char *name;
		assert(asprintf(&name, "fluid-grid-%04d.ppm", t+1) >= 0);
		std::ofstream f(name);
		free(name);
		double dpeak = 1;
		//unsigned short peak = ~0;
		unsigned short peak = 255;
		f << "P3\n" << N << ' ' << M << "\n" << peak << "\n";
		for (int j = M-1; j >=0; --j){
			for (int i = 0; i < N; ++i)
				for (const double& val : {.5*(dx[i][j]+dx[i+1][j]), p[i][j]*2-dpeak, .5*(dy[i][j]+dy[i][j+1])}) // rgb
				//for (const double& val : {p[i][j]*2-dpeak, p[i][j]*2-dpeak, p[i][j]*2-dpeak}) // grey
					f << (unsigned short)(std::max(0., std::min((double)peak, round((val+1)/2/dpeak*peak)))) << ' ';
			f << '\n';
		}
	}
	return 0;
} // vim: set ts=4 sw=4 noet:

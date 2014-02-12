#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "util.hpp"
#include <sparse_matrix.h>
#include <pcg_solver.h>

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

double sample(const grid& a, double x, double y){
	double pi = std::max(0., std::min(nextafter(a   .size()-1, 0), x)), // coords
	       pj = std::max(0., std::min(nextafter(a[0].size()-1, 0), y));
	int ii = floor(pi), ij = floor(pj);
	double s = pi-ii, t = pj-ij;
	return (1-s)*((1-t)*a[ii  ][ij]+t*a[ii  ][ij+1])
	       +  s *((1-t)*a[ii+1][ij]+t*a[ii+1][ij+1]);
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
		for (int j = 1; j+1 < a0[i].size(); ++j)
			a[i][j] = sample(a0, i-sample(dx, i+ox, j), j-sample(dy, i, j+oy));
	set_boundary(a, 0, type);
	return std::move(a);
}

void advect(std::vector<double>& mx, std::vector<double>& my, const grid& dx, const grid& dy){
	for (int i = 0; i < mx.size(); ++i){
		const double cur_dx = sample(dx, mx[i], my[i]),
		             cur_dy = sample(dy, mx[i], my[i]);
		mx[i] += cur_dx;
		my[i] += cur_dy;
	}
}

void update_state(std::vector<std::vector<bool> >& state, const std::vector<double> mx, const std::vector<double> my){
	for (std::vector<bool>& row : state)
		for (auto cell : row)
			cell = false;
	for (int i = 0; i < mx.size(); ++i)
		state[std::max(0, std::min((int)state.size()-1, (int)mx[i]))][std::max(0, std::min((int)state[0].size()-1, (int)my[i]))] = true;
}

void project(grid& dx, grid& dy, const std::vector<std::vector<bool> >& state){
	set_boundary(dx, 0, BOUNDARY_VERTICAL);
	set_boundary(dy, 0, BOUNDARY_HORIZONTAL);

	grid div = make_grid(dy.size(), dx[0].size()), p = div; // divergence = del dot v
	for (int i = 0; i < div.size(); ++i)
		for (int j = 0; j < div[i].size(); ++j)
			div[i][j] = (dx[i+1][j]-dx[i][j]+dy[i][j+1]-dy[i][j])/2;

	std::vector<std::vector<unsigned int> >row(p.size(), std::vector<unsigned int>(p[0].size(), -1)); // row for i, j
	unsigned int count = 0;
	for (int i = 0; i < p.size(); ++i)
		for (int j = 0; j < p[i].size(); ++j)
			if (state[i][j])
				row[i][j] = count++;
	std::vector<double>rhs(count); // rhs for row
	SparseMatrix<double>mat(count);

	for (int i = 0; i < p.size(); ++i)
		for (int j = 0; j < p[i].size(); ++j){
			if (!state[i][j])
				continue; // p = 0 by default
			const unsigned int ij = row[i][j];
			std::vector<unsigned int>indices;
			indices.push_back(ij);
			int neighbours = 0;
			if (i > 0){
				if (state[i-1][j])
					indices.push_back(row[i-1][j]);
				++neighbours;
			}
			if (i+1 < p.size()){
				if (state[i+1][j])
					indices.push_back(row[i+1][j]);
				++neighbours;
			}
			if (j > 0){
				if (state[i][j-1])
					indices.push_back(row[i][j-1]);
				++neighbours;
			}
			if (j+1 < p[i].size()){
				if (state[i][j+1])
					indices.push_back(row[i][j+1]);
				++neighbours;
			}
			std::vector<double>values(indices.size(), -1);
			values[0] = neighbours;
			rhs[ij] = div[i][j]*neighbours;
			mat.add_sparse_row(ij, indices, values);
		}

	std::vector<double>result(count);
	double residual;
	int iterations;
	bool success = PCGSolver<double>().solve(mat, rhs, result, residual, iterations);
	std::cerr << "residual = " << residual << " iterations = " << iterations << " success = " << success << std::endl;
	rpc("check_symmetric", mat);

	for (int i = 0; i < p.size(); ++i)
		for (int j = 0; j < p[i].size(); ++j)
			if (state[i][j])
				p[i][j] = result[row[i][j]]; // divisor = neighbours

	for (int i = 1; i+1 < dx.size(); ++i)
		for (int j = 0; j < dx[i].size(); ++j)
			dx[i][j] += p[i][j]-p[i-1][j];
	for (int i = 0; i < dy.size(); ++i)
		for (int j = 1; j+1 < dy[i].size(); ++j)
			dy[i][j] += p[i][j]-p[i][j-1];
}

int main(){
	int N = 50, M = 50;
	double mu = .1, g = -.05;
	// coordinates: math-style
	grid dx = make_grid(N+1, M), dy = make_grid(N, M+1), fx = dx, fy = dy;
	std::vector<double> mx = std::vector<double>(), my = std::vector<double>();
	std::vector<std::vector<bool> > state = std::vector<std::vector<bool> >(N, std::vector<bool>(M));
	//for (int i = 4*(M/4); i < 4*(M/2); ++i)
	//	for (int j = 4*(N/4); j < 4*(3*N/4); ++j){
	for (int i = 4*0; i < 4*M; ++i)
		for (int j = 4*0; j < 4*(N/2); ++j){
			mx.push_back(i*.25+.125*rand()/RAND_MAX);
			my.push_back(j*.25+.125*rand()/RAND_MAX);
		}
	update_state(state, mx, my);
	//for (int j = 1; j < M/2; ++j)
	//	fx[3][j] = .2; // wind on the left
	for (int t = 0; t < 100; ++t){
		// add forces
		dx += fx;
		dy += fy;
		for (int i = 0; i < dy.size(); ++i)
			for (int j = 0; j < dy[0].size(); ++j)
				dy[i][j] += g * state[i][j];
		project(dx, dy, state);

		// advection
		dx = advection(dx, dx, dy,  0, .5, BOUNDARY_VERTICAL);
		dy = advection(dy, dx, dy, .5,  0, BOUNDARY_HORIZONTAL);
		rpc("max_abs", std::string("dx"), dx, std::string("dy"), dy);
		//rpc("print_sum", std::string("state"), state);
		advect(mx, my, dx, dy);
		update_state(state, mx, my);

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
				for (const double& val : {.5*(dx[i][j]+dx[i+1][j]), state[i][j]*2-dpeak, .5*(dy[i][j]+dy[i][j+1])}) // rgb
				//for (const double& val : {p[i][j]*2-dpeak, p[i][j]*2-dpeak, p[i][j]*2-dpeak}) // grey
					f << (unsigned short)(std::max(0., std::min((double)peak, round((val+1)/2/dpeak*peak)))) << ' ';
			f << '\n';
		}
	}
	return 0;
} // vim: set ts=4 sw=4 noet:

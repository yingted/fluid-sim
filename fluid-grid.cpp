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

grid&& diffusion(const grid& a0, double mu, double boundary, int type){
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
		mx[i] += .5*(cur_dx+sample(dx, mx[i]+cur_dx, my[i]+cur_dy));
		my[i] += .5*(cur_dy+sample(dy, mx[i]+cur_dx, my[i]+cur_dy));
	}
}

void update_phi(grid& phi, const std::vector<double>& mx, const std::vector<double>& my, const int padding=2){
	for (auto& row : phi)
		for (auto& cell : row)
			cell = hypot(padding, padding);
	const double radius=1.02*sqrt(.5);
	for (int i = 0; i < mx.size(); ++i)
		for (int j = std::max(0, ((int)mx[i])-padding); j <= std::min((int)phi.size()-1, ((int)mx[i])+padding); ++j)
			for (int k = std::max(0, ((int)my[i])-padding); k <= std::min((int)phi[0].size()-1, ((int)my[i])+padding); ++k)
				phi[j][k] = std::min(phi[j][k], hypot(mx[i]-(j+.5), my[i]-(k+.5))-radius);
	rpc("update_phi", phi, mx, my);
}

void project(grid& dx, grid& dy, const grid& phi){
	set_boundary(dx, 0, BOUNDARY_VERTICAL);
	set_boundary(dy, 0, BOUNDARY_HORIZONTAL);

	grid div = make_grid<double>(dy.size(), dx[0].size()), p = div; // divergence = del dot v
	for (int i = 0; i < div.size(); ++i)
		for (int j = 0; j < div[i].size(); ++j)
			div[i][j] = dx[i+1][j]-dx[i][j]+dy[i][j+1]-dy[i][j]; // TODO check boundary?

	std::vector<std::vector<unsigned int> >row(p.size(), std::vector<unsigned int>(p[0].size(), -1)); // row for i, j
	unsigned int count = 0;
	for (int i = 0; i < p.size(); ++i)
		for (int j = 0; j < p[i].size(); ++j)
			if (phi[i][j] < 0)
				row[i][j] = count++;
	std::vector<double>rhs(count); // rhs for row
	SparseMatrix<double>mat(count);

#define THETA(i2,j2) (\
	(phi[i][j]<0)==(phi[(i2)][(j2)]<0)?\
		1:\
		std::max(1e-2,\
			phi[i][j]<0?\
				phi[(i2)][(j2)]/(phi[(i2)][(j2)]-phi[i][j]):\
				phi[i][j]/(phi[i][j]-phi[(i2)][(j2)])\
		)\
)
	for (int i = 0; i < p.size(); ++i)
		for (int j = 0; j < p[i].size(); ++j){
			if (!(phi[i][j] < 0))
				continue; // p = 0 by default
			const unsigned int ij = row[i][j];
			std::vector<unsigned int>indices;
			indices.push_back(ij);
			double neighbours = 0;
#define CHECK(i2,j2) do{\
	if (phi[(i2)][(j2)] < 0){\
		indices.push_back(row[(i2)][(j2)]);\
		neighbours+=1;\
	}else\
		neighbours+=1/THETA((i2),(j2));\
}while(0)
			if (i > 0)
				CHECK(i-1, j);
			if (i+1 < p.size())
				CHECK(i+1, j);
			if (j > 0)
				CHECK(i, j-1);
			if (j+1 < p[i].size())
				CHECK(i, j+1);
#undef CHECK
			std::vector<double>values(indices.size(), -1);
			values[0] = neighbours;
			rhs[ij] = div[i][j];
			mat.add_sparse_row(ij, indices, values);
		}

	std::vector<double>result(count);
	double residual;
	int iterations;
	PCGSolver<double>solver;
#ifdef NDEBUG
	solver.set_solver_parameters(1e-5, 1000, .97, .25);
#endif
	bool success = !count || solver.solve(mat, rhs, result, residual, iterations);
	std::cerr << "residual = " << residual << " iterations = " << iterations << " success = " << success << std::endl;
	rpc("check_symmetric", mat);

	for (int i = 0; i < p.size(); ++i)
		for (int j = 0; j < p[i].size(); ++j)
			if (phi[i][j] < 0)
				p[i][j] = result[row[i][j]]; // divisor = neighbours

#define DU(i2,j2) ((p[i][j]-p[(i2)][(j2)])/THETA((i2),(j2)))
	for (int i = 1; i+1 < dx.size(); ++i)
		for (int j = 0; j < dx[i].size(); ++j)
			dx[i][j] += DU(i-1, j);
	for (int i = 0; i < dy.size(); ++i)
		for (int j = 1; j+1 < dy[i].size(); ++j)
			dy[i][j] += DU(i, j-1);
#undef DU
#undef THETA
}

void dilate(grid& data, std::vector<std::vector<bool> >& mask, const double boundary=0, unsigned char layers=3){
	std::vector<std::pair<int, int> >cur, next;
	for (int i = 0; i < mask.size(); ++i)
		for (int j = 0; j < mask[i].size(); ++j)
			if (!mask[i][j] && (
					(i > 0 && mask[i-1][j]) ||
					(j > 0 && mask[i][j-1]) ||
					(i+1 < mask.size() && mask[i+1][j]) ||
					(j+1 < mask[i].size() && mask[i][j+1])
				))
				next.push_back(std::make_pair(i, j));

	for (; layers--; swap(cur, next)){
		for (const std::pair<int, int>& p : cur){
			if (mask[p.first][p.second])
				continue;
			double sum = 0;
			int count = 0;
#define CHECK(i,j) do{\
	if (mask[p.first+(i)][p.second+(j)]){\
		sum += data[p.first+(i)][p.second+(j)];\
		++count;\
	}else\
		next.push_back(std::make_pair(p.first+(i),p.second+(j)));\
}while(0)
			if (p.first)
				CHECK(-1, 0);
			if (p.second)
				CHECK(0, -1);
			if (p.first+1 < mask.size())
				CHECK(+1, 0);
			if (p.second+1 < mask[p.first].size())
				CHECK(0, +1);
#undef CHECK
			data[p.first][p.second] = sum/count;
		}
		for (const std::pair<int, int>& p : cur)
			mask[p.first][p.second] = true;
	}

	for (int i = 0; i < mask.size(); ++i)
		for (int j = 0; j < mask[i].size(); ++j)
			if (!mask[i][j])
				data[i][j] = boundary;
}

int main(){
#ifdef NDEBUG
	int N = 500, M = 500, T = 1000;
	double g = -.005;
#else
	int N = 30, M = 30, T = 100;
	double g = -.05;
#endif
	// coordinates: math-style
	grid dx = make_grid<double>(N+1, M), dy = make_grid<double>(N, M+1), fx = dx, fy = dy;
	std::vector<double> mx = std::vector<double>(), my = std::vector<double>();
	grid phi = make_grid<double>(N, M);
	//for (int i = 4*(M/4); i < 4*(M/2); ++i)
	//	for (int j = 4*(N/4); j < 4*(3*N/4); ++j){
	for (int i = 4*0; i < 4*M; ++i)
		for (int j = 4*0; j < 4*(N/2); ++j){
	//for (int i = 4*(M/2); i <= 4*(M/2); ++i)
	//	for (int j = 4*(N/2); j <= 4*(N/2); ++j){
			mx.push_back(i*.25+.125*rand()/RAND_MAX);
			my.push_back(j*.25+.125*rand()/RAND_MAX);
			if (my.back() > .25*N+.25*mx.back()*N/M){
				mx.pop_back();
				my.pop_back();
			}
		}
	update_phi(phi, mx, my);
	//for (int j = 1; j < M/2; ++j)
	//	fx[3][j] = .2; // wind on the left
	for (int t = 0; t < T; ++t){
		// add forces
		{
			dx += fx;
			dy += fy;
			for (int i = 0; i < dy.size(); ++i)
				for (int j = 0; j < dy[0].size(); ++j)
					if ((j > 0 && phi[i][j-1] < 0) || (j < phi[0].size() && phi[i][j] < 0))
						dy[i][j] += g; // flow in
		}

		// velocity boundary conditions
		{
			std::vector<std::vector<bool> > mask_dx = make_grid<bool>(N+1, M), mask_dy = make_grid<bool>(N, M+1);
			for (int i = 0; i < phi.size(); ++i)
				for (int j = 0; j < phi[0].size(); ++j)
					if (phi[i][j] < 0){
						mask_dx[i][j] = mask_dx[i+1][j] = true; // give 1 cell leeway
						mask_dy[i][j] = mask_dy[i][j+1] = true; // for boundary condition
					}
			for (int i = 0; i < dx.size(); ++i)
				for (int j = 0; j < dx.size(); ++j)
					if (!mask_dx[i][j])
						dx[i][j] = 0;
			for (int i = 0; i < dy.size(); ++i)
				for (int j = 0; j < dy.size(); ++j)
					if (!mask_dy[i][j])
						dy[i][j] = 0;
			project(dx, dy, phi);
			//rpc("deciles", std::string("dx"), dx);
			//rpc("deciles", std::string("dy"), dy);
			dilate(dx, mask_dx);
			dilate(dy, mask_dy);
		}

		// advection
		{
			grid ndx = advection(dx, dx, dy,  0, .5, BOUNDARY_VERTICAL);
			dy = advection(dy, dx, dy, .5,  0, BOUNDARY_HORIZONTAL);
			dx = std::move(ndx);
			rpc("max_abs", std::string("dx"), dx, std::string("dy"), dy);
			advect(mx, my, dx, dy);
			update_phi(phi, mx, my);
		}

		// write output
		{
			char *name;
			{
				int bytes = asprintf(&name, "fluid-grid-%04d.ppm", t+1);
				assert(bytes >= 0);
			}
			std::ofstream f(name);
			free(name);
			double dpeak = 1;
			//unsigned short peak = ~0;
			unsigned short peak = 255;
			f << "P3\n" << N << ' ' << M << "\n" << peak << "\n";
			for (int j = M-1; j >=0; --j){
				for (int i = 0; i < N; ++i)
					for (const double& val : {.5*(dx[i][j]+dx[i+1][j]), phi[i][j]*2-dpeak, .5*(dy[i][j]+dy[i][j+1])}) // rgb
					//for (const double& val : {p[i][j]*2-dpeak, p[i][j]*2-dpeak, p[i][j]*2-dpeak}) // grey
						f << (unsigned short)(std::max(0., std::min((double)peak, round((val+1)/2/dpeak*peak)))) << ' ';
				f << '\n';
			}
		}
	}
	return 0;
} // vim: set ts=4 sw=4 noet:

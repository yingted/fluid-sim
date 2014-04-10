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
#include <queue>
#include "util.hpp"
#include <sparse_matrix.h>
#include <pcg_solver.h>

#define BOUNDARY_NONE (0)
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

template<typename T>
int clamp(T x, int start, int stop){
	return floor(std::max<T>(start, std::min<T>(nextafter(stop, start), x)));
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
	return a;
}

grid advection(const grid& a0, const grid& dx, const grid& dy, double ox, double oy, int type){
	grid a = a0;
	for (int i = 0; i < a0.size(); ++i)
		for (int j = 0; j < a0[i].size(); ++j)
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

double phi_theta(double a, double b){
	return a < 0 ?
		b < 0 ?
			1 :
			a/(a-b) :
		b < 0 ?
			b/(b-a) :
			0;
}

#define THETA(i2,j2) (std::max(1e-2, phi_theta(phi[i][j], phi[(i2)][(j2)])))
void interpolate_surface(const grid& phi, std::vector<double>& bx, std::vector<double>& by){
	bx.clear();
	by.clear();
	for (int i = 0; i < phi.size(); ++i)
		for (int j = 0; j < phi[i].size(); ++j)
			if (phi[i][j] < 0){
#define CHECK(i2,j2) do{\
	if(!(phi[(i2)][(j2)] < 0)){\
		const double theta = THETA((i2),(j2));\
		bx.push_back(i*(1-theta)+(i2)*theta+.5);\
		by.push_back(j*(1-theta)+(j2)*theta+.5);\
	}\
}while(0)
				if (i > 0)
					CHECK(i-1, j);
				if (i+1 < phi.size())
					CHECK(i+1, j);
				if (j > 0)
					CHECK(i, j-1);
				if (j+1 < phi[i].size())
					CHECK(i, j+1);
#undef CHECK
			}
}

void check_phi(const grid& phi){
	return;
	for (int i = 0; i < phi.size(); ++i)
		for (int j = 0; j < phi[0].size(); ++j)
			for (int k = 0; k < phi.size(); ++k)
				for (int l = 0; l < phi[0].size(); ++l)
					assert(fabs(phi[i][j]-phi[k][l]) <= hypot(i-k, j-l)+2+1e-5);
}

void redistance(grid& phi, const std::vector<double>& bx, const std::vector<double>& by){
	const grid phi0 = phi;
	std::vector<std::vector<bool> >seen = make_grid<bool>(phi.size(), phi[0].size());
	std::vector<std::vector<int> >ancestor = make_grid<int>(phi.size(), phi[0].size());
	typedef std::pair<double, std::pair<int, int> >vertex_t;
	std::priority_queue<vertex_t, std::vector<vertex_t>, std::greater<vertex_t> >q;
	for (int i = 0; i < phi.size(); ++i)
		for (int j = 0; j < phi[0].size(); ++j)
			phi[i][j] = std::numeric_limits<double>::infinity();
	for (int i = 0; i < bx.size(); ++i){
		const int x_lo = clamp(bx[i]-.5, -1, phi.size()),
		          y_lo = clamp(by[i]-.5, -1, phi[0].size());
		for (int cx = max(0, x_lo); cx <= std::min<int>(phi.size()-1, x_lo+1); ++cx)
			for (int cy = max(0, y_lo); cy <= std::min<int>(phi[0].size()-1, y_lo+1); ++cy){
				const double d = hypot(cx+.5-bx[i], cy+.5-by[i]);
				if (d < fabs(phi[cx][cy])){
					phi[cx][cy] = copysign(d, sample(phi0, cx, cy));
					ancestor[cx][cy] = i;
					q.push(std::make_pair(d, std::make_pair(cx, cy)));
				}
			}
	}
	while (!q.empty()){
		const double d = q.top().first;
		const int cx = q.top().second.first,
		          cy = q.top().second.second;
		q.pop();
		if (d != fabs(phi[cx][cy]))
			continue;
		seen[cx][cy] = true;
		const int anc = ancestor[cx][cy];
		const double x = bx[anc],
		             y = by[anc];
#define CHECK(i2,j2) do{\
	const int i = (i2), j = (j2);\
	if (!seen[i][j]){\
		const double d2 = hypot(i+.5-x, j+.5-y);\
		if (d2 < fabs(phi[i][j])){\
			phi[i][j] = copysign(d2, sample(phi0, i, j));\
			ancestor[i][j] = anc;\
			q.push(std::make_pair(d2, std::make_pair(i, j)));\
		}\
	}\
}while(0)
		if (cx)
			CHECK(cx-1, cy);
		if (cx+1 < phi.size())
			CHECK(cx+1, cy);
		if (cy)
			CHECK(cx, cy-1);
		if (cy+1 < phi[0].size())
			CHECK(cx, cy+1);
#undef CHECK
	}
	check_phi(phi);
}

void project(const grid& solid_phi, grid& dx, grid& dy, const grid& phi){
	set_boundary(dx, 0, BOUNDARY_VERTICAL);
	set_boundary(dy, 0, BOUNDARY_HORIZONTAL);

#define THETA_X(i2,j2) (1-phi_theta(solid_phi[(i2)][(j2)], solid_phi[(i2)][(j2)+1]))
#define THETA_Y(i2,j2) (1-phi_theta(solid_phi[(i2)][(j2)], solid_phi[(i2)+1][(j2)]))
	grid div = make_grid<double>(dy.size(), dx[0].size()), p = div; // divergence = del dot v
	for (int i = 0; i < div.size(); ++i)
		for (int j = 0; j < div[i].size(); ++j)
			div[i][j] = THETA_X(i+1, j)*dx[i+1][j]-THETA_X(i, j)*dx[i][j]+THETA_Y(i, j+1)*dy[i][j+1]-THETA_Y(i, j)*dy[i][j]; // TODO check boundary?

	std::vector<std::vector<unsigned int> >row(p.size(), std::vector<unsigned int>(p[0].size(), -1)); // row for i, j
	unsigned int count = 0;
	for (int i = 0; i < p.size(); ++i)
		for (int j = 0; j < p[i].size(); ++j)
			if (phi[i][j] < 0)
				row[i][j] = count++;
	std::vector<double>rhs(count); // rhs for row
	SparseMatrix<double>mat(count);

	for (int i = 0; i < p.size(); ++i)
		for (int j = 0; j < p[i].size(); ++j){
			if (!(phi[i][j] < 0))
				continue; // p = 0 by default
			const unsigned int ij = row[i][j];
			std::vector<unsigned int>indices;
			std::vector<double>values;
			double neighbours = 0;
#define CHECK(i2,j2,w) do{\
	double coef = (w);\
	if (coef){\
		if (phi[(i2)][(j2)] < 0){\
			values.push_back(-coef);\
			indices.push_back(row[(i2)][(j2)]);\
		}else\
			coef /= THETA((i2),(j2));\
		neighbours += coef;\
	}\
}while(0)
			if (i > 0)
				CHECK(i-1, j, THETA_X(i, j));
			if (i+1 < p.size())
				CHECK(i+1, j, THETA_X(i+1, j));
			if (j > 0)
				CHECK(i, j-1, THETA_Y(i, j));
			if (j+1 < p[i].size())
				CHECK(i, j+1, THETA_Y(i, j+1));
#undef CHECK
			values.push_back(neighbours);
			indices.push_back(ij);
			rhs[ij] = div[i][j];
			if (neighbours)
				mat.add_sparse_row(ij, indices, values);
			//std::cerr << "mat[" << i << "][" << j << "]: values = " << values << ", rhs = " << rhs[ij] << std::endl;
			assert(values == values);
		}

	std::vector<double>result(count);
	double residual;
	int iterations;
	PCGSolver<double>solver;
#ifdef NDEBUG
	solver.set_solver_parameters(1e-5, 1000, .97, .25);
#endif
	rpc("max_abs", std::string("dx"), dx, std::string("dy"), dy, std::string("div"), div, std::string("rhs"), rhs);
	//std::cerr << "mat = " << mat << " rhs = " << rhs << std::endl;
	assert(rhs == rhs);
	bool success = !count || solver.solve(mat, rhs, result, residual, iterations);
	std::cerr << "residual = " << residual << " iterations = " << iterations << " success = " << success << std::endl;
	rpc("check_symmetric", mat);

	for (int i = 0; i < p.size(); ++i)
		for (int j = 0; j < p[i].size(); ++j)
			if (phi[i][j] < 0){
				p[i][j] = result[row[i][j]]; // divisor = neighbours
				//std::cerr << "p[" << i << "][" << j << "] = " << p[i][j] << std::endl;
			}
#define DU(i2,j2) ((p[i][j]-p[(i2)][(j2)])/THETA((i2),(j2)))
	for (int i = 1; i+1 < dx.size(); ++i)
		for (int j = 0; j < dx[i].size(); ++j)
			dx[i][j] += DU(i-1, j);
	for (int i = 0; i < dy.size(); ++i)
		for (int j = 1; j+1 < dy[i].size(); ++j)
			dy[i][j] += DU(i, j-1);
#undef DU
	//std::cerr << "dx = " << dx << std::endl;
	//std::cerr << "dy = " << dy << std::endl;
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

template<typename T>
void mask(T& a, const bool& mask_a, const T val){
	if (!mask_a)
		a = val;
}

template<typename T, typename B, typename V>
void mask(std::vector<T>& a, const std::vector<B>& mask_a, const V val){ // rotate by pi
	for (int i = 0; i < a.size(); ++i)
		mask(a[i], mask_a[i], val);
}

template<typename T>
void flip(std::vector<std::vector<T> >& a){ // rotate by pi
	reverse(a.begin(), a.end());
	for (std::vector<T>& row : a)
		reverse(row.begin(), row.end());
}

template<typename T>
void linear(T& a, T m, T b=0){ // transform by mx+b
	a = m*a + b;
}

template<typename T, typename E>
void linear(std::vector<E>& a, T m, T b=0){
	for (E& row : a)
		linear(row, m, b);
}

int main(){
	//srand(time(NULL));
#ifdef NDEBUG
	const int N = 500, M = 500, T = 1000, redistance_period = 1;
	const double gx = 0, gy = -.005;
#else
	const int N = 50, M = 50, T = 200, redistance_period = 1;
	//const int N = 10, M = 10, T = 1, redistance_period = 1;
	const double gx = 0, gy = -.05, mu = .1;
#endif
	// coordinates: math-style
	grid dx = make_grid<double>(N+1, M), dy = make_grid<double>(N, M+1), fx = dx, fy = dy;
	grid phi = make_grid<double>(N, M), solid_phi = make_grid<double>(N+1, M+1);
	for (int i = 0; i <= M; ++i)
		for (int j = 0; j <= N; ++j)
			solid_phi[i][j] = min(M, N)/2-hypot(i+.5-M/2, j+.5-N/2);
	check_phi(solid_phi);
	for (int i = 0; i < M; ++i)
		for (int j = 0; j < N; ++j)
			phi[i][j] = std::max<double>(-solid_phi[i][j],
				//j-N/2
				//j-N/2+.25*i
				i-M/2
			);
	{
		std::vector<double>bx, by;
		interpolate_surface(phi, bx, by);
		rpc("draw", solid_phi, dx, dy, phi, bx, by);
	}
	//for (int j = 1; j < M/2; ++j)
	//	fx[3][j] = .2; // wind on the left
	for (int t = 0; t < T; ++t){
		// add forces
		{
			dx += fx;
			dy += fy;
			for (int i = 1; i+1 < dx.size(); ++i)
				for (int j = 0; j < dx[0].size(); ++j)
					if (phi[i-1][j] < 0 || phi[i][j] < 0) // ignore theta = epsilon
						dx[i][j] += gx*THETA(i-1, j); // flow in
			for (int i = 0; i < dy.size(); ++i)
				for (int j = 1; j+1 < dy[0].size(); ++j)
					if (phi[i][j-1] < 0 || phi[i][j] < 0) // ignore theta = epsilon
						dy[i][j] += gy*THETA(i, j-1); // flow in
		}

		// velocity boundary conditions
		{
			std::vector<std::vector<bool> > mask_dx = make_grid<bool>(N+1, M), mask_dy = make_grid<bool>(N, M+1);
			for (int i = 0; i < phi.size(); ++i)
				for (int j = 0; j < phi[0].size(); ++j)
					if (phi[i][j] < 0 && (
							solid_phi[i][j] >= 0 ||
							solid_phi[i][j+1] >= 0 ||
							solid_phi[i+1][j] >= 0 ||
							solid_phi[i+1][j+1] >= 0
						)){
						mask_dx[i][j] = mask_dx[i+1][j] = true;
						mask_dy[i][j] = mask_dy[i][j+1] = true;
					}

			mask(dx, mask_dx, 0.);
			mask(dy, mask_dy, 0.);
			project(solid_phi, dx, dy, phi);
			rpc("max_abs", std::string("dx"), dx, std::string("dy"), dy);
			//rpc("deciles", std::string("dx"), dx);
			//rpc("deciles", std::string("dy"), dy);
			dilate(dx, mask_dx);
			dilate(dy, mask_dy);

			grid ndx = dx;
			for (int i = 0; i < dx.size(); ++i)
				for (int j = 0; j < dx[0].size(); ++j)
					if (THETA_X(i, j) == 0 && mask_dx[i][j]){
						const int i_lo = max(0, i-1), i_hi = min(((int)dx.size())-1, i+1);
						double plus = solid_phi[i_hi][j+1]-solid_phi[i_lo][j],
						      minus = solid_phi[i_hi][j]-solid_phi[i_lo][j+1],
						         nx = (plus+minus)/(i_hi-i_lo),
						         ny = plus-minus,
						         ux = dx[i][j],
						         uy = 0,
						         nn = nx*nx+ny*ny,
						 neighbours = 0;
#define CHECK(i2,j2) do\
	if (mask_dy[(i2)][(j2)]){\
		uy += dy[(i2)][(j2)];\
		neighbours += 1;\
	}\
while(0)
						CHECK(i_lo, j);
						CHECK(i_lo, j+1);
						CHECK(i_hi-1, j);
						CHECK(i_hi-1, j+1);
#undef CHECK
						uy /= neighbours;
						if (nn)
							ndx[i][j] = ux-(nx*ux+ny*uy)/nn*nx;
					}
			for (int i = 0; i < dy.size(); ++i)
				for (int j = 0; j < dy[0].size(); ++j)
					if (THETA_Y(i, j) == 0 && mask_dy[i][j]){
						const int j_lo = max(0, j-1), j_hi = min(((int)dy[0].size())-1, j+1);
						double plus = solid_phi[i+1][j_hi]-solid_phi[i][j_lo],
						      minus = solid_phi[i+1][j_lo]-solid_phi[i][j_hi],
						         nx = plus+minus,
						         ny = (plus-minus)/(j_hi-j_lo),
						         ux = 0,
						         uy = dy[i][j],
						         nn = nx*nx+ny*ny,
						 neighbours = 0;
#define CHECK(i2,j2) do\
	if (mask_dx[(i2)][(j2)]){\
		ux += dx[(i2)][(j2)];\
		neighbours += 1;\
	}\
while(0)
						CHECK(i, j_lo);
						CHECK(i+1, j_lo);
						CHECK(i, j_hi-1);
						CHECK(i+1, j_hi-1);
#undef CHECK
						ux /= neighbours;
						if (nn)
							dy[i][j] = uy-(nx*ux+ny*uy)/nn*ny;
					}
			dx = std::move(ndx);
		}

		// advection
		{
			std::vector<double>bx, by;
			grid ndx = advection(dx, dx, dy,  0, .5, BOUNDARY_VERTICAL);
			dy = advection(dy, dx, dy, .5,  0, BOUNDARY_HORIZONTAL);
			dx = std::move(ndx);
			//rpc("max_abs", std::string("dx"), dx, std::string("dy"), dy);
			phi = advection(phi, dx, dy, .5, .5, BOUNDARY_NONE);
			interpolate_surface(phi, bx, by);
			if ((t+1)%redistance_period == 0)
				redistance(phi, bx, by);
			rpc("draw", solid_phi, dx, dy, phi, bx, by);
		}

#if 0
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
#endif
	}
#if 0
flip(dx);
flip(dy);
flip(phi);
linear(g, -1.);
linear(dx, -1.);
linear(dy, -1.);
#endif
	return 0;
} // vim: set ts=4 sw=4 noet:

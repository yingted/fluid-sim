#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <algorithm>
#include <queue>
#include <map>
#include "util.hpp"
#include <sparse_matrix.h>
#include <pcg_solver.h>

double phi_theta(double a, double b){
	return a < 0 ?
		b < 0 ?
			1 :
			a/(a-b) :
		b < 0 ?
			b/(b-a) :
			0;
}

bool barycentric(double px, double py, double qx, double qy, double sx, double sy,
	double& D, double& Dpl, double& Dql, double& Drl){ // of sx wrt (r = 0), p, q
	D   = py*qx-px*qy; // determinant
	Dpl = qx*sy-qy*sx;
	Dql = py*sx-px*sy;
	Drl = D-Dpl-Dql;
	return D < 0 ? // check if D?l/D \in [0, 1] without division
		D <= Dpl && Dpl <= 0 &&
		D <= Dql && Dql <= 0: // rl condition checks if s is too far, and should always be true
		0 <= Dpl && Dpl <= D &&
		0 <= Dql && Dql <= D;
}

struct quad{ // NULL is the infinite cell
	quad *neighbour[4]; // right up left down
	quad *child[4]; // NE NW SW SE
	const quad *parent;
	const int index;
	const double r, x, y;
	const static int cos[4], sin[4];
	double u[4], solid_phi, phi, dx[8], dy[8], cell_dx, cell_dy;
	void copy_from(const quad *o){
		solid_phi = o->solid_phi;
		phi = o->phi;
		for (int i = 0; i < 4; ++i)
			u[i] = o->u[i];
		for (int i = 0; i < 8; ++i){
			dx[i] = o->dx[i];
			dy[i] = o->dy[i];
		}
		cell_dx = o->cell_dx;
		cell_dy = o->cell_dy;
	}
	quad ghost_cell(){
		quad ret(nan("ghost"), nan("ghost"), nan("ghost"));
		ret.solid_phi = -solid_phi;
		ret.phi = phi;
		for (int i = 0; i < 4; ++i)
			ret.u[i] = u[i];
		for (int i = 0; i < 8; ++i){
			ret.dx[i] = -dx[i];
			ret.dy[i] = -dy[i];
		}
		ret.cell_dx = -cell_dx;
		ret.cell_dy = -cell_dy;
		return ret;
	}
	quad(double x, double y, double r) : parent(NULL), index(-1), x(x), y(y), r(r), neighbour(), child(){}
	quad(quad *parent, int index) : parent(parent), index(index), neighbour(), child(),
		r(.5*parent->r),
		x(parent->x+r*(cos[index]-sin[index])), // add pi/4
		y(parent->y+r*(cos[index]+sin[index])){
		assert(parent != NULL);
		assert(0 <= index && index < 4);
		for (int i = 0; i < 4; ++i){
			assert((i == index || i == (index+1)%4) == ((cos[index]-sin[index])*cos[i]+(cos[index]+sin[index])*sin[i] > 0));
			if (i == index || i == (index+1)%4){ // positive dot product, shared wall
				neighbour[i] = parent->neighbour[i]; // larger neighbour
				if (neighbour[i] && neighbour[i]->child[0]) // child 1+i+(i-index)
					neighbour[i] = neighbour[i]->child[(3*index+2*i+1)%4];
			}else{ // same size neighbour
				if (i == (index+2)%4) // next
					neighbour[i] = parent->child[(index+1)%4];
				else if (i == (index+3)%4) // prev
					neighbour[i] = parent->child[(index+3)%4];
				assert(!neighbour[i] || neighbour[i]->r == r);
			}
			if (neighbour[i] && neighbour[i]->r == r)
				neighbour[i]->neighbour[(i+2)%4] = this;
		}
	}
	void check_relations(){
		assert(!parent == (index == -1));
		if (parent){
			assert(0 <= index && index < 4);
			assert(parent->child[index] == this);
			assert(r == .5*parent->r);
			assert(x == parent->x+r*(cos[index]-sin[index]));
			assert(y == parent->y+r*(cos[index]+sin[index]));
		}
		for (int i = 0; i < 4; ++i)
			if (neighbour[i]){
				assert(neighbour[i]->r >= r);
				if (neighbour[i]->r > r){
					assert(r == .5*neighbour[i]->r);
					for (int j = 0; j < 4; ++j)
						assert(neighbour[i]->child[j] == NULL);
				}
				const double big = neighbour[i]->r+r,
				           small = neighbour[i]->r-r;
				for (int j = 0; j < 4; ++j)
					for (int k = -1; k <= 1; k += 2)
						if (neighbour[i]->x == x+(cos[j]*big-sin[j]*small*k) &&
							neighbour[i]->y == y+(sin[j]*big+cos[j]*small*k))
							goto found;
				assert(false); // not found
found:;
				assert(neighbour[i]->neighbour[(i+2)%4] == (neighbour[i]->r == r ? this : parent));
			}
		for (int i = 0; i < 4; ++i){
			assert(!child[0] == !child[i]);
			if (child[i])
				child[i]->check_relations();
		}
	}
	bool contains(double px, double py)const{
		return x-r <= px && px <= x+r && y-r <= py && py <= y+r;
	}
	const quad *query(double px, double py)const{
		if (!child[0])
			return this;
		int i;
		if (px >= x)
			if (py >= y)
				i = 0;
			else
				i = 3;
		else if (py >= y)
			i = 1;
		else
			i = 2;
		return child[i]->query(px, py);
	}
	quad *query(double px, double py){
		return const_cast<quad*>(static_cast<const quad&>(*this).query(px, py));
	}
	double dist(double px, double py)const{
		return hypot(px-x, py-y);
	}
	double dist2(double px, double py)const{
		return (px-x)*(px-x)+(py-y)*(py-y);
	}
	double power(double px, double py)const{
		return dist2(px, py)-2*r*r;
	}
	const quad *query_nn(double px, double py)const{
		assert(contains(px, py));
		if (child[0])
			return query(px, py)->query_nn(px, py);
		double best = dist2(px, py);
		const quad *ret = this;
		auto test = [&best, &ret, px, py](const quad *const n){
			const double d2 = n->dist2(px, py);
			if (d2 < best){
				best = d2;
				ret = n;
			}
		};
		for (int i = 0; i < 4; ++i){
			const quad *const n = neighbour[i];
			if (n){
				test(n);
				if (n->child[0]){
					test(n->child[(i+1)%4]);
					test(n->child[(i+2)%4]);
				}
			}
		}
		return ret;
	}
	const quad *query_power(double px, double py)const{
		if (child[0])
			return query(px, py)->query_power(px, py);
		double best = power(px, py);
		const quad * ret = this;
		for (int i = 0; i < 4; ++i){
			const quad *const n = neighbour[i];
			if (n && !n->child[0]){
				const double d2 = n->power(px, py);
				if (d2 < best){
					best = d2;
					ret = n;
				}
			}
		}
		assert(best <= 0);
		return ret;
	}
	double query_sample(double sx, double sy, size_t offset)const{
		return query(sx, sy)->sample(sx, sy, offset);
	}
	void query_sample_u(double sx, double sy, double& dx, double& dy)const{
		query_power(sx, sy)->sample_u(sx, sy, dx, dy);
	}
	double& field(size_t offset){
		return *(double*)(((char*)(this))+offset);
	}
	double sample(double sx, double sy, size_t offset)const{
		assert(!child[0]);
		sx = max(-r, min(r, sx-x));
		sy = max(-r, min(r, sy-y));
		static double ret = nan("not found");
		const bool failed = (const_cast<quad*>(this))->visit_my_vertices([sx, sy, offset](quad *c, quad *p, quad *pq, quad *q){
			double lx = -1, ly = -1;
			if (!pq)
				goto barycentric;
			if ((p->x == c->x && q->x == pq->x) || (q->x == c->x && p->x == pq->x)){
				lx = sx/(pq->x-c->x);
				const quad *const s = p->x == c->x ? p : q,
					       *const t = p->x != c->x ? p : q;
				if (0 <= lx && lx <= 1)
					ly = (sy-lx*(t->y-c->y))/((1-lx)*(s->y-c->y)+lx*(pq->y-t->y));
			}else if((p->y == c->y && q->y == pq->y) || (q->y == c->y && p->y == pq->y)){
				ly = sy/(pq->y-c->y);
				const quad *const s = p->y == c->y ? p : q,
					       *const t = p->y != c->y ? p : q;
				if (0 <= ly && ly <= 1)
					lx = (sx-ly*(t->x-c->x))/((1-ly)*(s->x-c->x)+ly*(pq->x-t->x));
			}else
				goto barycentric;
			if (0 <= lx && lx <= 1 &&
				0 <= ly && ly <= 1){
				ret = lx*(1-ly)*p->field(offset)+lx*ly*q->field(offset)+(1-lx)*ly*pq->field(offset)+(1-lx)*(1-ly)*c->field(offset);
				return false;
			}
			return true;
barycentric:
			double D, Dpl, Dql, Drl;
			if (barycentric(p->x-c->x, p->y-c->y, q->x-c->x, q->y-c->y, sx, sy, D, Dpl, Dql, Drl)){
				ret = (Dpl*p->field(offset)+Dql*q->field(offset)+Drl*c->field(offset))/D;
				return false;
			}
			return true;
		});
		assert(!failed);
		return ret;
	}
	static void face_endpoints(quad *c, quad *q, double& px, double& py, double& qx, double& qy){
		if(q->r < c->r){
			face_endpoints(q, c, qx, qy, px, py);
			qx += q->x-c->x;
			px += q->x-c->x;
			qy += q->y-c->y;
			py += q->y-c->y;
			return;
		}
		double cx = q->x-c->x,
			   cy = q->y-c->y,
			  dcx = -cy, // rotate by pi/2
			  dcy =  cx;
		if (q->r > c->r){
			dcx /= 3;
			dcy /= 3;
			if (fabs(cx) < fabs(cy)){
				cx = 0;
				cy = 2*cy/9;
			}else if (fabs(cy) < fabs(cx)){
				cx = 2*cx/9;
				cy = 0;
			}else
				assert(false);
		}else{
			cx *= .5;
			cy *= .5;
			dcx *= .5;
			dcy *= .5;
		}
		px = cx-dcx;
		py = cy-dcy;
		qx = cx+dcx;
		qy = cy+dcy;
	}
	void sample_u(double sx, double sy, double& dx, double& dy)const{
		assert(!child[0]);
		sx -= x;
		sy -= y;
		dx = dy = nan("not found");
		int i = 0;
		const bool failed = (const_cast<quad*>(this))->visit_my_vertices([sx, sy, &i, &dx, &dy](quad *c, quad *p, quad*, quad *q){
			double D, Dpl, Dql, Drl, px, py, qx, qy;
			face_endpoints(c, q, px, py, qx, qy);
			if (barycentric(px, py, qx, qy, sx, sy, D, Dpl, Dql, Drl)){
				assert(D || !(Dpl || Dql || Drl));
				dx = (Dpl*c->dx[i]+Dql*c->dx[(i+1)%8]+Drl*c->cell_dx)/(D ? D : 1);
				dy = (Dpl*c->dy[i]+Dql*c->dy[(i+1)%8]+Drl*c->cell_dy)/(D ? D : 1);
				assert(!isnan(dx));
				assert(!isnan(dy));
				return false;
			}
			++i;
			return true;
		});
		assert(!isnan(dx));
		assert(!isnan(dy));
		assert(!failed);
	}
	void split(std::function<void(quad*)>cb=NULL){ // XXX check shared edge condition
		for (int i = 0; i < 4; ++i)
			if (neighbour[i] && neighbour[i]->r > r)
				neighbour[i]->split(cb);
		for (int i = 0; i < 4; ++i)
			assert(!neighbour[i] || neighbour[i]->r == r);
		for (int i = 0; i < 4; ++i){
			assert(!child[i]);
			if (cb)
				cb(child[i] = new quad(this, i));
		}
	}
	void merge(){
		for (int i = 0; i < 4; ++i)
			assert(child[i]);
		assert(!"merge not implemented");
	}
	double nx(int i)const{ // not dividing by radius sum
		assert(neighbour[i]);
		assert(!neighbour[i]->child[0]);
		return (neighbour[i]->x-x)/(neighbour[i]->r > r ? 1.5 : 1);
	}
	double ny(int i)const{
		assert(neighbour[i]);
		assert(!neighbour[i]->child[0]);
		return (neighbour[i]->y-y)/(neighbour[i]->r > r ? 1.5 : 1);
	}
	double n(int i)const{
		return neighbour[i] ? hypot(nx(i), ny(i)) : 2*r;
	}
	double theta(int i)const{
		assert(neighbour[i]);
		return phi_theta(phi, neighbour[i]->phi);
	}
	bool visit_neighbours(std::function<bool(quad*, quad*, double, double&)>cb){
		assert(this);
		assert(cb);
		for (int i = 0; i < 4; ++i)
			if (!neighbour[i])
				continue;
			else if(!neighbour[i]->child[0]){
				double old = u[i];
				bool quit = !cb(this, neighbour[i], n(i), u[i]);
				neighbour[i]->u[(i+2)%4] -= u[i]-old;
				if (quit)
					return false;
			}else
				for (int j = 1; j < 3; ++j){
					quad *const n = neighbour[i]->child[(i+j)%4];
					double tmp = -n->u[(i+2)%4];
					bool quit = !cb(this, n, n->n((i+2)%4), tmp);
					u[i] -= -tmp-n->u[(i+2)%4];
					n->u[(i+2)%4] = -tmp;
					if (quit)
						return false;
				}
		return true;
	}
	double div(){
		double ret = 0;
		visit_neighbours([&ret](quad *p, quad *q, double n, double& u){
			const double theta = phi_theta(p->phi, q->phi);
			if (theta)
				ret += u*n*(1-phi_theta(p->solid_phi, q->solid_phi))/max(1e-2, theta);
			return true;
		});
		return ret;
	}
	double div_old(){
		double ret = 0;
		for (int i = 0; i < 4; ++i)
			if (!neighbour[i])
				continue;
			else if(!neighbour[i]->child[0])
				ret += n(i)*u[i];
			else
				ret -= neighbour[i]->child[(i+1)%4]->n((i+2)%4)*neighbour[i]->child[(i+1)%4]->u[(i+2)%4]
					+  neighbour[i]->child[(i+2)%4]->n((i+2)%4)*neighbour[i]->child[(i+2)%4]->u[(i+2)%4];
		return ret;
	}
	bool visit_cells(std::function<bool(quad*)>cb){
		assert(this);
		assert(cb);
		if (child[0]){
			for (int i = 0; i < 4; ++i)
				if (!child[i]->visit_cells(cb))
					return false;
		}else
			return cb(this);
		return true;
	}
	bool visit_faces(std::function<bool(quad*, int)>cb){
		assert(this);
		assert(cb);
		return visit_cells([cb](quad *p){
			for (int j = 0; j < 4; ++j){
				quad *const q = p->neighbour[j];
				if (q &&
					!q->child[0] &&
					(p->r < q->r || (p->r == q->r && p < q)) &&
					!cb(p, j))
					return false;
			}
			return true;
		});
	}
	bool visit_my_vertices(std::function<bool(quad*, quad*, quad*, quad*)>cb){
		assert(!child[0]);
		char ghost_buf[3*sizeof(quad)];
		memset(ghost_buf, 0, sizeof(ghost_buf));
		quad *ghost = (quad*)ghost_buf;
		const quad ghost_values = ghost_cell();
		int nr_ghost = 0;
		quad *p = NULL, *q = NULL, *pq = NULL;
		for (int i = 0; i < 9; ++i){
			const int j = i/2%4;
			quad *const n = neighbour[j];
			if (!(n && n->child[0]))
				++i;
#define IF_TRY_GHOST(q,n,k) \
if (!(n)->neighbour[(k)]){\
	const double macro_gr = (n)->r; \
	new(&ghost[nr_ghost]) quad((n)->x+2*macro_gr*cos[(k)], (n)->y+2*macro_gr*sin[(k)], macro_gr);\
	ghost[nr_ghost].copy_from(&ghost_values);\
	(q) = &ghost[nr_ghost];\
	nr_ghost = (nr_ghost+1)%(sizeof(ghost_buf)/sizeof(ghost[0]));\
}
			IF_TRY_GHOST(q, this, j)
			else
				q = n->child[0] ? n->child[(j+2-i%2)%4] : n; // (j+1+(1-i)%2)%4 in python
			assert(q);
			if (p){
				if (r != q->r)
					pq = NULL;
				if (!pq && r == p->r && (q->y-y)*cos[j]-(q->x-x)*sin[j] == q->r-r){
					quad *const pred = q->neighbour[(j+3)%4];
					if (pred && pred->r == q->r){
						pq = pred;
						if (pq->child[0])
							pq = pq->child[(j+1)%4];
						if (q->r != pq->r)
							pq = NULL;
					}else IF_TRY_GHOST(pq, q, (j+3)%4);
					assert(!pq || q->r == pq->r);
					assert(!pq || fabs(pq->x-x) == pq->r+r);
					assert(!pq || fabs(pq->y-y) == pq->r+r);
				}
				if (!cb(this, p, pq, q))
					return false;
			}
			p = q;
			pq = NULL;
			if ((p->y-y)*cos[j]-(p->x-x)*sin[j] == r-p->r){
				quad *const succ = p->neighbour[(j+1)%4];
				if (succ && succ->r == p->r){
					pq = succ;
					if (pq->child[0])
						pq = pq->child[(j+2)%4];
					if (p->r != pq->r)
						pq = NULL;
				}else IF_TRY_GHOST(pq, p, (j+1)%4);
				assert(!pq || p->r == pq->r);
				assert(!pq || fabs(pq->x-x) == pq->r+r);
				assert(!pq || fabs(pq->y-y) == pq->r+r);
			}
#undef IF_TRY_GHOST
		}
		return true;
	}
};
const int quad::cos[4] = {1, 0, -1, 0}, quad::sin[4] = {0, 1, 0, -1};

void reconstruct_velocity(quad *root){
	static std::map<std::pair<const quad*, const quad*>, double>unx, uny, nxny, nx2, ny2;
	unx.clear();
	uny.clear();
	nxny.clear();
	nx2.clear();
	ny2.clear();
	root->visit_faces([&](const quad *const p, const int j){
		quad *const q = p->neighbour[j];
		assert(!p->child[0] && !q->child[0]);
		const std::pair<const quad*, const quad*>pq = std::make_pair(p, q), qp = std::make_pair(q, p);
		const static double epsilon = .001; // XXX adjust
		double    u = p->u[j],
				 nx = p->nx(j),
				 ny = p->ny(j),
				  n = hypot(nx, ny),
				 n2 = nx*nx+ny*ny,
			  w_inv = .25*n2+epsilon*epsilon,
			divisor = w_inv*w_inv;
		nx /= n;
		ny /= n;
#define FLOW(a,q) do{\
(a)[pq] += (q);\
(a)[qp] += (q);\
}while(0)
		FLOW(unx, u*nx/divisor); // cell x cell => edge
		FLOW(uny, u*ny/divisor); // don't premultiply
		FLOW(nxny, nx*ny/divisor);
		FLOW(nx2, nx*nx/divisor);
		FLOW(ny2, ny*ny/divisor);
#undef FLOW
		return true;
	});
	root->visit_cells([&](quad *const n){
		double cell_dx_n = 0, cell_dy_n = 0;
		int i = 0;
		const bool success = n->visit_my_vertices([&](quad *p, const quad *q, const quad *r, const quad *s){
			double my_unx = 0, my_uny = 0, my_nxny = 0, my_nx2 = 0, my_ny2 = 0;
#define VISIT(p,q) do{\
const std::pair<const quad*, const quad*>pq = std::make_pair(p, q);\
my_unx += unx[pq];\
my_uny += uny[pq];\
my_nxny += nxny[pq];\
my_nx2 += nx2[pq];\
my_ny2 += ny2[pq];\
}while(0)
			VISIT(p, q); // edge x edge => vertex
			if (r){
				VISIT(q, r);
				VISIT(r, s);
			}else
				VISIT(q, s);
			VISIT(s, p);
#undef VISIT
			const double det = my_nxny*my_nxny-my_nx2*my_ny2;
			if (det){
				p->dx[i] = (my_nxny*my_uny-my_ny2*my_unx)/det; // vertex* => cell
				p->dy[i] = (my_nxny*my_unx-my_nx2*my_uny)/det;
			}else if(my_nx2 || my_nxny || my_ny2){ // always true
				double A, B, C; // (A, B, C) = (nx2, nxny, unx)+(nxny, ny2, uny)
				if (my_nxny < 0 && my_nx2){ // can cancel in A, so first-second
					A = my_nx2-my_nxny; // XXX check if this is the regularized solution
					B = my_nxny-my_ny2;
					C = my_unx-my_uny;
				}else{ // cannot cancel in A, so first+second
					A = my_nx2+my_nxny;
					B = my_nxny+my_ny2;
					C = my_unx+my_uny;
				}
				const double coef = C/(A*A+B*B);
				p->dx[i] = A*coef;
				p->dy[i] = B*coef;
			}else{
				assert(!(my_unx || my_uny)); // corner
				p->dx[i] = 0;
				p->dy[i] = 0;
			}
			cell_dx_n += p->dx[i  ];
			cell_dy_n += p->dy[i++];
			assert(i <= 8);
			if (i != 8){
				p->dx[i] = p->dx[0];
				p->dy[i] = p->dy[0];
			}
			return true;
		});
		assert(success);
		n->cell_dx = cell_dx_n/i;
		n->cell_dy = cell_dy_n/i;
		return true;
	});
}

void project(std::vector<quad*>& a){
	static std::map<quad*, size_t>row;
	static std::vector<double>rhs;
	static std::vector<double>result;
	row.clear();
	rhs.clear();
	result.clear();

	for (quad *n : a)
		if (!n->child[0])
			for (int j = 0; j < 4; ++j)
				if (!n->neighbour[j])
					n->u[j] = 0;

	for (quad *n : a)
		if (!n->child[0] && n->phi < 0){
			row[n] = rhs.size();
			rhs.push_back(n->div());
		}
	SparseMatrix<double>mat(rhs.size());
	result.resize(rhs.size());
	for (quad *n : a)
		if (!n->child[0] && n->phi < 0)
			n->visit_neighbours([&mat](quad *p, quad *q, double n, double& u){
				const double theta = phi_theta(p->phi, q->phi);
				if (!theta)
					return true;
				double w = n*(1-phi_theta(p->solid_phi, q->solid_phi));
				const int pr = row[p];
				if (q->phi < 0)
					mat.add_to_element(pr, row[q], -w);
				else
					w /= max(1e-2, theta);
				mat.add_to_element(pr, pr, w);
				return true;
			});
	rpc("check_symmetric", mat);
	double residual;
	int iterations;
	PCGSolver<double>solver;
#ifdef NDEBUG
	solver.set_solver_parameters(1e-5, 1000, .97, .25);
#endif
	//std::cerr << "mat = " << mat << " rhs = " << rhs << std::endl;
	assert(rhs == rhs);
	bool success = !rhs.size() || solver.solve(mat, rhs, result, residual, iterations);
	std::cerr << "residual = " << residual << " iterations = " << iterations << " success = " << success << std::endl;
	for (quad *n : a)
		if (!n->child[0] && n->phi < 0)
			n->visit_neighbours([](quad *p, quad *q, double n, double& u){
				double pp = row.count(p) ? result[row[p]] : 0,
					   qp = row.count(q) ? result[row[q]] : 0;
				if (phi_theta(p->phi, q->phi) && p < q) // pointer comparison
					u += qp-pp;
				return true;
			});
#ifndef NDEBUG
	for (quad *n : a)
		if (!n->child[0] && n->phi < 0)
			assert(fabs(n->div()) <= 1e-2);
#endif
}

void advect_velocity(quad *root){
	static std::vector<double>nu;
	nu.clear();
	root->visit_faces([root](quad *const p, int j){
		bool is_slanted_face = p->neighbour[j]->r > p->r;
		double dx, dy, nx = p->nx(j), ny = p->ny(j),
			cx = p->x+p->r*quad::cos[j]*(1-is_slanted_face/3.),
			cy = p->y+p->r*quad::sin[j]*(1-is_slanted_face/3.);
		p->sample_u(cx, cy, dx, dy);
		root->query_sample_u(
			max(root->x-root->r, min(root->x+root->r, cx-dx)),
			max(root->y-root->r, min(root->y+root->r, cy-dy)), dx, dy);
		nu.push_back((dx*nx+dy*ny)/hypot(nx, ny));
		return true;
	});

	static int i;
	i = 0;
	root->visit_faces([](quad *const p, int j){
		p->neighbour[j]->u[(j+2)%4] = -(p->u[j] = nu[i++]);
		return true;
	});
}

void advect_phi(quad *root, std::vector<quad*>& a){
	static std::vector<double>phi;
	phi.clear();
	phi.resize(a.size());
	static_assert(std::is_standard_layout<quad>::value, "cannot use offsetof");
	for (int i = 0; i < a.size(); ++i)
		if (!a[i]->child[0]){
			double dx, dy;
			root->query_sample_u(a[i]->x, a[i]->y, dx, dy);
			phi[i] = root->query_sample(a[i]->x-dx, a[i]->y-dy, offsetof(quad, phi));
		}
	for (int i = 0; i < a.size(); ++i)
		if (!a[i]->child[0])
			a[i]->phi = phi[i];
}

void advect(quad *root, std::vector<quad*>& a){
	advect_velocity(root);
	advect_phi(root, a);
}
void interpolate_surface(quad *root, std::vector<double>& bx, std::vector<double>& by){
	bx.clear();
	by.clear();
	root->visit_faces([&bx, &by](quad *p, int j){
		quad *q = p->neighbour[j];
		if ((p->phi < 0) == (q->phi < 0))
			return true;
		if (q->phi < 0)
			swap(p, q);
		const double theta = phi_theta(p->phi, q->phi);
		bx.push_back((1-theta)*p->x+theta*q->x);
		by.push_back((1-theta)*p->y+theta*q->y);
		return true;
	});
}

template<>
void _print_array_contents<quad*>(std::ostream& os, quad *const& elt){
	os << "{\"phi\":" << elt->phi << ",\"solid_phi\":" << elt->solid_phi << ",\"x\":" << elt->x << ",\"y\":" << elt->y << ",\"r\":" << elt->r << ",\"leaf\":" << !elt->child[0] << "}";
}

int main(){
	const double gx = 0, gy = -.05, T = 50;
	quad *root = new quad(0, 0, 1);
	static std::vector<double>bx, by;
	std::vector<quad*>a;

	a.push_back(root);
	for (int i = 0; i < a.size(); ++i){
		quad *const c = a[i];
		//std::cerr << c->x << ", " << c->y << ", " << c->r << std::endl;
		//c->solid_phi = 1-hypot(c->x, c->y);
		//c->phi = max(-c->solid_phi, c->x+.25*c->y);
		c->solid_phi = 1;
		c->phi = c->y;
		for (int i = 0; i < 4; ++i)
			c->u[i] = 0;
		for (int i = 0; i < 8; ++i)
			c->dx[i] = c->dy[i] = 0;
		c->cell_dx = c->cell_dy = 0;
		if (c->r < (min(fabs(c->solid_phi), fabs(c->phi)) < 1.42e-1 ? 1e-2 : 1e-1))
			continue;
		c->split([&a](quad *n){
			a.push_back(n);
		});
	}
	root->check_relations();
	interpolate_surface(root, bx, by);
	rpc("draw_quad", a, bx, by);

	for (int t = 0; t < T; ++t){
		// forces
		root->visit_faces([gx, gy](quad *const p, int j){
			quad *const q = p->neighbour[j];
			if (p->phi < 0 || q->phi < 0){
				//const double flow = (gx*p->nx(j)+gy*p->ny(j))*p->theta(j);
				const double flow = (gx*p->nx(j)+gy*p->ny(j))*p->theta(j)/p->n(j);
				p->u[j] += flow;
				q->u[(j+2)%4] -= flow;
			}
			return true;
		});

		project(a);

		reconstruct_velocity(root);

		advect(root, a);

		a.clear();
		a.push_back(root);
		for (int i = 0; i < a.size(); ++i)
			if (a[i]->child[0])
				for (int j = 0; j < 4; ++j)
					a.push_back(a[i]->child[j]);
		root->check_relations();
		interpolate_surface(root, bx, by);
		rpc("draw_quad", a, bx, by);
	}
	return 0;
}

// old code

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

#if 0
int main(){
	for (int t = 0; t < T; ++t){
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
			interpolate_surface(phi, bx, by);
			if ((t+1)%redistance_period == 0)
				redistance(phi, bx, by);
			rpc("draw", solid_phi, dx, dy, phi, bx, by);
		}
	}
	return 0;
} // vim: set ts=4 sw=4 noet:
#endif

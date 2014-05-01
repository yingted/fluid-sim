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

double solid_phi(double x, double y);

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
	double u[4], phi, dx[8], dy[8], cell_dx, cell_dy;
	void copy_from(const quad *o){
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
#ifdef NDEBUG
		return;
#endif
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
		return;
		sx += x;
		sy += y;
		const double sp = solid_phi(sx, sy);
		if (sp < 0){
			const static double eps = pow(2, -((53+1)/3));
			const double eps_x = max(1e-6, fabs(sx))*eps,
			             eps_y = max(1e-6, fabs(sy))*eps,
			                rx = (solid_phi(sx+eps_x, sy)-sp)/eps_x,
			                ry = (solid_phi(sx, sy+eps_y)-sp)/eps_y,
			                 r = (dx*rx+dy*ry)/(rx*rx+ry*ry);
			dx -= r*rx;
			dy -= r*ry;
		}
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
			if (phi_theta(p->phi, q->phi)){
				double px, py, qx, qy;
				face_endpoints(p, q, px, py, qx, qy);
				assert(n == hypot(qx-px, qy-py));
				ret += u*n*(1-phi_theta(solid_phi(px, py), solid_phi(qx, qy)));
			}
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

	SparseMatrix<double>mat(0);
	for (quad *n : a)
		if (!n->child[0])
			n->visit_neighbours([&mat](quad *p, quad *q, double n, double& u){
				double theta = phi_theta(p->phi, q->phi);
				double px, py, qx, qy;
				quad::face_endpoints(p, q, px, py, qx, qy);
				double solid = phi_theta(solid_phi(px, py), solid_phi(qx, qy)), w = 1-solid;
				if (!theta || !w){
					u = 0;
					return true;
				}

				int pr, qr;
				if (row.count(p))
					pr = row[p];
				else{
					pr = rhs.size();
					row[p] = pr;
					rhs.push_back(p->div());
				}
				if (row.count(q))
					qr = row[q];
				else{
					qr = row.size();
					row[q] = qr;
					rhs.push_back(q->div());
				}
				assert(row.size() == rhs.size());
				mat.resize(rhs.size());

				if (solid) // partial solid face
					theta = 1; // activate cells
				w *= n/hypot(q->x-p->x, q->y-p->y);
				if (theta == 1) // both active
					mat.add_to_element(pr, qr, -w);
				else if (!(p->phi < 0)) // one inactive, p
					return true;
				else // one inactive, q
					w /= max(1e-6, theta);
				mat.add_to_element(pr, pr, w);
				return true;
			});
#ifndef NDEBUG
	for (std::pair<quad *, size_t>e : row)
		assert(rhs[e.second] == e.first->div());
#endif

	rpc("check_symmetric", mat);
	double residual;
	int iterations;
	PCGSolver<double>solver;
#ifdef NDEBUG
	solver.set_solver_parameters(1e-5, 1000, .97, .25);
#endif
	//std::cerr << "mat = " << mat << " rhs = " << rhs << std::endl;
	assert(rhs == rhs);
	result.resize(rhs.size());
	bool success = !rhs.size() || solver.solve(mat, rhs, result, residual, iterations);
	std::cerr << "cells = " << rhs.size();
	if (rhs.size())
		std::cerr << " residual = " << residual << " iterations = " << iterations << " success = " << success;
	std::cerr << std::endl;
	for (quad *n : a)
		if (!n->child[0])
			n->visit_neighbours([](quad *p, quad *q, double n, double& u){
				double pp = row.count(p) ? result[row[p]] : 0,
				       qp = row.count(q) ? result[row[q]] : 0,
				    theta = phi_theta(p->phi, q->phi);
				double px, py, qx, qy;
				quad::face_endpoints(p, q, px, py, qx, qy);
				double solid = phi_theta(solid_phi(px, py), solid_phi(qx, qy)), w = 1-solid;
				if (!theta || !w)
					return true;
				if (solid) // partial solid face
					theta = 1; // activate cells
				if (p < q) // pointer comparison
					u += (qp-pp)/(max(1e-6, theta)*hypot(q->x-p->x, q->y-p->y));
				return true;
			});
#ifndef NDEBUG
	double worst = 0;
	for (quad *n : a)
		if (!n->child[0] && n->phi < 0)
			worst = max(worst, fabs(n->div()));
			//assert(fabs(n->div()) <= 1e-6);
	std::cerr << "worst = " << worst << std::endl;
#endif
}

void advect_velocity(quad *root){
	static std::vector<double>nu;
	nu.clear();
	root->visit_faces([root](quad *const p, int j){
		bool is_slanted_face = p->neighbour[j]->r > p->r;
		double dx, dy, dx2, dy2, nx = p->nx(j), ny = p->ny(j),
			cx = p->x+p->r*quad::cos[j]*(1-is_slanted_face/3.),
			cy = p->y+p->r*quad::sin[j]*(1-is_slanted_face/3.);
		p->sample_u(cx, cy, dx, dy);
		root->query_sample_u(
			max(root->x-root->r, min(root->x+root->r, cx-dx)),
			max(root->y-root->r, min(root->y+root->r, cy-dy)), dx2, dy2);
		nu.push_back(.5*((dx+dx2)*nx+(dy+dy2)*ny)/hypot(nx, ny));
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
	typedef std::pair<double, quad*>cell_t;
	static std::priority_queue<cell_t, std::vector<cell_t>, std::greater<cell_t> >q;
	static std::map<quad*, int>ancestor;
	static std::map<quad*, double>dist;
	ancestor.clear();
	dist.clear();
	bx.clear();
	by.clear();
#define NEIGH(p,i) do{\
	const double d2 = hypot(p->x-bx[i], p->y-by[i]);\
	if (!dist.count(p) || d2 < dist[p]){\
		dist[p] = d2;\
		ancestor[p] = i;\
		q.push(std::make_pair(d2, p));\
	}\
}while(0)

	root->visit_faces([&bx, &by](quad *u, int j){
		quad *v = u->neighbour[j];
		if ((u->phi < 0) == (v->phi < 0))
			return true;
		if (v->phi < 0)
			swap(u, v);
		const double theta = phi_theta(u->phi, v->phi);
		bx.push_back((1-theta)*u->x+theta*v->x);
		by.push_back((1-theta)*u->y+theta*v->y);
		NEIGH(u, bx.size()-1);
		NEIGH(v, bx.size()-1);
		return true;
	});
	while(!q.empty()){ // redistance
		const double d = q.top().first;
		quad *u = q.top().second;
		q.pop();
		if (d != dist[u])
			continue;
		const int anc = ancestor[u];
		double dx_n = 0, dy_n = 0, valid_area = 0;
#define SEEN (dist.count(v) && (dist[v] < d || (dist[v] == d && v < u)))
		u->visit_neighbours([&, d](quad *u, quad *v, double area, double& flow){
			if (!(u->phi < 0) && (SEEN || v->phi < 0)){ // air and valid
				const double nx = v->x-u->x, ny = v->y-u->y, n = hypot(nx, ny);
				valid_area += area;
//std::cerr << "(" << nx/n*(area*flow) << ", " << ny/n*(area*flow) << ") ";
				dx_n += nx/n*(area*flow);
				dy_n += ny/n*(area*flow);
			}
			return true;
		});
//if (valid_area)
//	std::cerr << "=> (" << dx_n << ", " << dy_n << ")/" << valid_area << " =>";
		u->visit_neighbours([=](quad *u, quad *v, double area, double& flow){
			const bool seen = SEEN;
			if (!(u->phi < 0) && !(seen || v->phi < 0)){ // air and invalid
				const double nx = v->x-u->x, ny = v->y-u->y, n = hypot(nx, ny);
//std::cerr << " " << (nx*dx_n+ny*dy_n)/(n*valid_area);
				flow = (nx*dx_n+ny*dy_n)/(n*valid_area);
			}
			if (!seen)
				NEIGH(v, anc);
			return true;
		});
//if (valid_area)
//	std::cerr << std::endl;
#undef SEEN
	}
	root->visit_cells([](quad *n){
		n->phi = copysign(dist[n], n->phi);
		return true;
	});
}

void extrapolate_solid(std::vector<quad*>& a){
	return;
	for (quad *p : a)
		if (!p->child[0])
			p->visit_neighbours([](quad *p, quad *q, double n, double& u){
				double px, py, qx, qy;
				quad::face_endpoints(p, q, px, py, qx, qy);
				u *= fabs(solid_phi(qx, qy)-solid_phi(px, py))/hypot(qx-px, qy-py);
				return true;
			});
}

template<>
void _print_array_contents<quad*>(std::ostream& os, quad *const& elt){
	double dx, dy;
	elt->query_sample_u(elt->x, elt->y, dx, dy);
	os << "{\"phi\":" << elt->phi << ",\"solid_phi\":" << solid_phi(elt->x, elt->y) << ",\"x\":" << elt->x << ",\"y\":" << elt->y << ",\"r\":" << elt->r << ",\"leaf\":" << !elt->child[0] << ",\"dx\":" << dx << ",\"dy\":" << dy << "}";
}

double solid_phi(double x, double y){
	return .9-hypot(x, y); // XXX put in quadtree
	//return y+2;
}

int main(){
	const double gx = 0, gy = -.05, T = 5;
	quad *root = new quad(0, 0, 1);
	static std::vector<double>bx, by;
	std::vector<quad*>a;

	a.push_back(root);
	for (int i = 0; i < a.size(); ++i){
		quad *const c = a[i];
		//std::cerr << c->x << ", " << c->y << ", " << c->r << std::endl;
		c->phi = max(-solid_phi(c->x, c->y), c->x+.25*c->y);
		//c->phi = c->x+.25*c->y;
		//c->phi = c->y;
		for (int i = 0; i < 4; ++i)
			c->u[i] = 0;
		for (int i = 0; i < 8; ++i)
			c->dx[i] = c->dy[i] = 0;
		c->cell_dx = c->cell_dy = 0;
		//if (c->r < (min(fabs(solid_phi(c->x, c->y)), fabs(c->phi)) < 1.42e-1 ? 1e-2 : 1e-1))
		if (c->r < 3e-2)
			continue;
		c->split([&a](quad *n){
			a.push_back(n);
		});
	}
	root->check_relations();
	interpolate_surface(root, bx, by);
	rpc("draw_quad", a, bx, by);

	for (int t = 0; t < T; ++t){
		root->visit_faces([gx, gy](quad *const p, int j){ // forces
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
		interpolate_surface(root, bx, by);
		extrapolate_solid(a);
		reconstruct_velocity(root);
		advect(root, a);
		root->check_relations();
		a.clear();
		a.push_back(root);
		for (int i = 0; i < a.size(); ++i)
			if (a[i]->child[0])
				for (int j = 0; j < 4; ++j)
					a.push_back(a[i]->child[j]);
		rpc("draw_quad", a, bx, by);
	}
	return 0;
} // vim: ts=4 sw=4

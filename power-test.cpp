// taken from CGAL-4.4:examples/Triangulation_3/regular_3.cpp
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/centroid.h>
#include <CGAL/barycenter.h>
#include <CGAL/bounding_box.h>
#include <CGAL/enum.h>
#include <CGAL/Origin.h>
#include <boost/random/variate_generator.hpp>
#include <boost/nondet_random.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <fstream>
#include <cassert>
#include <vector>
#include <queue>
#include <set>

#ifdef EXACT
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel   K;
#else
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
#endif

typedef K::Vector_3                                         Vector_3;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;

typedef Traits::RT                                          Weight;
typedef Traits::Bare_point                                  Point;
typedef Traits::Weighted_point                              Weighted_point;
typedef Traits::Iso_cuboid_3                                Iso_cuboid_3;

typedef CGAL::Regular_triangulation_3<Traits>               Rt;

typedef Rt::Finite_edges_iterator                           Finite_edges_iterator;
typedef Rt::Cell_circulator                                 Cell_circulator;
typedef Rt::Cell_handle                                     Cell_handle;

class rand_bool{
public:
	virtual bool operator()() = 0;
	virtual void report(){};
	virtual void feed(bool (*cb)(rand_bool&)){
		while (cb(*this) && next());
	}
protected:
	virtual bool next() = 0; // cannot be called after returning false
};

class exhaustive_rng : public rand_bool{
private:
	std::vector<bool>history;
	int pos;
public:
	virtual bool operator()(){
		if (pos == history.size())
			history.push_back(false);
		return history[pos++];
	}
	virtual void feed(bool (*cb)(rand_bool&)){
		history.clear();
		pos = 0;
		rand_bool::feed(cb);
	}
	virtual void report(){
		std::cerr << "rng: ";
		for (int i = 0; i < pos; ++i)
			std::cerr << history[i];
		std::cerr << std::endl;
	}
protected:
	virtual bool next(){
		assert (pos == history.size());
		while (!history.empty() && history.back())
			history.pop_back();
		if (history.empty())
			return false;
		history.pop_back();
		history.push_back(true);
		pos = 0;
		return true;
	}
};

class random_rng : public rand_bool{
private:
	int count, stop;
	const boost::bernoulli_distribution<>dist;
	boost::random_device dev;
	boost::variate_generator<boost::random_device&, boost::bernoulli_distribution<> >rng;
public:
	random_rng(int n) : count(0), stop(n), dev(), dist(), rng(dev, dist){}
	virtual bool operator()(){
		return rng();
	}
protected:
	virtual bool next(){
		return ++count < stop;
	}
};

const Iso_cuboid_3 domain(-1, -1, -1, 1, 1, 1);
double branch_radius_cutoff;
std::map<const char*, int>hit_counter;
std::map<K::FT, int>area_counter;

void rand_pts(std::set<Weighted_point>& pts, Point c, double r, rand_bool& rng){
	std::queue<std::pair<Point, double> >q;
	q.push(std::make_pair(c, r));
	while (!q.empty()){
		c = q.front().first;
		r = q.front().second;
		q.pop();
		if (r <= branch_radius_cutoff)
			goto leaf;
		{
			const Weight parent_weight(3*((2*r)*(2*r)));
			for (int i = -3; i <= 3; i += 2)
				for (int j = -3; j <= 3; j += 2)
					for (int k = -3; k <= 3; k += 2)
						if (i*i+j*j+k*k == 11 || i*i+j*j+k*k == 19){
							const Point p = c+r*Vector_3(i, j, k);
							if (!domain.has_on_unbounded_side(p) && pts.count(Weighted_point(p, parent_weight)))
								goto leaf;
						}
		}
		if (!rng()){
leaf:
			const Weight w(3*(r*r));
			pts.insert(Weighted_point(c, w));
			for (int i = -1; i <= 1; ++i) // insert ghost cells to check boundary conditions
				for (int j = -1; j <= 1; ++j)
					for (int k = -1; k <= 1; ++k)
						if (i*i+j*j+k*k == 1){
							const Point p = c+(2*r)*Vector_3(i, j, k);
							if (domain.has_on_unbounded_side(p))
								pts.insert(Weighted_point(p, w)); // ghost cell outside
						}
			continue;
		}
		const double child_r = .5*r;
		for (int i = -1; i <= 1; i += 2)
			for (int j = -1; j <= 1; j += 2)
				for (int k = -1; k <= 1; k += 2)
					q.push(std::make_pair(c+child_r*Vector_3(i, j, k), child_r));
	}
}

time_t start, last;
#ifdef OUTPUT
std::ofstream obj(OUTPUT);
#endif

int check_uniqueness(const Point p) { // 0 <= x <= y <= z
	auto c_0_x = CGAL::compare(0, p.x());
	auto c_0_y = CGAL::compare(0, p.y());
	auto c_0_z = CGAL::compare(0, p.z());
	auto c_x_y = CGAL::compare(p.x(), p.y());
	auto c_y_z = CGAL::compare(p.y(), p.z());
	if (c_0_x == CGAL::EQUAL || c_0_y == CGAL::EQUAL || c_0_z == CGAL::EQUAL) // x=0, y=0, z=0
		return 0;
	if (c_x_y == CGAL::EQUAL || c_y_z == CGAL::EQUAL || p.z() == p.x()) // x=y, y=z, z=x
		return 0;
	if (c_0_x == CGAL::LARGER)
		return -1;
	if (c_x_y == CGAL::LARGER)
		return -1;
	if (c_y_z == CGAL::LARGER)
		return -1;
	return 1;
}

const static Point ORIGIN(CGAL::ORIGIN);

struct symmetric_cmp {
	bool operator()(const Weighted_point a, const Weighted_point b) const {
		switch (CGAL::compare(a.weight(), b.weight())) {
			case CGAL::SMALLER:
				return false;
			case CGAL::LARGER:
				return true;
		}
		switch (CGAL::compare_distance_to_point(ORIGIN, a.point(), b.point())) {
			case CGAL::SMALLER:
				return true;
			case CGAL::LARGER:
				return false;
		}
		return false;
	}
};

bool test_octree(rand_bool& rng){
	std::set<Weighted_point>pts;
	rand_pts(pts, ORIGIN, 1, rng);
	//std::cout << "cells (including ghost): " << pts.size() << std::endl;
	//for (std::set<Weighted_point>::const_iterator it = pts.begin(); it != pts.end(); ++it)
	//	std::cout << *it << std::endl;
	{
		std::multiset<Weighted_point, symmetric_cmp> pts_m(pts.begin(), pts.end());
		for (auto it = pts_m.begin(); it != pts_m.end();) {
			const auto it_e = pts_m.upper_bound(*it);
			switch (check_uniqueness(CGAL::centroid(it, it_e, CGAL::Dimension_tag<0>()))) {
				case -1:
					++hit_counter["pruned centroid"];
					goto end;
				case 1:
					goto triangulate;
			}
			it = it_e;
		}
	}
triangulate:
	{
		Rt T;
		ptrdiff_t num_inserted = T.insert(pts.begin(), pts.end());
		assert(pts.size() == num_inserted); // check for hidden XXX need to change traits
		assert(T.is_valid());
		//std::cout << "faces: " << T.number_of_finite_edges() << std::endl;
		for (Finite_edges_iterator it = T.finite_edges_begin(); it != T.finite_edges_end(); ++it){
			Weighted_point a = it->first->vertex(it->second)->point(),
						   b = it->first->vertex(it->third)->point();
			switch (domain.has_on_unbounded_side(a)+domain.has_on_unbounded_side(b)){
				case 0:
					++hit_counter["interior faces"];
					break;
				case 1:
					++hit_counter["border faces"]; // XXX need to check this
#ifdef BORDER
					break;
#endif
					continue;
				case 2:
					++hit_counter["ghost faces"];
					continue;
			}

			Vector_3 area2(0, 0, 0);
			Cell_circulator u = T.incident_cells(*it), end = u; // loop through vertices
			assert(u != 0);
#ifdef OUTPUT
			static int vertices = 0;
			int old_vertices = vertices;
			int count = hit_counter["octrees"], wrap = 10;
			Point origin = CGAL::ORIGIN+(5*Vector_3(count/wrap/wrap, count/wrap%wrap, count%wrap));
#endif
			Vector_3 u_vec = T.dual(u)-CGAL::ORIGIN;
			do{
#ifdef OUTPUT
				const Point o = origin+u_vec;
#ifdef EXACT_OUTPUT
				const CGAL::Gmpq x = o.x().exact(), y = o.y().exact(), z = o.z().exact();
				CGAL::Gmpz w0 = integral_division(x.denominator(), gcd(x.denominator(), y.denominator()))*y.denominator(),
							w = integral_division(w0, gcd(w0, z.denominator()))*z.denominator();
				obj << "v " << (x*w).numerator() << ' ' << (y*w).numerator() << ' ' << (z*w).numerator() << ' ' << w << '\n';
#else
				obj << "v " << o << '\n';
#endif
				++vertices;
#endif
				++u;
				Vector_3 v_vec = T.dual(u)-CGAL::ORIGIN;
				area2 = area2+cross_product(u_vec, v_vec);
				u_vec = v_vec;
			}while(u != end);
#ifdef OUTPUT
			obj << 'f';
			for (int i = old_vertices+1; i <= vertices; ++i)
				obj << ' ' << i;
			obj << '\n';
#endif
			Weight w = CGAL::min(a.weight(), b.weight()); // normalize by smaller weight
			++area_counter[area2.squared_length()/(16*w*w/9)];
		}
	}
end:
	{
		const int tested_count = ++hit_counter["octrees"];
		const static int period = 5;
		time_t now = time(NULL);
		if (now-last >= period){
			last = now;
			std::cerr << tested_count << " octrees in " << now-start << " s" << std::endl;
			rng.report();
		}
	}
	return true;
}


int main(int argc, char *argv[]){
	if (argc != 2 && argc != 3){
		std::cerr << "usage: " << argv[0] << " depth [random runs]" << std::endl;
		return 1;
	}
	branch_radius_cutoff = pow(.5, atof(argv[1]));
	std::cout << "testing with branch node cutoff radius " << branch_radius_cutoff << std::endl;
	start = last = time(NULL);
	if (argc == 2)
		exhaustive_rng().feed(test_octree);
	else
		random_rng(atoi(argv[2])).feed(test_octree);
	std::cout << "finished in " << time(NULL)-start << " s" << std::endl;
	std::cout << "hits:" << std::endl;
	for (std::map<const char*, int>::const_iterator it = hit_counter.begin(); it != hit_counter.end(); ++it)
		printf("%10d %s\n", it->second, it->first);
	std::cout << "(area/r^2)^2, r <= R: " << area_counter.size() << " areas" << std::endl;
	for (std::map<Weight, int>::const_iterator it = area_counter.begin(); it != area_counter.end(); ++it){
		printf("%10d ", it->second);
		std::cout << it->first << std::endl;
	}
	return 0;
}

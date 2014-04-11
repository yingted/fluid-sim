// taken from CGAL-4.4:examples/Triangulation_3/regular_3.cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <cassert>
#include <vector>
#include <sstream>
#include <queue>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Vector_3                                         Vector_3;
typedef K::Iso_cuboid_3                                     Iso_cuboid_3;

typedef CGAL::Regular_triangulation_euclidean_traits_3<K>   Traits;

typedef Traits::RT                                          Weight;
typedef Traits::Bare_point                                  Point;
typedef Traits::Weighted_point                              Weighted_point;

typedef CGAL::Regular_triangulation_3<Traits>               Rt;

typedef Rt::Finite_edges_iterator                           Finite_edges_iterator;
typedef Rt::Vertex_iterator                                 Vertex_iterator;
typedef Rt::Vertex_handle                                   Vertex_handle;
typedef Rt::Facet                                           Facet;
typedef Rt::Segment                                         Segment;
typedef Rt::Cell_handle                                     Cell_handle;

struct rand_bool{
	std::vector<bool>history;
	int pos;
	bool operator()(){
		if (pos == history.size())
			history.push_back(false);
		return history[pos++];
	}
	bool next(){
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
	void feed(bool (*cb)(rand_bool&)){
		history.clear();
		pos = 0;
		while(cb(*this) && this->next());
	}
};

const Iso_cuboid_3 domain(-1, -1, -1, 1, 1, 1);
double branch_radius_cutoff;
std::map<const char*, int>hit_counter;

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
			const Weight parent_weight((2*r)*(2*r));
			for (int i = -3; i <= 3; i += 2)
				for (int j = -3; j <= 3; j += 2)
					for (int k = -3; k <= 3; k += 2)
						if (i*i+j*j+k*k == 11){
							const Point p = c+r*Vector_3(i, j, k);
							if (!domain.has_on_unbounded_side(p) && pts.count(Weighted_point(p, parent_weight)))
								goto leaf;
						}
		}
		if (!rng()){
leaf:
			const Weight w(r*r);
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

bool test_octree(rand_bool& rng){
	std::set<Weighted_point>pts;
	rand_pts(pts, Point(0, 0, 0), 1, rng);
	//std::cout << "cells (including ghost): " << pts.size() << std::endl;
	//for (std::set<Weighted_point>::const_iterator it = pts.begin(); it != pts.end(); ++it)
	//	std::cout << *it << std::endl;
	{
		std::map<Weight, int>weight_count;
		for (std::set<Weighted_point>::const_iterator it = pts.begin(); it != pts.end(); ++it)
			++weight_count[it->weight()];
		for (std::map<Weight, int>::const_iterator it = weight_count.begin(); it != weight_count.end(); ++it)
			std::cout << it->first << ": " << it->second << std::endl;
		std::cout << "total: " << pts.size() << std::endl;
		std::cout << std::endl;
	}
	Rt T;
	ptrdiff_t num_inserted = T.insert(pts.begin(), pts.end());
	assert(pts.size() == num_inserted); // check for hidden
	assert(T.is_valid());
	//std::cout << "faces: " << T.number_of_finite_edges() << std::endl;
	for (Finite_edges_iterator it = T.finite_edges_begin(); it != T.finite_edges_end(); ++it){
		Weighted_point a = it->first->vertex(it->second)->point(),
		               b = it->first->vertex(it->third)->point();
		if (domain.has_on_unbounded_side(a) && domain.has_on_unbounded_side(b))
			continue;
		if (a.weight() < b.weight())
			std::swap(a, b); // make a larger
		const Vector_3 d = b-a;
		const double R2 = a.weight(), r2 = b.weight(), d2 = d.squared_length();
#define COUNT(w) ((d.x()*d.x() == (w))+(d.y()*d.y() == (w))+(d.z()*d.z() == (w)))
		//std::cout << a << ", " << b << std::endl;
		const double r = sqrt(r2), R = sqrt(R2), rR = r+R;
		if (d2 == 3*rR*rR) // single point
			continue;
		switch (COUNT(rR*rR)){
			case 0: // 3d intersection volume
				assert(false);
			case 1: // 2d plane
				if (R2 == r2){
					++hit_counter["case 1: same level"];
					assert(COUNT(0) == 2);
				}else if (R2 == 4*r2){
					++hit_counter["case 2: unit \"knight's move\""];
					assert(COUNT(r2) == 2);
				}else
					assert(false);
				break;
			case 2: // 1d (small slice)
				if (R2 == 4*r2){
					++hit_counter["case 3: abbb"];
					assert(COUNT(r2) == 1);
					if (d[0]*d[0] == rR*rR)
						assert(pts.count(Weighted_point(a+Vector_3(d[0]-copysign(R, d[0]), d[1], d[2]), r2)));
					if (d[1]*d[1] == rR*rR)
						assert(pts.count(Weighted_point(a+Vector_3(d[0], d[1]-copysign(R, d[1]), d[2]), r2)));
					if (d[2]*d[2] == rR*rR)
						assert(pts.count(Weighted_point(a+Vector_3(d[0], d[1], d[2]-copysign(R, d[2])), r2)));
				}else if (R2 == r2){
					++hit_counter["case 4: axab"];
					assert(COUNT(0) == 1);
					//std::cout << a << ", " << b << std::endl;
					int adjacent_large =
						+ (d[0] && pts.count(Weighted_point(a+Vector_3(d[0]-copysign(R, d[0]), d[1], d[2]), R2)))
						+ (d[1] && pts.count(Weighted_point(a+Vector_3(d[0], d[1]-copysign(R, d[1]), d[2]), R2)))
						+ (d[2] && pts.count(Weighted_point(a+Vector_3(d[0], d[1], d[2]-copysign(R, d[2])), R2))),
						adjacent_small =
						+ (d[0] && pts.count(Weighted_point(a+Vector_3(d[0]-copysign(.5*R, d[0]), d[1], d[2]), .25*R2)))
						+ (d[1] && pts.count(Weighted_point(a+Vector_3(d[0], d[1]-copysign(.5*R, d[1]), d[2]), .25*R2)))
						+ (d[2] && pts.count(Weighted_point(a+Vector_3(d[0], d[1], d[2]-copysign(.5*R, d[2])), .25*R2)));
					assert(adjacent_small);
					assert(adjacent_large+adjacent_large == 2);
				}else
					assert(false);
				break;
			case 3: // 0d (must be no area)
				assert(false);
				break;
		}
#undef COUNT
	}
	const int tested_count = ++hit_counter["octrees"];

	{	
		const static int period = 5;
		time_t now = time(NULL);
		if (now-last >= period){
			last = now;
			std::cerr << tested_count << " octrees in " << now-start << " s" << std::endl;
			std::cerr << "rng: ";
			for (int i = 0; i < rng.pos; ++i)
				std::cerr << rng.history[i];
			std::cerr << std::endl;
		}
	}
	return true;
}


int main(int argc, char *argv[]){
	if (argc != 2){
		std::cerr << "usage: " << argv[0] << " depth" << std::endl;
		return 1;
	}
	branch_radius_cutoff = pow(.5, atof(argv[1]));
	std::cout << "testing with branch node cutoff radius " << branch_radius_cutoff << std::endl;
	start = last = time(NULL);
	rand_bool().feed(test_octree);
	std::cout << "finished in " << time(NULL)-start << " s" << std::endl;
	for (std::map<const char*, int>::const_iterator it = hit_counter.begin(); it != hit_counter.end(); ++it)
		printf("%10d %s\n", it->second, it->first);
	return 0;
}

// taken from CGAL-4.4:examples/Triangulation_3/regular_3.cpp
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <cassert>
#include <vector>
#include <sstream>

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

void rand_pts(std::set<Weighted_point>& pts, const Point c, const double r, rand_bool& rng){
	if (r <= branch_radius_cutoff)
		goto leaf;
	{
		const Weight parent_weight = (2*r)*(2*r);
		for (int i = -3; i <= 3; i += 2)
			for (int j = -3; j <= 3; j += 2)
				for (int k = -3; k <= 3; k += 2)
					if (i*i+j*j+k*k == 11 && pts.count(Weighted_point(c+r*Vector_3(i, j, k), parent_weight)))
						goto leaf;
	}
	if (!rng()){
leaf:
		const Weight w = r*r;
		pts.insert(Weighted_point(c, w));
		for (int i = -1; i <= 1; ++i) // insert ghost cells to check boundary conditions
			for (int j = -1; j <= 1; ++j)
				for (int k = -1; k <= 1; ++k)
					if (i*i+j*j+k*k == 1){
						const Point p = c+(2*r)*Vector_3(i, j, k);
						if (domain.has_on_unbounded_side(p))
							pts.insert(Weighted_point(p, w)); // ghost cell outside
					}
		return;
	}
	const double child_r = .5*r;
	for (int i = -1; i <= 1; i += 2)
		for (int j = -1; j <= 1; j += 2)
			for (int k = -1; k <= 1; k += 2)
				rand_pts(pts, c+child_r*Vector_3(i, j, k), child_r, rng);
}

int tested_count;

bool test_octree(rand_bool& rng){
	std::set<Weighted_point>pts;
	rand_pts(pts, Point(0, 0, 0), 1, rng);
	//std::cout << "cells (including ghost): " << pts.size() << std::endl;
	//for (std::set<Weighted_point>::const_iterator it = pts.begin(); it != pts.end(); ++it)
	//	std::cout << *it << std::endl;
	Rt T;
	ptrdiff_t num_inserted = T.insert(pts.begin(), pts.end());
	assert(pts.size() == num_inserted); // check for hidden
	assert(T.is_valid());
	//std::cout << "faces: " << T.number_of_finite_edges() << std::endl;
	for (Finite_edges_iterator it = T.finite_edges_begin(); it != T.finite_edges_end(); ++it){
		const Weighted_point a = it->first->vertex(it->second)->point(),
		                     b = it->first->vertex(it->third)->point();
		if(domain.has_on_unbounded_side(a) && domain.has_on_unbounded_side(b))
			continue;
		const double R2 = CGAL::max(a.weight(), b.weight()), r2 = CGAL::min(a.weight(), b.weight());
		const Vector_3 d = b-a;
#define COUNT(w) ((d.x()*d.x() == (w))+(d.y()*d.y() == (w))+(d.z()*d.z() == (w)))
		//std::cout << a << ", " << b << std::endl;
		const double rR = sqrt(r2)+sqrt(R2);
		switch (COUNT(rR*rR)){
			case 0: // 3d intersection volume
				assert(false);
			case 1: // 2d plane
				break;
			default: // 1d, 0d (no area)
				continue;
		}
		if (R2 == r2) // same level
			assert(COUNT(0)+COUNT(4*r2) == 3); // COUNT(4*r2) == 1 for equal cell
		else if (R2 == 4*r2){
			const int far = COUNT(9*r2); // far == 1 for smaller cell
			assert(far);
			assert(COUNT(r2)+far == 3);
		}
	}
	++tested_count;
	//std::cout << '.';
	//std::cout.flush();
	return true;
}


int main(int argc, char *argv[]){
	if (argc != 2){
		std::cerr << "usage: " << argv[0] << " depth" << std::endl;
		return 1;
	}
	branch_radius_cutoff = pow(.5, atof(argv[1]));
	std::cout << "testing with branch node cutoff radius " << branch_radius_cutoff << std::endl;
	rand_bool().feed(test_octree);
	//std::cout << std::endl;
	std::cout << "tested " << tested_count << " octrees" << std::endl;
	return 0;
}

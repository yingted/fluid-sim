#include <iostream>
#include <vector>
#include <sstream>
#include <cassert>
#include "util.hpp"

int main(){
	const int N = 2, M = 2;
	grid a = make_grid(N, M);
	a[0][0] = 3.14;
	auto orig = a;
	std::string s = "[[3.14,0],[0,0]]";
	std::stringstream ss;
	ss << a;
	assert(ss.str() == s);
	ss >> a;
	assert(a == orig);
	return 0;
}

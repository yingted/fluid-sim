#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <vector>
#include <string>

template<typename E>
std::ostream& operator<<(std::ostream& os, const std::vector<E>& vec){
	char next = '[';
	for (const E& e : vec){
		os << next << e;
		next = ',';
	}
	if (next != ',') // eat trailing comma
		os << next;
	return os << "]";
}

std::istream& operator>>(std::istream& is, char ch){
	if (is.peek() == ch)
		is.ignore();
	else
		is.setstate(std::istream::failbit);
	return is;
}

template<typename E>
std::istream& operator>>(std::istream& is, std::vector<E>& vec){
	std::vector<E> v;
	bool fail = false;
	if (is >> '['){
		do{
			E e;
			fail |= !(is >> e);
			v.push_back(e);
		}while (is >> ',');
		if (is.rdstate() == std::istream::failbit)
			is.clear();
		is >> ']';
		if (fail){
			std::cout << "fail" << std::endl;
			is.setstate(std::istream::failbit);
		}
	}
	vec = v;
	return is;
}

typedef std::vector<std::vector<double> > grid;

extern inline grid make_grid(const size_t N, const size_t M, double val=0){
	return std::move(std::vector<std::vector<double> > (N, std::vector<double> (M, val)));
}

#endif

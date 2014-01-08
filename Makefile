util-test: util-test.cpp util.hpp
	g++ -std=c++11 -o util-test{,.cpp}
smoke-grid: smoke-grid.cpp util.hpp
	g++ -std=c++11 -o smoke-grid{,.cpp}
clean:
	rm -f smoke-grid-*.ppm util-test smoke-grid
smoke-grid-play: smoke-grid
	./smoke-grid
	ffplay -framerate 10 -vf scale=-1:500 smoke-grid-%04d.ppm
.PHONY: clean smoke-grid-play

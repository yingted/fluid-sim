FFMPEG_FLAGS := -framerate 10
CFLAGS := -std=c++11 -Ipcgsolver
util-test: util-test.cpp util.hpp
	g++ $(CFLAGS) -o util-test{,.cpp}
util-test-run: util-test
	bash -c "./util-test | diff <(echo -n '{\"method\":\"test\",\"params\":[\"string\"]}') -"
smoke-grid: smoke-grid.cpp util.hpp
	g++ $(CFLAGS) -o smoke-grid{,.cpp}
fluid-grid: fluid-grid.cpp util.hpp pcgsolver
	g++ $(CFLAGS) -o fluid-grid{,.cpp} -L/usr/lib64/atlas -lcblas
clean:
	rm -f {smoke,fluid}-grid{,-*.ppm,.ogg} util-test
smoke-grid-0001.ppm: smoke-grid
	./smoke-grid
smoke-grid-play: smoke-grid-0001.ppm
	ffplay $(FFMPEG_FLAGS) -vf scale=-1:500 smoke-grid-%04d.ppm
fluid-grid-0001.ppm: fluid-grid
	set -o pipefail; ./fluid-grid | ./rpc.py
fluid-grid-play: fluid-grid-0001.ppm
	ffplay $(FFMPEG_FLAGS) -vf scale=-1:500 fluid-grid-%04d.ppm
%.ogg: %-0001.ppm
	ffmpeg $(FFMPEG_FLAGS) -i $(patsubst %.ogg,%-%04d.ppm,$@) $@
pcgsolver:
	curl http://www.cs.ubc.ca/~rbridson/fluidsimulation/pcgsolver.tar.gz | tar xz
.PHONY: clean smoke-grid-play fluid-grid-play util-test-run

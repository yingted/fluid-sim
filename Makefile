FFMPEG_FLAGS := -framerate 10
FFMPEG_CONVERT := -c:v libx264 -qp 0 -y
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
	rm -f {smoke,fluid}-grid{,-gl}{,-*.ppm,.mkv,.trace} util-test
smoke-grid-0001.ppm: smoke-grid
	./smoke-grid
smoke-grid-play: smoke-grid-0001.ppm
	ffplay $(FFMPEG_FLAGS) -vf scale=-1:500 smoke-grid-%04d.ppm
#fluid-grid-0001.ppm: fluid-grid
#	set -o pipefail; ./fluid-grid | ./rpc.py
fluid-grid.trace fluid-grid-0001.ppm: fluid-grid rpc.py apitrace
	set -o pipefail; ./fluid-grid | apitrace/build/apitrace trace -o fluid-grid.trace ./rpc.py
fluid-grid-play: fluid-grid-0001.ppm
	ffplay $(FFMPEG_FLAGS) -vf scale=-1:500 fluid-grid-%04d.ppm
fluid-grid-gl.mkv: fluid-grid.trace apitrace
	apitrace/build/glretrace -s - $< | ffmpeg -f image2pipe -vcodec ppm -i pipe: $(FFMPEG_CONVERT) $@
%.mkv: %-0001.ppm
	ffmpeg $(FFMPEG_FLAGS) -i $(patsubst %.mkv,%-%04d.ppm,$@) $(FFMPEG_CONVERT) $@
apitrace:
	git clone git://github.com/apitrace/apitrace.git
	cd apitrace && cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release && make -C build
pcgsolver:
	curl http://www.cs.ubc.ca/~rbridson/fluidsimulation/pcgsolver.tar.gz | tar xz
power-test: power-test.cpp
	g++ $(CFLAGS) -o power-test{,.cpp} -lCGAL -lgmp -lboost_random -O3
power-test-out: power-test.cpp
	g++ -std=c++11 -Ipcgsolver -o power-test{-out,.cpp} -lCGAL -lboost_random -lgmp -D'OUTPUT="out.obj"' -DEXACT -lmpfr -lboost_thread-mt -O3
.PHONY: clean smoke-grid-play fluid-grid-play util-test-run

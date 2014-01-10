FFMPEG_FLAGS := -framerate 10
util-test: util-test.cpp util.hpp
	g++ -std=c++11 -o util-test{,.cpp}
smoke-grid: smoke-grid.cpp util.hpp
	g++ -std=c++11 -o smoke-grid{,.cpp}
clean:
	rm -f smoke-grid{,-*.ppm,.ogg} util-test
smoke-grid-0000.ppm: smoke-grid
	./smoke-grid
smoke-grid-play: smoke-grid
	ffplay $(FFMPEG_FLAGS) -vf scale=-1:500 smoke-grid-%04d.ppm
%.ogg: %-0000.ppm
	ffmpeg $(FFMPEG_FLAGS) -i smoke-grid-%04d.ppm $@
.PHONY: clean smoke-grid-play

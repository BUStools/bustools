RELEASE_OS ?= local
RELEASE_VERSION ?= local

.PHONY : build_linux build_mac build_windows compile_release clean

build_linux:
	mkdir -p build
	cmake -S . -B build
	# This will fail.
	- make -C build -j
	g++ -std=c++11 -O3 -DNDEBUG -static-libgcc -static-libstdc++ -rdynamic \
	build/src/CMakeFiles/bustools.dir/bustools_main.cpp.o \
	-o build/src/bustools build/src/libbustools_core.a -lpthread

compile_release:
	mkdir -p release/bustools
	cp -rf build/src/bustools release/bustools/
	cp -rf LICENSE release/bustools/
	cp -rf README.md release/bustools/
	tar -czvf release/bustools_${RELEASE_OS}-${RELEASE_VERSION}.tar.gz -C release bustools

clean:
	rm -rf build
	rm -rf cache
	rm -rf CMakeFiles
	rm -rf release

RELEASE_OS ?= local
RELEASE_VERSION ?= local

.PHONY : build compile_release_linux compile_release_mac compile_release_windows clean

build:
	mkdir -p build
	cd build \
	&& cmake .. -DLINK=static \
	&& make -j

compile_release_linux compile_release_mac:
	mkdir -p release/bustools
	cp -rf build/src/bustools release/bustools/
	cp -rf LICENSE release/bustools/
	cp -rf README.md release/bustools/
	cd release \
	&& tar -czvf bustools_${RELEASE_OS}-${RELEASE_VERSION}.tar.gz bustools

compile_release_windows:
	mkdir -p release/bustools
	cp -rf build/src/bustools.exe release/bustools/
	cp -rf LICENSE release/bustools/
	cp -rf README.md release/bustools/
	cd release \
	&& zip -r bustools_${RELEASE_OS}-${RELEASE_VERSION}.zip bustools

clean:
	rm -rf build
	rm -rf cache
	rm -rf CMakeFiles
	rm -rf release

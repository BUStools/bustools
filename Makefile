.PHONY : build_linux build_macos build_windows

build_linux:
	mkdir build
	cd build
	cmake ..
	make -j

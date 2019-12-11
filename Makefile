.PHONY : build_linux build_macos build_windows

build_linux:
	cmake -S . -B build
	make -C build

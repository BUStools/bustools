#!/bin/bash

# Run this on a Mac

git clone https://github.com/BUStools/bustools
cd bustools && mkdir build && cd build
cmake ..
make
mv src/bustools ../
cd ../../
tar --no-xattrs --exclude='._*' -czvf bustools_mac-master.tar.gz bustools/bustools bustools/README.md bustools/LICENSE



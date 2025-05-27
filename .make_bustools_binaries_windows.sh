#!/bin/bash

# Setup docker

# sudo apt install docker.io  # version 20.10.25-0ubuntu1~20.04.1
# sudo groupadd docker
# sudo usermod -aG docker ${USER}

# Exit and log back in


docker run --rm dockcross/windows-static-x64:20221217-6afd127 > ./dockcross-windows-static-x64
chmod +x ./dockcross-windows-static-x64
./dockcross-windows-static-x64 bash -c "rm -rf zlib-1.3.1.tar.gz && wget http://www.zlib.net/zlib-1.3.1.tar.gz && tar -xvzf zlib-1.3.1.tar.gz && cd zlib-1.3.1 && ./configure --static && make && cd .. && rm -rf bustools && git clone https://github.com/BUStools/bustools && cd bustools && mkdir build && cd build && cmake ..  -DZLIB_LIBRARY=/work/zlib-1.3.1/libz.a -DZLIB_INCLUDE_DIR=/work/zlib-1.3.1/ && make && mv src/bustools.exe ../ && cd ../../"

rm -rf bustools_windows-master.zip
zip bustools_windows-master.zip bustools/bustools.exe bustools/README.md bustools/LICENSE


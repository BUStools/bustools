#!/bin/bash

# Setup docker

# sudo apt install docker.io  # version 20.10.25-0ubuntu1~20.04.1
# sudo groupadd docker
# sudo usermod -aG docker ${USER}

# Exit and log back in

uid=$(id -u $USER|tr -d "\n" )
gid=$(id -g $USER|tr -d "\n" )
docker run --rm dockbuild/centos7:latest > ./dockbuild
chmod +x dockbuild
./dockbuild
docker run --rm dockbuild/centos7-devtoolset7-gcc7:latest > dockbuild-centos7-devtoolset7-gcc7-latest
docker run -ti -v ./:/work -e BUILDER_UID=$uid -e BUILDER_GID=$gid -e BUILDER_USER=$USER -e BUILDER_GROUP=$USER --platform linux dockbuild/centos7-devtoolset7-gcc7:latest \
	bash -c "rm -rf bustools && git clone https://github.com/BUStools/bustools && cd bustools && mkdir build && cd build && cmake .. && make && mv src/bustools ../ && cd ../../"

rm -rf bustools_linux-master.tar.gz
tar --no-xattrs --exclude='._*' -czvf bustools_linux-master.tar.gz bustools/bustools bustools/README.md bustools/LICENSE



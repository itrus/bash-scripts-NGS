#!/bin/bash
wget https://github.com/andersen-lab/ivar/archive/refs/heads/master.zip
unzip ./master.zip
cd ./ivar-master/
./autogen.sh
./configure
make
sudo make install
ivar version

#!/bin/bash

git clone https://github.com/ssloy/ultimaille.git

mkdir build
cd build 
cmake ..
make -j4



#!/bin/bash

# Install dependencies
yum install -y gcc-gfortran cmake autoconf automake libtool

ls
cd /io
ls
./bootstrap
./configure --prefix=$PWD

make install

mkdir -p /io/build
cp -R bin /io/build/
cp -R share /io/build/

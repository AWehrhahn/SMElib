#!/bin/bash

# Install dependencies
yum install -y gcc-gfortran cmake autoconf automake libtool

cd /io
./bootstrap
./configure --prefix=$PWD

make install

for PYBIN in /opt/python/*/bin; do
    echo "${PYBIN}"
    "${PYBIN}/pip" install -r /io/test/requirements.txt
    "${PYBIN}/pytest"
done

mkdir -p /io/build
cp -R lib /io/build/
cp -R share /io/build/
